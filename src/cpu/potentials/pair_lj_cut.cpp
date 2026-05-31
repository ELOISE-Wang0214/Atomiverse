/// @file cpu/potentials/pair_lj_cut.cpp
/// @brief CPU PairLJCut — LJ 12-6 shifted, xsimd ket-loop, OMP force reduction.

#include "atomiverse/cpu/potentials/pair_lj_cut.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/md/mpi/domain.h"
#include "atomiverse/cpu/md/pair_list.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

// xsimd vector batch<double> for the LJ inner kernel (Track B.2).
// Linked directly to atomiverse via target_link_libraries; -march=native is
// already gated to CXX in the top-level CMakeLists. We use the *default*
// architecture (NEON64 on Apple silicon, AVX2/AVX-512 on x86-64); on
// architectures without xsimd support the batch_t::size collapses to 1 and
// the loop degenerates into the scalar code path with no measurable
// overhead.
#include <xsimd/xsimd.hpp>

#include "atomiverse/cpu/box_simd.h"

namespace ATV {
namespace cpu {

PairLJCut::PairLJCut(double epsilon, double sigma, double cutoff)
    : Pair(static_cast<Real>(cutoff)), epsilon_(epsilon), sigma_(sigma) {
    sigma6_ = sigma * sigma * sigma * sigma * sigma * sigma;
    cutoff_sq_ = cutoff * cutoff;

    // Compute V_shift = V(r_cutoff). The base's cutoff_ holds the
    // canonical Å value; the derived caches (cutoff_sq_, v_shift_,
    // sigma6_) are computed once from it.
    double r2inv = double(1) / cutoff_sq_;
    double r6inv = r2inv * r2inv * r2inv * sigma6_;
    v_shift_ = double(4) * epsilon_ * (r6inv * r6inv - r6inv);
}

PairLJCut& PairLJCut::set_special_bonds(double s12, double s13, double s14) & {
    s12_ = s12;
    s13_ = s13;
    s14_ = s14;
    return *this;
}

PairLJCut PairLJCut::set_special_bonds(double s12, double s13, double s14) && {
    s12_ = s12;
    s13_ = s13;
    s14_ = s14;
    return std::move(*this);
}

void PairLJCut::build(const BuildContext& ctx) {
    n_threads_ = 1;
#ifdef _OPENMP
    n_threads_ = omp_get_max_threads();
#endif
    int n = ctx.frame.atoms.size();
    fx_thr_.resize(n_threads_);
    fy_thr_.resize(n_threads_);
    fz_thr_.resize(n_threads_);
    for (int t = 0; t < n_threads_; ++t) {
        fx_thr_[t].resize(n, double(0));
        fy_thr_[t].resize(n, double(0));
        fz_thr_[t].resize(n, double(0));
    }

    // Tag-indexed special-bonds exclusion: active when any special-bonds
    // tier is fully excluded (scale 0). The map itself is (re)built from
    // the live per-atom topology by `refresh_exclusion_map`, called on
    // every neighbor-list rebuild via `rebuild_filtered_list` — NOT cached
    // once here. A build-time-only map goes stale the moment multi-rank
    // `AtomMigration` moves an atom into a rank whose snapshot never held
    // that tag, dropping its exclusions and feeding a ~1.5 Å bonded pair
    // to the LJ kernel (MSU-PE multi-rank NPT collapse). We still build it
    // once here so a no-migration / single-`compute()` caller (e.g.
    // `Evaluate`) has a valid map before the first filtered walk.
    any_excl_ = (s12_ == 0.0) || (s13_ == 0.0) || (s14_ == 0.0);
    excluded_by_tag_.clear();
    if (!any_excl_)
        return;
    const auto* dd = ctx.runtime.dd;
    const int n_owners = (dd && ctx.runtime.mpi_size > 1)
                             ? dd->n_local
                             : ctx.frame.atoms.nlocal;
    refresh_exclusion_map(ctx.frame, n_owners);
}

void PairLJCut::refresh_exclusion_map(const Frame& f, int n_owners) {
    // Rebuild the tag-keyed exclusion lists from the CURRENT per-atom
    // `nspecial` / `special` topology of owner rows `[0, n_owners)`.
    // `special(i, ·)` holds GLOBAL TAGS, populated by `materialize` /
    // the LAMMPS-data reader via `Atom::build_special_list`, and travels
    // with the atom across migration (pack/unpack/move_row). Iterating
    // owners and recording each owner's own partners makes
    // `excluded_by_tag_[tag_i]` complete for every owner `i` the filtered
    // walk consults — correct under any partition.
    excluded_by_tag_.clear();
    const bool has_excl_12 = (s12_ == 0.0);
    const bool has_excl_13 = (s13_ == 0.0);
    const bool has_excl_14 = (s14_ == 0.0);
    if (!(has_excl_12 || has_excl_13 || has_excl_14))
        return;

    auto add_pair = [&](std::int64_t a, std::int64_t b) {
        if (a == b)
            return;
        excluded_by_tag_[a].push_back(b);
        excluded_by_tag_[b].push_back(a);
    };

    const auto& at = f.atoms;
    for (int i = 0; i < n_owners; ++i) {
        const std::int64_t ti = at.tag[i];
        const int n12 = at.nspecial(i, 0);
        const int n13 = at.nspecial(i, 1);
        const int n14 = at.nspecial(i, 2);
        if (has_excl_12) {
            for (int s = 0; s < n12; ++s)
                add_pair(ti, at.special(i, s));
        }
        if (has_excl_13) {
            for (int s = n12; s < n13; ++s)
                add_pair(ti, at.special(i, s));
        }
        if (has_excl_14) {
            for (int s = n13; s < n14; ++s)
                add_pair(ti, at.special(i, s));
        }
    }

    for (auto& kv : excluded_by_tag_) {
        auto& list = kv.second;
        std::sort(list.begin(), list.end());
        list.erase(std::unique(list.begin(), list.end()), list.end());
    }
}

void PairLJCut::rebuild_filtered_list(const Frame& f,
                                      const md::PairList& pairs) {
    // Walk the upstream CSR once, drop excluded pairs, write a dense CSR
    // into `filtered_*`. Runs once per NL rebuild (~every `every` steps)
    // amortizing the `tag_ptr[j]` gather + `binary_search` that the old
    // in-kernel path paid on every launch × every pair.
    const int n = pairs.inum;
    // Refresh the exclusion map from the live topology of the current
    // owners before filtering. `pairs.inum` is the authoritative owned
    // count (trimmed to `dd->n_local` in multi-rank). This is what keeps
    // special-bonds exclusions correct after an `AtomMigration` reshuffles
    // ownership — the stale build-time map was the MSU-PE collapse cause.
    refresh_exclusion_map(f, n);
    filtered_inum_ = n;
    filtered_numneigh_.assign(n, 0);
    filtered_firstneigh_.assign(n, 0);
    filtered_neighbors_.clear();
    filtered_neighbors_.reserve(pairs.total_pairs);

    const int* numneigh = pairs.numneigh.data();
    const int* firstneigh = pairs.firstneigh.data();
    const int* neighbors = pairs.neighbors.data();
    const std::int64_t* tag = f.atoms.tag.data();

    for (int i = 0; i < n; ++i) {
        const std::int64_t i_tag = tag[i];
        auto it = excluded_by_tag_.find(i_tag);
        const std::vector<std::int64_t>* my_excl =
            (it != excluded_by_tag_.end()) ? &it->second : nullptr;

        const int off = firstneigh[i];
        const int jn = numneigh[i];
        filtered_firstneigh_[i] = static_cast<int>(filtered_neighbors_.size());
        int kept = 0;
        for (int jj = 0; jj < jn; ++jj) {
            const int j = neighbors[off + jj];
            if (my_excl) {
                const std::int64_t j_tag = tag[j];
                if (std::binary_search(my_excl->begin(), my_excl->end(), j_tag))
                    continue;
            }
            filtered_neighbors_.push_back(j);
            ++kept;
        }
        filtered_numneigh_[i] = kept;
    }
}

void PairLJCut::compute(const LaunchContext& ctx) {
    auto& f = ctx.frame;
    if (!ctx.runtime.pair_list) {
        throw std::runtime_error(
            "PairLJCut::compute: runtime.pair_list is null. Did you add a "
            "NeighborList component to the force phase ahead of LJ?");
    }
    const auto& pairs = *ctx.runtime.pair_list;
    const auto& box = f.box;
    int n_force = f.atoms.size();
    const auto* dd = ctx.runtime.dd;
    if (dd && ctx.runtime.mpi_size > 1) {
        n_force = dd->n_local + dd->n_ghost;
    }
    double energy_total = 0;

    // Exclusion-filtered CSR (built once per NL rebuild). See
    // `rebuild_filtered_list` for the invariant. When no exclusions are
    // configured we iterate the upstream CSR directly — no copy, no
    // allocation, no per-pair branch.
    const bool use_filtered = any_excl_;
    if (use_filtered && pairs.version != filtered_version_) {
        rebuild_filtered_list(f, pairs);
        filtered_version_ = pairs.version;
    }
    const int k_inum = use_filtered ? filtered_inum_ : pairs.inum;
    const int* __restrict k_ilist = use_filtered ? nullptr // contiguous 0..n-1
                                                 : pairs.ilist.data();
    const int* __restrict k_numneigh =
        use_filtered ? filtered_numneigh_.data() : pairs.numneigh.data();
    const int* __restrict k_firstneigh =
        use_filtered ? filtered_firstneigh_.data() : pairs.firstneigh.data();
    const int* __restrict k_neighbors =
        use_filtered ? filtered_neighbors_.data() : pairs.neighbors.data();

    // Pair virial accumulators (Voigt order: xx yy zz xy xz yz). Reduced
    // across OpenMP threads via OMP reduction; cross-rank reduction is
    // handled by ReduceEnergy at end_of_step.
    double v_xx = 0, v_yy = 0, v_zz = 0, v_xy = 0, v_xz = 0, v_yz = 0;

#pragma omp parallel reduction(+ : energy_total) reduction(+ : v_xx)           \
    reduction(+ : v_yy) reduction(+ : v_zz) reduction(+ : v_xy)                \
    reduction(+ : v_xz) reduction(+ : v_yz)
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        auto& fxt = fx_thr_[tid];
        auto& fyt = fy_thr_[tid];
        auto& fzt = fz_thr_[tid];
        std::fill(fxt.begin(), fxt.begin() + n_force, double(0));
        std::fill(fyt.begin(), fyt.begin() + n_force, double(0));
        std::fill(fzt.begin(), fzt.begin() + n_force, double(0));

        // Hoist box edges + LJ constants to locals so the inner loop keeps
        // them in registers instead of reloading the Mat3 / member fields
        // on every pair.
        const double lx = box(0, 0), ly = box(1, 1), lz = box(2, 2);
        const double half_lx = double(0.5) * lx;
        const double half_ly = double(0.5) * ly;
        const double half_lz = double(0.5) * lz;
        const double rc_sq = cutoff_sq_;
        const double sig6 = sigma6_;
        const double eps4 = double(4) * epsilon_;
        const double eps24 = double(24) * epsilon_;
        const double vshift = v_shift_;

        // Raw pointers to the per-thread force arenas; indexing via
        // `operator[]` on xtensor goes through a strided stride check on
        // every access even for contiguous 1-D storage. `.data()` drops to
        // a plain `double*` the compiler keeps in a register.
        const double* __restrict x_ptr = f.atoms.x.data();
        const double* __restrict y_ptr = f.atoms.y.data();
        const double* __restrict z_ptr = f.atoms.z.data();
        double* __restrict fxt_ptr = fxt.data();
        double* __restrict fyt_ptr = fyt.data();
        double* __restrict fzt_ptr = fzt.data();

        // Track B.2 — xsimd ket-loop vectorization of the inner pair loop.
        //
        // Strategy: pack W=batch_t::size neighbors into vector registers,
        // run all distance / cutoff / polynomial / force math through
        // SIMD (NEON64 W=2 on Apple silicon, AVX2 W=4, AVX-512 W=8).
        // Excluded / out-of-cutoff lanes are masked to zero contributions
        // via xsimd::select — branch-free.
        //
        // Bit-identicality: per-lane SIMD ops on NEON / AVX use plain
        // (non-FMA) multiplies and adds (xsimd 11.2.0 routes operator*/+
        // to vmulq/vaddq, not vfmaq), and division uses vdivq_f64 (full
        // IEEE divide). The horizontal accumulation into per-atom force
        // (`fi_x`), per-thread force buffers (`fxt_ptr[j]`), energy and
        // virial is performed by extracting lanes one at a time and
        // adding them in the *same* linear order as the original scalar
        // loop. This preserves the exact summation tree, so per-pair
        // contributions are bit-identical iff the per-lane intermediates
        // (`r2`, `fpair`, `fxp`, etc.) are bit-identical to scalar.
        //
        // The one residual ULP risk is `-ffp-contract=on` (clang default):
        // the *scalar* baseline may emit a fused multiply-add for
        // `dx*dx + dy*dy + dz*dz` written in a single statement, while
        // xsimd's batch operator+ stays unfused. We measure and report
        // the actual delta; if it exceeds 1e-12 the spec acceptance gate
        // (ac-012) will need a tolerance refinement (caller's call).
        using batch_t = xsimd::batch<double>;
        constexpr std::size_t W = batch_t::size;

        const batch_t v_lx(lx), v_ly(ly), v_lz(lz);
        const batch_t v_half_lx(half_lx), v_half_ly(half_ly),
            v_half_lz(half_lz);
        const batch_t v_neg_half_lx(-half_lx), v_neg_half_ly(-half_ly),
            v_neg_half_lz(-half_lz);
        // Floor-form min-image constants. The legacy single-`select` pair
        // form below only handled |Δ| ≤ L; cross-PBC ghosts at ±2L are
        // silently dropped by the cutoff filter. The floor form
        // `Δ -= L · round(Δ/L)` is correct for arbitrary |Δ| in one
        // SIMD op without a branch (see Box::min_image fix).
        const batch_t v_inv_lx(1.0 / lx), v_inv_ly(1.0 / ly),
            v_inv_lz(1.0 / lz);
        const batch_t v_rc_sq(rc_sq);
        const batch_t v_sig6(sig6);
        const batch_t v_eps4(eps4), v_eps24(eps24), v_vshift(vshift);
        const batch_t v_one(1.0), v_two(2.0), v_zero(0.0);

        // Stack-resident packing buffers reused per ket iteration.
        // alignas(64) keeps store_aligned legal on AVX-512 boxes; on NEON
        // 16-byte alignment suffices but 64 is cheap and harmless.
        alignas(64) double xj_buf[16], yj_buf[16], zj_buf[16];
        alignas(64) double dx_buf[16], dy_buf[16], dz_buf[16];
        alignas(64) double fxp_buf[16], fyp_buf[16], fzp_buf[16];
        alignas(64) double r2_buf[16], r6inv_buf[16];
        alignas(64) std::int32_t mask_buf[16];

#pragma omp for schedule(static)
        for (int ii = 0; ii < k_inum; ++ii) {
            // Filtered path: i == ii (neighbors are already pruned of
            // excluded partners). Upstream path: `k_ilist[ii]` gives the
            // owner atom index when NL builds a sparse subset.
            const int i = use_filtered ? ii : k_ilist[ii];
            const int jnum = k_numneigh[i];
            const int off = k_firstneigh[i];
            // Atom-i state hoisted into registers for the inner loop.
            // `fi_*` accumulates Σ_j F_ij once; we store back to `fxt_ptr[i]`
            // a single time at the end -- removes an `fxt[i] += ...` memory
            // RMW per pair (~50 ns bubble per iter on Apple Silicon).
            const double xi = x_ptr[i];
            const double yi = y_ptr[i];
            const double zi = z_ptr[i];
            const batch_t v_xi(xi), v_yi(yi), v_zi(zi);
            double fi_x = 0, fi_y = 0, fi_z = 0;

            int jj = 0;
            // ── SIMD ket loop: W neighbors per iteration ────────────
            for (; jj + static_cast<int>(W) <= jnum;
                 jj += static_cast<int>(W)) {
                // Manual gather: pack the next W neighbor coords into the
                // aligned stack buffers. j-indices are non-contiguous (NL
                // is a permutation), so a true SIMD gather instruction
                // would be needed for 1-cycle load -- on NEON64 there is
                // no double-precision gather, so the scalar-load loop is
                // the right thing here.
                for (std::size_t lane = 0; lane < W; ++lane) {
                    const int j =
                        k_neighbors[off + jj + static_cast<int>(lane)];
                    xj_buf[lane] = x_ptr[j];
                    yj_buf[lane] = y_ptr[j];
                    zj_buf[lane] = z_ptr[j];
                }
                batch_t v_xj = batch_t::load_aligned(xj_buf);
                batch_t v_yj = batch_t::load_aligned(yj_buf);
                batch_t v_zj = batch_t::load_aligned(zj_buf);

                batch_t v_dx = v_xi - v_xj;
                batch_t v_dy = v_yi - v_yj;
                batch_t v_dz = v_zi - v_zj;
                // Branch-free minimum image: select(dx > +L/2, dx-L, dx);
                // then select(dx < -L/2, dx+L, dx). Hoisted constants
                // (v_lx, v_half_lx, ...) keep the compiler from reloading
                // box members across the inner loop -- swapping this for
                // box.min_image_batch costs ~7% per-step (member loads
                // through `this` are not hoisted under aliasing analysis).
                // Floor-form min-image (single shot, handles |Δ| up to
                // any integer multiple of L). The legacy single-`select`
                // pair form was correct for `|Δ| ≤ L` but silently left
                // `|Δ| = L` for cross-PBC ghosts at ±2L (then the
                // cutoff filter dropped legitimate pairs).
                v_dx = v_dx - v_lx * xsimd::round(v_dx * v_inv_lx);
                v_dy = v_dy - v_ly * xsimd::round(v_dy * v_inv_ly);
                v_dz = v_dz - v_lz * xsimd::round(v_dz * v_inv_lz);
                // Suppress unused-variable warnings for the legacy
                // half-L locals which are now reserved for the v3
                // hot-path microbench (`bench_md_step` v_half_*
                // sensitivity).
                (void)v_half_lx;
                (void)v_neg_half_lx;
                (void)v_half_ly;
                (void)v_neg_half_ly;
                (void)v_half_lz;
                (void)v_neg_half_lz;

                // r2 with parenthesized addition tree to keep order
                // matching scalar `(dx*dx + dy*dy) + dz*dz`.
                batch_t v_dxsq = v_dx * v_dx;
                batch_t v_dysq = v_dy * v_dy;
                batch_t v_dzsq = v_dz * v_dz;
                batch_t v_r2 = (v_dxsq + v_dysq) + v_dzsq;

                // Cutoff mask. Lanes outside the cutoff get fpair=0 so
                // they contribute nothing to forces / energy / virial; the
                // scalar code reaches those lanes via its `continue`.
                auto mask = v_r2 < v_rc_sq;
                // Guard against r2 == 0 (would-be NaN from 1/r2). All
                // lanes inside the cutoff are guaranteed > 0 (no atom
                // self-pair in the half list); padding lanes outside
                // cutoff are masked anyway. We replace masked-out r2 with
                // 1.0 to keep the divide finite.
                batch_t v_r2_safe = xsimd::select(mask, v_r2, v_one);

                batch_t v_r2inv = v_one / v_r2_safe;
                batch_t v_r6inv = ((v_r2inv * v_r2inv) * v_r2inv) * v_sig6;
                batch_t v_fpair =
                    v_eps24 * v_r2inv * (v_two * v_r6inv * v_r6inv - v_r6inv);
                v_fpair = xsimd::select(mask, v_fpair, v_zero);

                batch_t v_fxp = v_fpair * v_dx;
                batch_t v_fyp = v_fpair * v_dy;
                batch_t v_fzp = v_fpair * v_dz;

                v_dx.store_aligned(dx_buf);
                v_dy.store_aligned(dy_buf);
                v_dz.store_aligned(dz_buf);
                v_fxp.store_aligned(fxp_buf);
                v_fyp.store_aligned(fyp_buf);
                v_fzp.store_aligned(fzp_buf);
                v_r6inv.store_aligned(r6inv_buf);
                // Mask: store as integer (1 inside cutoff, 0 outside) so
                // the scalar accumulation loop below can branch on it.
                // batch_bool stores as the underlying SIMD compare result;
                // we materialize via a select-to-int.
                batch_t v_mask_d = xsimd::select(mask, v_one, v_zero);
                alignas(64) double mask_d_buf[16];
                v_mask_d.store_aligned(mask_d_buf);

                // ── Scalar accumulation in original linear order ──
                // Lane-by-lane add into fi_*, fxt[j], v_xx..v_yz, energy.
                // This is identical FP order to the scalar loop: each
                // surviving pair contributes the same per-pair triple
                // (fxp, fyp, fzp) the scalar code would have computed.
                for (std::size_t lane = 0; lane < W; ++lane) {
                    if (mask_d_buf[lane] == 0.0)
                        continue;
                    const int j =
                        k_neighbors[off + jj + static_cast<int>(lane)];
                    const double dx_l = dx_buf[lane];
                    const double dy_l = dy_buf[lane];
                    const double dz_l = dz_buf[lane];
                    const double fxp_l = fxp_buf[lane];
                    const double fyp_l = fyp_buf[lane];
                    const double fzp_l = fzp_buf[lane];
                    const double r6inv_l = r6inv_buf[lane];

                    fi_x += fxp_l;
                    fi_y += fyp_l;
                    fi_z += fzp_l;
                    fxt_ptr[j] -= fxp_l;
                    fyt_ptr[j] -= fyp_l;
                    fzt_ptr[j] -= fzp_l;

                    v_xx += dx_l * fxp_l;
                    v_yy += dy_l * fyp_l;
                    v_zz += dz_l * fzp_l;
                    v_xy += dx_l * fyp_l;
                    v_xz += dx_l * fzp_l;
                    v_yz += dy_l * fzp_l;

                    energy_total +=
                        eps4 * (r6inv_l * r6inv_l - r6inv_l) - vshift;
                }
                (void)mask_buf;
            }

            // ── Scalar tail (last jnum % W neighbors) ──────────────
            for (; jj < jnum; ++jj) {
                const int j = k_neighbors[off + jj];
                double dx = xi - x_ptr[j];
                double dy = yi - y_ptr[j];
                double dz = zi - z_ptr[j];
                // Inlined branchy min-image: matches the SIMD ket-loop form
                // above and keeps the hoisted half_lx / lx locals in
                // registers across the entire scalar tail. (Substituting
                // box.min_image here costs the same ~7% as the SIMD path.)
                // Multi-wrap to handle cross-PBC ghosts at ±2L on small
                // Cartesian grids (see Box::min_image fix).
                while (dx > half_lx)
                    dx -= lx;
                while (dx < -half_lx)
                    dx += lx;
                while (dy > half_ly)
                    dy -= ly;
                while (dy < -half_ly)
                    dy += ly;
                while (dz > half_lz)
                    dz -= lz;
                while (dz < -half_lz)
                    dz += lz;
                const double r2 = dx * dx + dy * dy + dz * dz;
                if (r2 >= rc_sq)
                    continue;

                const double r2inv = double(1) / r2;
                const double r6inv = r2inv * r2inv * r2inv * sig6;
                const double fpair =
                    eps24 * r2inv * (double(2) * r6inv * r6inv - r6inv);

                const double fxp = fpair * dx;
                const double fyp = fpair * dy;
                const double fzp = fpair * dz;
                fi_x += fxp;
                fi_y += fyp;
                fi_z += fzp;
                fxt_ptr[j] -= fxp;
                fyt_ptr[j] -= fyp;
                fzt_ptr[j] -= fzp;

                v_xx += dx * fxp;
                v_yy += dy * fyp;
                v_zz += dz * fzp;
                v_xy += dx * fyp;
                v_xz += dx * fzp;
                v_yz += dy * fzp;

                energy_total += eps4 * (r6inv * r6inv - r6inv) - vshift;
            }
            fxt_ptr[i] += fi_x;
            fyt_ptr[i] += fi_y;
            fzt_ptr[i] += fi_z;
        }

// Reduce per-thread forces into frame
#pragma omp for schedule(static)
        for (int i = 0; i < n_force; ++i) {
            for (int t = 0; t < n_threads_; ++t) {
                f.atoms.fx[i] += fx_thr_[t][i];
                f.atoms.fy[i] += fy_thr_[t][i];
                f.atoms.fz[i] += fz_thr_[t][i];
            }
        }
    }

    f.scalars.energy += energy_total;
    f.scalars.energy_pair += energy_total;
    f.scalars.stress[0] += v_xx;
    f.scalars.stress[1] += v_yy;
    f.scalars.stress[2] += v_zz;
    f.scalars.stress[3] += v_xy;
    f.scalars.stress[4] += v_xz;
    f.scalars.stress[5] += v_yz;
}

void PairLJCut::destroy(const DestroyContext& ctx) {
    fx_thr_.clear();
    fy_thr_.clear();
    fz_thr_.clear();
    excluded_by_tag_.clear();
    filtered_neighbors_.clear();
    filtered_numneigh_.clear();
    filtered_firstneigh_.clear();
    filtered_inum_ = 0;
    filtered_version_ = 0;
    any_excl_ = false;
}

} // namespace cpu
} // namespace ATV
