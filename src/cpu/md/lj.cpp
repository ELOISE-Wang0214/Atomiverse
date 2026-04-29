/// @file cpu/md/lj.cpp
/// @brief CPU LJ 12-6 (shifted) with OpenMP per-thread force reduction.

#include "atomiverse/cpu/md/lj.h"
#include "atomiverse/cpu/md/mpi/domain.h"
#include "atomiverse/cpu/md/pair_list.h"
#include "atomiverse/core/mat3.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ATV {
namespace cpu {
namespace md {

static inline void minimum_image(double& dx, double& dy, double& dz, const Mat3& box) {
    double lx = box(0, 0), ly = box(1, 1), lz = box(2, 2);
    if (dx >  lx * double(0.5)) dx -= lx;
    if (dx < -lx * double(0.5)) dx += lx;
    if (dy >  ly * double(0.5)) dy -= ly;
    if (dy < -ly * double(0.5)) dy += ly;
    if (dz >  lz * double(0.5)) dz -= lz;
    if (dz < -lz * double(0.5)) dz += lz;
}

LJ::LJ(double epsilon, double sigma, double cutoff)
    : epsilon_(epsilon), sigma_(sigma), cutoff_(cutoff) {
    sigma6_ = sigma * sigma * sigma * sigma * sigma * sigma;
    cutoff_sq_ = cutoff * cutoff;

    // Compute V_shift = V(r_cutoff)
    double r2inv = double(1) / cutoff_sq_;
    double r6inv = r2inv * r2inv * r2inv * sigma6_;
    v_shift_ = double(4) * epsilon_ * (r6inv * r6inv - r6inv);
}

LJ& LJ::set_special_bonds(double s12, double s13, double s14) & {
    s12_ = s12;
    s13_ = s13;
    s14_ = s14;
    return *this;
}

LJ LJ::set_special_bonds(double s12, double s13, double s14) && {
    s12_ = s12;
    s13_ = s13;
    s14_ = s14;
    return std::move(*this);
}

void LJ::build(const cpu::BuildContext& ctx) {
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

    // Build tag-indexed exclusion map when any special-bonds scale is 0.
    // Frame topology stores GLOBAL TAGS (tag == local index in single-rank
    // by default; differs once partition compacts atoms in multi-rank).
    // Per-tag exclusion is correct under any partition: at launch time we
    // read the partner's tag from frame.atoms.tag[j] and look it up in
    // the owner's excluded set.
    const bool has_excl_12 = (s12_ == 0.0);
    const bool has_excl_13 = (s13_ == 0.0);
    const bool has_excl_14 = (s14_ == 0.0);
    const bool any_excl = has_excl_12 || has_excl_13 || has_excl_14;

    excluded_by_tag_.clear();
    if (!any_excl) return;

    auto add_pair = [&](std::int64_t a, std::int64_t b) {
        if (a == b) return;
        excluded_by_tag_[a].push_back(b);
        excluded_by_tag_[b].push_back(a);
    };

    const auto& f = ctx.frame;
    if (has_excl_12) {
        const int nb = f.bonds.size();
        for (int b = 0; b < nb; ++b) {
            add_pair(static_cast<std::int64_t>(f.bonds.atomi[b]),
                     static_cast<std::int64_t>(f.bonds.atomj[b]));
        }
    }
    if (has_excl_13) {
        const int na = f.angles.size();
        for (int a = 0; a < na; ++a) {
            add_pair(static_cast<std::int64_t>(f.angles.atomi[a]),
                     static_cast<std::int64_t>(f.angles.atomk[a]));
        }
    }
    if (has_excl_14) {
        const int nd = f.dihedrals.size();
        for (int d = 0; d < nd; ++d) {
            add_pair(static_cast<std::int64_t>(f.dihedrals.atomi[d]),
                     static_cast<std::int64_t>(f.dihedrals.atoml[d]));
        }
    }

    for (auto& kv : excluded_by_tag_) {
        auto& list = kv.second;
        std::sort(list.begin(), list.end());
        list.erase(std::unique(list.begin(), list.end()), list.end());
    }
}

void LJ::rebuild_filtered_list(const cpu::Frame& f, const PairList& pairs) {
    // Walk the upstream CSR once, drop excluded pairs, write a dense CSR
    // into `filtered_*`. Runs once per NL rebuild (~every `every` steps)
    // amortizing the `tag_ptr[j]` gather + `binary_search` that the old
    // in-kernel path paid on every launch × every pair.
    const int n = pairs.inum;
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
        const int jn  = numneigh[i];
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

void LJ::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    if (!ctx.runtime.pair_list) {
        throw std::runtime_error(
            "LJ::launch: runtime.pair_list is null. Did you add a "
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
    const bool use_filtered = !excluded_by_tag_.empty();
    if (use_filtered && pairs.version != filtered_version_) {
        rebuild_filtered_list(f, pairs);
        filtered_version_ = pairs.version;
    }
    const int  k_inum       = use_filtered ? filtered_inum_
                                           : pairs.inum;
    const int* __restrict k_ilist     = use_filtered
                                           ? nullptr  // contiguous 0..n-1
                                           : pairs.ilist.data();
    const int* __restrict k_numneigh  = use_filtered
                                           ? filtered_numneigh_.data()
                                           : pairs.numneigh.data();
    const int* __restrict k_firstneigh= use_filtered
                                           ? filtered_firstneigh_.data()
                                           : pairs.firstneigh.data();
    const int* __restrict k_neighbors = use_filtered
                                           ? filtered_neighbors_.data()
                                           : pairs.neighbors.data();

    // Pair virial accumulators (Voigt order: xx yy zz xy xz yz). Reduced
    // across OpenMP threads via OMP reduction; cross-rank reduction is
    // handled by ReduceEnergy at end_of_step.
    double v_xx = 0, v_yy = 0, v_zz = 0, v_xy = 0, v_xz = 0, v_yz = 0;

    #pragma omp parallel reduction(+:energy_total) \
        reduction(+:v_xx) reduction(+:v_yy) reduction(+:v_zz) \
        reduction(+:v_xy) reduction(+:v_xz) reduction(+:v_yz)
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
            double fi_x = 0, fi_y = 0, fi_z = 0;

            for (int jj = 0; jj < jnum; ++jj) {
                const int j = k_neighbors[off + jj];
                double dx = xi - x_ptr[j];
                double dy = yi - y_ptr[j];
                double dz = zi - z_ptr[j];
                if (dx >  half_lx) dx -= lx;
                if (dx < -half_lx) dx += lx;
                if (dy >  half_ly) dy -= ly;
                if (dy < -half_ly) dy += ly;
                if (dz >  half_lz) dz -= lz;
                if (dz < -half_lz) dz += lz;
                const double r2 = dx * dx + dy * dy + dz * dz;
                if (r2 >= rc_sq) continue;

                const double r2inv = double(1) / r2;
                const double r6inv = r2inv * r2inv * r2inv * sig6;
                const double fpair = eps24 * r2inv * (double(2) * r6inv * r6inv - r6inv);

                const double fxp = fpair * dx;
                const double fyp = fpair * dy;
                const double fzp = fpair * dz;
                fi_x += fxp;
                fi_y += fyp;
                fi_z += fzp;
                fxt_ptr[j] -= fxp;
                fyt_ptr[j] -= fyp;
                fzt_ptr[j] -= fzp;

                // Pair virial contribution: σ_αβ += r_α · F_β. Half-list
                // counts each pair once with the local owner as i, so this
                // captures the full pairwise virial when summed across
                // ranks (same convention as energy).
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

    f.energy += energy_total;
    f.stress[0] += v_xx;
    f.stress[1] += v_yy;
    f.stress[2] += v_zz;
    f.stress[3] += v_xy;
    f.stress[4] += v_xz;
    f.stress[5] += v_yz;
}

void LJ::destroy(const cpu::DestroyContext& ctx) {
    fx_thr_.clear();
    fy_thr_.clear();
    fz_thr_.clear();
    excluded_by_tag_.clear();
    filtered_neighbors_.clear();
    filtered_numneigh_.clear();
    filtered_firstneigh_.clear();
    filtered_inum_ = 0;
    filtered_version_ = 0;
}

} // namespace md
} // namespace cpu
} // namespace ATV
