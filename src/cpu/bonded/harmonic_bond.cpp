/// @file cpu/bonded/harmonic_bond.cpp
/// @brief CPU harmonic bond stretching implementation.

#include "atomiverse/cpu/bonded/harmonic_bond.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/md/mpi/domain.h"

#include <algorithm>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ATV {
namespace cpu {
namespace bonded {

static inline void minimum_image(double& dx, double& dy, double& dz, const Mat3& box) {
    double lx = box(0, 0), ly = box(1, 1), lz = box(2, 2);
    if (dx >  lx * double(0.5)) dx -= lx;
    if (dx < -lx * double(0.5)) dx += lx;
    if (dy >  ly * double(0.5)) dy -= ly;
    if (dy < -ly * double(0.5)) dy += ly;
    if (dz >  lz * double(0.5)) dz -= lz;
    if (dz < -lz * double(0.5)) dz += lz;
}

HarmonicBond::HarmonicBond(std::vector<double> K, std::vector<double> r0)
    : K_(std::move(K)), r0_(std::move(r0)) {}

void HarmonicBond::build(const cpu::BuildContext& ctx) {
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
}

void HarmonicBond::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    const int n_bonds = static_cast<int>(f.bonds.atomi.size());
    // Bond virial accumulators (Voigt xx yy zz xy xz yz). Bonds
    // contribute significantly to the stress tensor for stiff-bond
    // systems (e.g. PE chains, K ≈ 15 eV/Å²): without this term
    // the NPT barostat sees a much-too-high isotropic pressure and
    // over-expands the box monotonically. Matches LAMMPS
    // `bond_harmonic::compute` which accumulates the same
    // `fpair * Δ_α * Δ_β` terms into its virial array.
    double v_xx = 0, v_yy = 0, v_zz = 0, v_xy = 0, v_xz = 0, v_yz = 0;
    if (n_bonds == 0) return;

    const auto* dd = ctx.runtime.dd;
    const bool multi_rank = (dd != nullptr) && (ctx.runtime.mpi_size > 1);
    int n_atoms = f.atoms.size();
    if (multi_rank) {
        n_atoms = dd->n_local + dd->n_ghost;
        if (local_bonds_generation_ != dd->neighbor_generation) {
            local_bonds_.clear();
            local_bonds_.reserve(static_cast<std::size_t>(n_bonds));
            for (int b = 0; b < n_bonds; ++b) {
                const std::int64_t ti = static_cast<std::int64_t>(f.bonds.atomi[b]);
                const std::int64_t tj = static_cast<std::int64_t>(f.bonds.atomj[b]);
                const std::int64_t owner_tag = (ti < tj) ? ti : tj;
                auto it_owner = dd->tag_to_local.find(owner_tag);
                if (it_owner != dd->tag_to_local.end() &&
                    it_owner->second < dd->n_local) {
                    local_bonds_.push_back(b);
                }
            }
            local_bonds_generation_ = dd->neighbor_generation;
        }
    }
    const int n_work = multi_rank ? static_cast<int>(local_bonds_.size())
                                  : n_bonds;
    if (n_work == 0) return;

    double energy_total = 0;

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
        std::fill(fxt.begin(), fxt.begin() + n_atoms, double(0));
        std::fill(fyt.begin(), fyt.begin() + n_atoms, double(0));
        std::fill(fzt.begin(), fzt.begin() + n_atoms, double(0));

        // Pointer hoist: drop xtensor's strided indexing on every access.
        const double* __restrict x_ptr = f.atoms.x.data();
        const double* __restrict y_ptr = f.atoms.y.data();
        const double* __restrict z_ptr = f.atoms.z.data();

        #pragma omp for schedule(static)
        for (int ib = 0; ib < n_work; ++ib) {
            const int b = multi_rank ? local_bonds_[static_cast<std::size_t>(ib)]
                                     : ib;
            int i = f.bonds.atomi[b];
            int j = f.bonds.atomj[b];
            int type = f.bonds.type[b];

            // Multi-rank: bonds.atomi/atomj store global tags. Bond owner =
            // rank holding min(tag) atom (LAMMPS convention). Skip bonds we
            // don't own, and skip bonds whose other atom is not within
            // halo skin (its tag won't be in tag_to_local).
            if (multi_rank) {
                const std::int64_t ti = static_cast<std::int64_t>(i);
                const std::int64_t tj = static_cast<std::int64_t>(j);
                auto it_i = dd->tag_to_local.find(ti);
                auto it_j = dd->tag_to_local.find(tj);
                if (it_i == dd->tag_to_local.end() ||
                    it_j == dd->tag_to_local.end()) continue;
                i = it_i->second;
                j = it_j->second;
            }

            // Pack the bond tuple's positions into stack-resident triples.
            // Three contiguous loads per atom let the compiler keep them in
            // registers and removes the strided indirection through xtensor's
            // `operator[]` from the math below — speeds up auto-vectorization
            // of the per-bond arithmetic. Bit-identical to the inlined form
            // (same loads, same order of FP operations).
            const double ri[3] = { x_ptr[i], y_ptr[i], z_ptr[i] };
            const double rj[3] = { x_ptr[j], y_ptr[j], z_ptr[j] };
            double dx = ri[0] - rj[0];
            double dy = ri[1] - rj[1];
            double dz = ri[2] - rj[2];
            minimum_image(dx, dy, dz, f.box);

            double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            double dr = r - r0_[type];

            // V = K(r - r0)^2, F = -2K(r - r0) * r_hat
            energy_total += K_[type] * dr * dr;

            double fpair = -double(2) * K_[type] * dr / r;
            double fx_val = fpair * dx;
            double fy_val = fpair * dy;
            double fz_val = fpair * dz;

            fxt[i] += fx_val;  fxt[j] -= fx_val;
            fyt[i] += fy_val;  fyt[j] -= fy_val;
            fzt[i] += fz_val;  fzt[j] -= fz_val;

            // Virial: σ_αβ += Δ_α · F_β. Per-bond contribution,
            // counted once by the owner rank (matches the LAMMPS
            // convention of accumulating in the rank whose `atomi`
            // or `atomj` is locally owned — the bond-ownership
            // filter above guarantees each bond is counted once).
            v_xx += dx * fx_val;
            v_yy += dy * fy_val;
            v_zz += dz * fz_val;
            v_xy += dx * fy_val;
            v_xz += dx * fz_val;
            v_yz += dy * fz_val;
        }

        #pragma omp for schedule(static)
        for (int i = 0; i < n_atoms; ++i) {
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

void HarmonicBond::destroy(const cpu::DestroyContext& ctx) {
    fx_thr_.clear();
    fy_thr_.clear();
    fz_thr_.clear();
    local_bonds_.clear();
    local_bonds_generation_ = ~std::uint64_t{0};
}

} // namespace bonded
} // namespace cpu
} // namespace ATV
