/// @file cpu/pw/pw_solver.cpp
/// @brief PW implementation — orchestrates the internal step helpers
///        through the plane-wave SCF self-consistency cycle.
///
/// Low-level callers construct `PW` explicitly and execute its
/// build/run lifecycle; workflow callers use `PWEvaluator`.
///
/// Problem assembly (Ewald, G-vectors, FFT grids, V_loc, per-k KB
/// projectors / atomic wavefunctions on q = k+G, initial density) is the
/// shared `internal::setup_pw_problem()` — the same routine the CUDA
/// backend uses, so the two backends stay numerically aligned. This file
/// owns only the CPU SCF loop.
///
/// MPI: when built with `ATV_HAS_MPI` and MPI is initialized, k-points
/// are distributed round-robin across `MPI_COMM_WORLD`; the band-energy
/// scalar and the real-space density accumulation are Allreduced each
/// SCF iteration, so all ranks hold the identical mixed density and
/// take identical convergence decisions.

#include "atomiverse/cpu/planewave/pw_solver.h"

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/log.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/planewave/internal/convergence.h"
#include "atomiverse/cpu/planewave/internal/davidson.h"
#include "atomiverse/cpu/planewave/internal/density_builder.h"
#include "atomiverse/cpu/planewave/internal/g_vector.h"
#include "atomiverse/cpu/planewave/internal/poisson.h"
#include "atomiverse/cpu/planewave/internal/pw_forces.h"
#include "atomiverse/cpu/planewave/internal/pw_setup.h"
#include "atomiverse/cpu/planewave/internal/pw_stress.h"
#include "atomiverse/cpu/planewave/internal/pw_workspace.h"
#include "atomiverse/cpu/planewave/internal/upf_reader.h"
#include "atomiverse/cpu/planewave/internal/veff_builder.h"
#include "atomiverse/cpu/planewave/internal/xc_builder.h"
#include "atomiverse/cpu/planewave/qe_units.h"

#ifdef ATV_HAS_MPI
#include <mpi.h>
#endif

#include <algorithm>
#include <cmath>
#include <complex>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace ATV {
namespace cpu {
namespace pw {

namespace {

/// @brief MPI layout for the k-point round-robin distribution.
struct MpiLayout {
    int rank = 0;
    int nprocs = 1;
    bool active = false; ///< true only when MPI is built AND initialized.
};

MpiLayout mpi_layout() {
    MpiLayout layout;
#ifdef ATV_HAS_MPI
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (initialized) {
        MPI_Comm_rank(MPI_COMM_WORLD, &layout.rank);
        MPI_Comm_size(MPI_COMM_WORLD, &layout.nprocs);
        layout.active = (layout.nprocs > 1);
    }
#endif
    return layout;
}

#ifdef ATV_HAS_MPI
MPI_Datatype mpi_real_type() {
    return sizeof(AccReal) == sizeof(double) ? MPI_DOUBLE : MPI_FLOAT;
}
#endif

} // namespace

// ── PW ────────────────────────────────────────────────────────────

namespace {
void ensure_helpers(PW& self,
                    std::unique_ptr<internal::PoissonSolver>& poisson,
                    std::unique_ptr<internal::XCBuilder>& xc,
                    std::unique_ptr<internal::VeffBuilder>& veff,
                    std::unique_ptr<internal::DensityBuilder>& density) {
    if (!poisson)
        poisson = std::make_unique<internal::PoissonSolver>();
    if (!xc)
        xc = std::make_unique<internal::XCBuilder>();
    if (!veff)
        veff = std::make_unique<internal::VeffBuilder>();
    if (!density)
        density = std::make_unique<internal::DensityBuilder>(
            AccReal(self.mix_beta), self.pulay_depth);
}
} // namespace

// Helpers are allocated on first build()/run() from current member parameters.

// Out-of-line where internal::* are complete types (required for the
// unique_ptr members over forward-declared helpers).
PW::PW() = default;
PW::~PW() = default;
PW::PW(PW&&) noexcept = default;
PW& PW::operator=(PW&&) noexcept = default;

void PW::build(const cpu::Frame& reference) {
    if (!poisson_) {
        poisson_ = std::make_unique<internal::PoissonSolver>();
        xc_ = std::make_unique<internal::XCBuilder>();
        veff_ = std::make_unique<internal::VeffBuilder>();
        density_ = std::make_unique<internal::DensityBuilder>(
            AccReal(mix_beta), pulay_depth);
    }
    if (!upf_path.empty()) {
        // Validate the pseudopotential up front (throws on parse error);
        // setup_pw_problem() re-reads it at each run().
        (void)internal::read_upf_v1(upf_path);
    } else {
        // Validate the toy valence table covers every atom type up front
        // (throws on unsupported Z — no silent 0-electron fallback).
        const int* __restrict sp_type = reference.span(kAtom_type);
        for (size_t i = 0; i < reference.nlocal(); ++i)
            (void)internal::valence_from_z(sp_type[i]);
    }
    built_ = true;
}

PWResult& PW::run(const cpu::Frame& frame) {
    // Lazy build — convenience for callers that skip explicit build().
    if (!built_)
        build(frame);

    result_ = PWResult{};

    const auto layout = mpi_layout();
    const bool log_rank0 = (layout.rank == 0);

    // ── Shared problem assembly (identical to the CUDA backend) ──
    internal::PWProblem prob;
    internal::setup_pw_problem(frame, upf_path, ecut_ry,
                               k_mesh, prob);
    internal::PWWorkspace& ws = prob.ws;
    auto& gv0 = prob.gv0;
    const int n_kpts = static_cast<int>(prob.kg.kpoints.size());
    const int n_bands = prob.n_bands;
    const AccReal vol_bohr = prob.vol_bohr;
    ws.mix_beta = AccReal(mix_beta);
    ws.pulay_depth = pulay_depth;

    if (log_rank0)
        ATV_INFO("[pw] start atoms={} elec={} bands={} ecut_ry={} {} "
                 "k_mesh={} ranks={}",
                 frame.nlocal(), prob.n_elec, n_bands, ecut_ry,
                 upf_path.empty() ? "toy" : "upf", k_mesh,
                 layout.nprocs);

    // ── Step-helper preparation (grid-sized scratch) ──
    xc_->prepare(ws);
    veff_->prepare(ws);
    density_->prepare(ws);
    density_->reset(ws);
    if (!davidson_ || davidson_bands_ != n_bands) {
        davidson_ = std::make_unique<internal::DavidsonSolver>(
            n_bands, AccReal(tol) * AccReal(0.1), 40);
        davidson_bands_ = n_bands;
    }
    davidson_->prepare(ws);
    internal::ConvergenceCheck convergence{AccReal(tol), 2};

    // ── SCF loop ──
    bool converged = false;
    for (int iter = 0; iter < max_iter; ++iter) {
        poisson_->step(ws);
        xc_->step(ws);
        veff_->step(ws);

        // Accumulate density and band energy over this rank's k-points.
        density_->begin(ws);
        AccReal e_band = 0;
        auto* saved_gvec = ws.gvec;
        for (int ik = layout.rank; ik < n_kpts; ik += layout.nprocs) {
            auto& gvk = prob.gv_k[ik];

            // Swap the workspace onto this k-point's G-sphere. T_G must
            // be sized to the k-sphere (its G-count can exceed Γ's).
            ws.gvec = &gvk;
            ws.n_g = gvk.count();
            ws.T_G.resize({static_cast<std::size_t>(ws.n_g)});
            for (int ig = 0; ig < ws.n_g; ++ig)
                ws.T_G[ig] = gvk.ekin[ig]; // 0.5*|k+G|², Cartesian k

            // Nonlocal/wfc tables live on q = k+G: bind this k's tables.
            if (!prob.kb_proj_k.empty())
                ws.kb_proj = prob.kb_proj_k[ik];
            if (!prob.atom_wfc_k.empty())
                ws.atom_wfc = prob.atom_wfc_k[ik];

            davidson_->step(ws);

            const AccReal weight = prob.kg.weights[ik];
            for (int b = 0; b < n_bands; ++b)
                e_band += weight * ws.eigenvalues[b] * AccReal(2);

            density_->accumulate_kpoint(ws, weight, n_bands);

            // Restore the density G-sphere.
            ws.gvec = saved_gvec;
            ws.n_g = gv0.count();
        }

#ifdef ATV_HAS_MPI
        if (layout.active) {
            MPI_Allreduce(MPI_IN_PLACE, density_->accumulation(), ws.nr,
                          mpi_real_type(), MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &e_band, 1, mpi_real_type(), MPI_SUM,
                          MPI_COMM_WORLD);
        }
#endif

        // FFT the (globally summed) density to G-space; Pulay mixing via
        // the shared internal::pulay_mix.
        density_->finish(ws);

        converged = convergence.step(ws);
        result_.iterations = iter + 1;

        const auto terms =
            internal::iteration_energy(ws, e_band, prob.ew_total, vol_bohr);
        result_.energy_ha = static_cast<double>(terms.total);
        result_.energy = result_.energy_ha / qe_units::rydberg_to_hartree;

        // Shared rank-0-gated formatter (same line as the CUDA loop).
        internal::log_iteration_energy(terms, ws, prob.ew_total, iter);

        if (converged) {
            result_.converged = true;
            if (log_rank0)
                ATV_INFO("[pw] converged iter={} energy_ry={:.8f}", iter + 1,
                         result_.energy);
            break;
        }
    }

    if (!result_.converged && log_rank0)
        ATV_WARN("[pw] not converged after {} iters energy_ry={:.8f}",
                 max_iter, result_.energy);

    return result_;
}

namespace {
PW copy_params(const PW& m) {
    PW p;
    p.upf_path = m.upf_path;
    p.ecut_ry = m.ecut_ry;
    p.k_mesh = m.k_mesh;
    p.tol = m.tol;
    p.max_iter = m.max_iter;
    p.mix_beta = m.mix_beta;
    p.pulay_depth = m.pulay_depth;
    return p;
}
} // namespace

std::vector<std::array<double, 3>> pw_forces_analytic(const cpu::Frame& frame,
                                                      const PW& method) {
    return internal::pw_forces_analytic_impl(frame, method);
}

std::array<double, 6> pw_stress(const cpu::Frame& frame, const PW& method) {
    return internal::pw_stress_analytic_impl(frame, method);
}

std::vector<std::array<double, 3>>
pw_forces_fd(const cpu::Frame& frame, const PW& method, double h_ang) {
    if (h_ang <= 0.0)
        throw std::invalid_argument("pw_forces_fd: h_ang must be > 0");
    const int n = frame.nlocal();
    if (n <= 0)
        throw std::invalid_argument("pw_forces_fd: frame has no atoms");

    // Working copy: only positions are perturbed; type/mass/box stay fixed.
    cpu::Frame work = frame;

    auto x = work.span(kAtom_x);
    auto y = work.span(kAtom_y);
    auto z = work.span(kAtom_z);
    StorageReal* axes[3] = {x.data(), y.data(), z.data()};

    std::vector<std::array<double, 3>> forces(static_cast<std::size_t>(n),
                                              {0.0, 0.0, 0.0});
    // Energy is in Ry; convert dE/dR to eV/Å for the Frame force convention.
    const double ry_to_ev = qe_units::rydberg_to_ev;

    for (int a = 0; a < n; ++a) {
        for (int k = 0; k < 3; ++k) {
            const StorageReal r0 = axes[k][a];

            axes[k][a] = r0 + static_cast<StorageReal>(h_ang);
            PW plus = copy_params(method);
            plus.build(work);
            const double e_plus = plus.run(work).energy; // Ry

            axes[k][a] = r0 - static_cast<StorageReal>(h_ang);
            PW minus = copy_params(method);
            minus.build(work);
            const double e_minus = minus.run(work).energy;

            axes[k][a] = r0;
            // F = -dE/dR with E in eV, R in Å
            forces[static_cast<std::size_t>(a)][static_cast<std::size_t>(k)] =
                -(e_plus - e_minus) * ry_to_ev / (2.0 * h_ang);
        }
    }
    return forces;
}

} // namespace pw
} // namespace cpu
} // namespace ATV
