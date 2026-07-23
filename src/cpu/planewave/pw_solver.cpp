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

    if (spin_polarized) {
        if (multiplicity < 1)
            throw std::invalid_argument(
                "PW::run: multiplicity must be >= 1 when spin_polarized");
    }

    const auto layout = mpi_layout();
    const bool log_rank0 = (layout.rank == 0);

    // ── Shared problem assembly (identical to the CUDA backend) ──
    internal::PWProblem prob;
    internal::setup_pw_problem(frame, upf_path, ecut_ry,
                               k_mesh, prob);
    internal::PWWorkspace& ws = prob.ws;
    auto& gv0 = prob.gv0;
    const int n_kpts = static_cast<int>(prob.kg.kpoints.size());
    const int n_elec = prob.n_elec;
    const AccReal vol_bohr = prob.vol_bohr;
    ws.mix_beta = AccReal(mix_beta);
    ws.pulay_depth = pulay_depth;

    // Collinear occupation split: N↑ − N↓ = multiplicity − 1, N↑ + N↓ = n_elec.
    int n_alpha = n_elec / 2;
    int n_beta = n_elec - n_alpha;
    if (spin_polarized) {
        const int mag = multiplicity - 1; // 2S
        if ((n_elec + mag) % 2 != 0 || mag > n_elec)
            throw std::invalid_argument(
                "PW::run: n_elec and multiplicity have incompatible parity "
                "(need N↑,N↓ non-negative integers)");
        n_alpha = (n_elec + mag) / 2;
        n_beta = (n_elec - mag) / 2;
        ws.nspin = 2;
        result_.magnetization = static_cast<double>(mag);
        const std::size_t ng = static_cast<std::size_t>(ws.n_g);
        ws.rho_spin_G.resize(2);
        ws.rho_spin_G_old.resize(2);
        ws.V_xc_spin_G.resize(2);
        for (int s = 0; s < 2; ++s) {
            ws.rho_spin_G[s].resize({ng});
            ws.rho_spin_G_old[s].resize({ng});
            ws.V_xc_spin_G[s].resize({ng});
        }
        // Seed spin densities by splitting the atomic total density.
        const AccReal fa =
            (n_elec > 0) ? AccReal(n_alpha) / AccReal(n_elec) : AccReal(0.5);
        const AccReal fb =
            (n_elec > 0) ? AccReal(n_beta) / AccReal(n_elec) : AccReal(0.5);
        for (int ig = 0; ig < ws.n_g; ++ig) {
            ws.rho_spin_G[0][ig] = fa * ws.rho_G[ig];
            ws.rho_spin_G[1][ig] = fb * ws.rho_G[ig];
            ws.rho_spin_G_old[0][ig] = ws.rho_spin_G[0][ig];
            ws.rho_spin_G_old[1][ig] = ws.rho_spin_G[1][ig];
        }
    } else {
        ws.nspin = 1;
        result_.magnetization = 0;
    }

    const int n_bands_a = std::max(n_alpha, 1);
    const int n_bands_b = std::max(n_beta, 1);
    const int n_bands =
        spin_polarized ? std::max(n_bands_a, n_bands_b) : prob.n_bands;
    // Restricted band count stays n_elec/2 from setup.
    const int n_bands_r = spin_polarized ? n_bands : prob.n_bands;

    if (log_rank0)
        ATV_INFO("[pw] start atoms={} elec={} bands={} ecut_ry={} {} "
                 "k_mesh={} ranks={} spin={} mult={} Na={} Nb={}",
                 frame.nlocal(), n_elec, n_bands_r, ecut_ry,
                 upf_path.empty() ? "toy" : "upf", k_mesh, layout.nprocs,
                 spin_polarized ? 2 : 1, multiplicity, n_alpha, n_beta);

    // ── Step-helper preparation (grid-sized scratch) ──
    xc_->prepare(ws);
    veff_->prepare(ws);
    density_->prepare(ws);
    density_->reset(ws);
    if (!davidson_ || davidson_bands_ != n_bands_r) {
        davidson_ = std::make_unique<internal::DavidsonSolver>(
            n_bands_r, AccReal(tol) * AccReal(0.1), 40);
        davidson_bands_ = n_bands_r;
    }
    davidson_->prepare(ws);
    internal::ConvergenceCheck convergence{AccReal(tol), 2};

    // ── SCF loop ──
    bool converged = false;
    for (int iter = 0; iter < max_iter; ++iter) {
        AccReal e_band = 0;
        auto* saved_gvec = ws.gvec;

        if (spin_polarized) {
            // Total density for Poisson / energy from spin channels.
            for (int ig = 0; ig < ws.n_g; ++ig)
                ws.rho_G[ig] = ws.rho_spin_G[0][ig] + ws.rho_spin_G[1][ig];

            poisson_->step(ws);
            xc_->step_spin(ws);

            const int n_occ[2] = {n_alpha, n_beta};
            // Solve max(N↑,N↓) states per spin; occupy only n_occ[s].
            if (davidson_bands_ != n_bands_r) {
                davidson_ = std::make_unique<internal::DavidsonSolver>(
                    n_bands_r, AccReal(tol) * AccReal(0.1), 40);
                davidson_bands_ = n_bands_r;
                davidson_->prepare(ws);
            }
            xt::xarray<std::complex<AccReal>> rho_out[2];
            for (int s = 0; s < 2; ++s) {
                // Veff_σ = V_loc + V_H + Vxc_σ
                for (int ig = 0; ig < ws.n_g; ++ig)
                    ws.V_xc_G[ig] = ws.V_xc_spin_G[s][ig];
                veff_->step(ws);

                density_->begin(ws);
                for (int ik = layout.rank; ik < n_kpts; ik += layout.nprocs) {
                    auto& gvk = prob.gv_k[ik];
                    ws.gvec = &gvk;
                    ws.n_g = gvk.count();
                    ws.T_G.resize({static_cast<std::size_t>(ws.n_g)});
                    for (int ig = 0; ig < ws.n_g; ++ig)
                        ws.T_G[ig] = gvk.ekin[ig];
                    if (!prob.kb_proj_k.empty())
                        ws.kb_proj = prob.kb_proj_k[ik];
                    if (!prob.atom_wfc_k.empty())
                        ws.atom_wfc = prob.atom_wfc_k[ik];

                    davidson_->step(ws);

                    const AccReal weight = prob.kg.weights[ik];
                    for (int b = 0; b < n_occ[s]; ++b)
                        e_band += weight * ws.eigenvalues[b]; // spin-orbital

                    density_->accumulate_kpoint(ws, weight, n_occ[s],
                                                AccReal(1));

                    ws.gvec = saved_gvec;
                    ws.n_g = gv0.count();
                }

#ifdef ATV_HAS_MPI
                if (layout.active) {
                    MPI_Allreduce(MPI_IN_PLACE, density_->accumulation(),
                                  ws.nr, mpi_real_type(), MPI_SUM,
                                  MPI_COMM_WORLD);
                }
#endif
                density_->extract(ws, rho_out[s]);
            }

#ifdef ATV_HAS_MPI
            if (layout.active) {
                MPI_Allreduce(MPI_IN_PLACE, &e_band, 1, mpi_real_type(),
                              MPI_SUM, MPI_COMM_WORLD);
            }
#endif

            // Linear mix each spin channel (simple, stable for molecules).
            for (int s = 0; s < 2; ++s) {
                for (int ig = 0; ig < ws.n_g; ++ig) {
                    ws.rho_spin_G_old[s][ig] = ws.rho_spin_G[s][ig];
                    ws.rho_spin_G[s][ig] =
                        ws.rho_spin_G_old[s][ig]
                        + AccReal(mix_beta)
                              * (rho_out[s][ig] - ws.rho_spin_G_old[s][ig]);
                }
            }
            for (int ig = 0; ig < ws.n_g; ++ig) {
                ws.rho_G_old[ig] = ws.rho_G[ig];
                ws.rho_G[ig] = ws.rho_spin_G[0][ig] + ws.rho_spin_G[1][ig];
            }

            // Convergence on total density change.
            AccReal dr = 0;
            for (int ig = 0; ig < ws.n_g; ++ig) {
                auto d = ws.rho_G[ig] - ws.rho_G_old[ig];
                dr += d.real() * d.real() + d.imag() * d.imag();
            }
            converged = (std::sqrt(dr) < AccReal(tol));
            result_.iterations = iter + 1;

            // Spin energy: e_one = e_band - 2 EH - ∫(ρ↑Vxc↑+ρ↓Vxc↓)
            auto terms = internal::iteration_energy(ws, e_band, prob.ew_total,
                                                    vol_bohr);
            // Replace restricted e_xc_int (used ρ·Vxc_avg) with spin integral.
            terms.e_one = e_band - AccReal(2.0) * terms.e_h - ws.e_xc_int_spin;
            terms.total =
                terms.e_one + terms.e_h + ws.E_xc + prob.ew_total + ws.E_aug;
            result_.energy_ha = static_cast<double>(terms.total);
            result_.energy =
                result_.energy_ha / qe_units::rydberg_to_hartree;
            internal::log_iteration_energy(terms, ws, prob.ew_total, iter);

        } else {
            // ── Restricted (spin-unpolarized) path — unchanged ──
            poisson_->step(ws);
            xc_->step(ws, xc);
            veff_->step(ws);

            density_->begin(ws);
            e_band = 0;
            for (int ik = layout.rank; ik < n_kpts; ik += layout.nprocs) {
                auto& gvk = prob.gv_k[ik];
                ws.gvec = &gvk;
                ws.n_g = gvk.count();
                ws.T_G.resize({static_cast<std::size_t>(ws.n_g)});
                for (int ig = 0; ig < ws.n_g; ++ig)
                    ws.T_G[ig] = gvk.ekin[ig];
                if (!prob.kb_proj_k.empty())
                    ws.kb_proj = prob.kb_proj_k[ik];
                if (!prob.atom_wfc_k.empty())
                    ws.atom_wfc = prob.atom_wfc_k[ik];

                davidson_->step(ws);

                const AccReal weight = prob.kg.weights[ik];
                for (int b = 0; b < n_bands_r; ++b)
                    e_band += weight * ws.eigenvalues[b] * AccReal(2);

                density_->accumulate_kpoint(ws, weight, n_bands_r);

                ws.gvec = saved_gvec;
                ws.n_g = gv0.count();
            }

#ifdef ATV_HAS_MPI
            if (layout.active) {
                MPI_Allreduce(MPI_IN_PLACE, density_->accumulation(), ws.nr,
                              mpi_real_type(), MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(MPI_IN_PLACE, &e_band, 1, mpi_real_type(),
                              MPI_SUM, MPI_COMM_WORLD);
            }
#endif

            density_->finish(ws);

            converged = convergence.step(ws);
            result_.iterations = iter + 1;

            const auto terms =
                internal::iteration_energy(ws, e_band, prob.ew_total, vol_bohr);
            result_.energy_ha = static_cast<double>(terms.total);
            result_.energy =
                result_.energy_ha / qe_units::rydberg_to_hartree;
            internal::log_iteration_energy(terms, ws, prob.ew_total, iter);
        }

        if (converged) {
            result_.converged = true;
            if (log_rank0)
                ATV_INFO("[pw] converged iter={} energy_ry={:.8f} mag={}",
                         iter + 1, result_.energy, result_.magnetization);
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
    p.spin_polarized = m.spin_polarized;
    p.multiplicity = m.multiplicity;
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
