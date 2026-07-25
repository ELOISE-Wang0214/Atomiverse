/// @file cpu/pw/internal/scf_loop.cpp
/// @brief Shared minimal restricted fixed-occupation plane-wave SCF loop.

#include "atomiverse/cpu/planewave/internal/scf_loop.h"

#include "atomiverse/core/log.h"
#include "atomiverse/cpu/planewave/internal/convergence.h"
#include "atomiverse/cpu/planewave/internal/davidson.h"
#include "atomiverse/cpu/planewave/internal/density_builder.h"
#include "atomiverse/cpu/planewave/internal/poisson.h"
#include "atomiverse/cpu/planewave/internal/veff_builder.h"
#include "atomiverse/cpu/planewave/internal/xc_builder.h"
#include "atomiverse/cpu/planewave/qe_units.h"

#ifdef ATV_HAS_MPI
#include <mpi.h>
#endif

namespace ATV {
namespace cpu {
namespace pw {
namespace internal {

namespace {

struct MpiLayout {
    int rank = 0;
    int nprocs = 1;
    bool active = false;
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

PWResult run_scf_fixed_occ(const cpu::Frame& frame, const PW& method,
                           PWProblem& prob, PoissonSolver& poisson,
                           XCBuilder& xc, VeffBuilder& veff,
                           DensityBuilder& density, DavidsonSolver& davidson,
                           const char* tag) {
    PWResult result{};
    const auto layout = mpi_layout();

    setup_pw_problem(frame, method.upf_path, method.ecut_ry, method.k_mesh,
                     prob);
    PWWorkspace& ws = prob.ws;
    auto& gv0 = prob.gv0;
    const int n_kpts = static_cast<int>(prob.kg.kpoints.size());
    const int n_bands = prob.n_bands;
    const AccReal vol_bohr = prob.vol_bohr;
    ws.mix_beta = AccReal(method.mix_beta);
    ws.pulay_depth = method.pulay_depth;

    xc.prepare(ws);
    veff.prepare(ws);
    density.prepare(ws);
    density.reset(ws);
    davidson.prepare(ws);
    ConvergenceCheck convergence{AccReal(method.tol), 2};

    for (int iter = 0; iter < method.max_iter; ++iter) {
        poisson.step(ws);
        xc.step(ws, method.xc);
        veff.step(ws);

        density.begin(ws);
        AccReal e_band = 0;
        auto* saved_gvec = ws.gvec;
        for (int ik = layout.rank; ik < n_kpts; ik += layout.nprocs) {
            auto& gvk = prob.gv_k[static_cast<std::size_t>(ik)];
            ws.gvec = &gvk;
            ws.n_g = gvk.count();
            ws.T_G.resize({static_cast<std::size_t>(ws.n_g)});
            for (int ig = 0; ig < ws.n_g; ++ig)
                ws.T_G[ig] = gvk.ekin[ig];

            if (!prob.kb_proj_k.empty())
                ws.kb_proj = prob.kb_proj_k[static_cast<std::size_t>(ik)];
            if (!prob.atom_wfc_k.empty())
                ws.atom_wfc = prob.atom_wfc_k[static_cast<std::size_t>(ik)];

            davidson.step(ws);

            const AccReal weight =
                prob.kg.weights[static_cast<std::size_t>(ik)];
            for (int b = 0; b < n_bands; ++b)
                e_band += weight * ws.eigenvalues[b] * AccReal(2);

            density.accumulate_kpoint(ws, weight, n_bands);

            ws.gvec = saved_gvec;
            ws.n_g = gv0.count();
        }

#ifdef ATV_HAS_MPI
        if (layout.active) {
            MPI_Allreduce(MPI_IN_PLACE, density.accumulation(), ws.nr,
                          mpi_real_type(), MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &e_band, 1, mpi_real_type(), MPI_SUM,
                          MPI_COMM_WORLD);
        }
#endif

        density.finish(ws);
        result.iterations = iter + 1;

        const auto terms =
            iteration_energy(ws, e_band, prob.ew_total, vol_bohr);
        result.energy_ha = static_cast<double>(terms.total);
        result.energy = result.energy_ha / qe_units::rydberg_to_hartree;

        if (convergence.step(ws)) {
            result.converged = true;
            if (layout.rank == 0)
                ATV_INFO("[pw] {}-scf converged iter={} energy_ry={:.8f}", tag,
                         iter + 1, result.energy);
            break;
        }
    }
    if (!result.converged && layout.rank == 0)
        ATV_WARN("[pw] {}-scf not converged after {} iters", tag,
                 method.max_iter);
    return result;
}

} // namespace internal
} // namespace pw
} // namespace cpu
} // namespace ATV
