/// @file cpu/pw/internal/scf_loop.h
/// @brief internal::run_scf_fixed_occ — minimal restricted fixed-occupation
///        plane-wave SCF shared by the analytic force / stress / DFPT paths.
///
/// Promoted from the (previously duplicated) file-local `run_scf` helpers in
/// `pw_forces.cpp` and `pw_stress.cpp` on their third caller (DFPT). Honors
/// `method.xc` (the historical `pw_forces.cpp` copy silently evaluated the
/// default functional). Facade-private: never on the public include path.
#pragma once

#include "atomiverse/cpu/planewave/internal/pw_setup.h"
#include "atomiverse/cpu/planewave/pw_solver.h"

namespace ATV {
namespace cpu {
namespace pw {
namespace internal {

class PoissonSolver;
class XCBuilder;
class VeffBuilder;
class DensityBuilder;
class DavidsonSolver;

/// @brief Run the restricted fixed-occupation SCF loop to convergence.
///
/// Assembles the problem via `setup_pw_problem` and iterates
/// Poisson → XC (`method.xc`) → Veff → per-k Davidson → density + Pulay mix.
/// On exit `prob` holds the final mixed density on the Γ sphere; k-points are
/// distributed round-robin over MPI ranks when MPI is initialized.
///
/// @param frame    Geometry [Å]; not mutated.
/// @param method   Plane-wave method members (`ecut_ry`, `k_mesh`, `xc`, …).
/// @param prob     Problem storage to fill (must outlive the caller's use).
/// @param poisson  Hartree step helper.
/// @param xc       XC step helper.
/// @param veff     Effective-potential step helper.
/// @param density  Density accumulation / Pulay mixing helper.
/// @param davidson Band solver (caller chooses band count / tolerance).
/// @param tag      Short label for the rank-0 convergence log line.
/// @return Energy [Ry]/[Hartree], convergence flag, iteration count.
PWResult run_scf_fixed_occ(const cpu::Frame& frame, const PW& method,
                           PWProblem& prob, PoissonSolver& poisson,
                           XCBuilder& xc, VeffBuilder& veff,
                           DensityBuilder& density, DavidsonSolver& davidson,
                           const char* tag);

} // namespace internal
} // namespace pw
} // namespace cpu
} // namespace ATV
