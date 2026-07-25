/// @file cpu/planewave/internal/restricted_scf.h
/// @brief Shared restricted fixed-occupation PW SCF loop for post-processors.
///
/// Extracted from `bands_dos.cpp` at its second caller
/// (`berry_polarization.cpp`). Runs the plain restricted SCF on
/// `method.k_mesh` and leaves the mixed density and a consistent
/// \(V_\mathrm{eff}\) on `prob` for NSCF diagonalization at arbitrary
/// k-lists (bands, DOS, Berry-phase strings).
#pragma once

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/planewave/internal/pw_setup.h"
#include "atomiverse/cpu/planewave/pw_solver.h"

namespace ATV {
namespace cpu {
namespace pw {
namespace internal {

/// @brief Restricted fixed-occupation SCF; leaves final density + V_eff on
///        `prob` (frozen for subsequent `rebuild_k_tables` + NSCF).
/// @param frame  Geometry: positions [Å], types (Z), simulation cell.
/// @param method Plane-wave method members (`ecut_ry`, `k_mesh`, `xc`, …).
/// @param prob   Problem to assemble and converge (must outlive NSCF use).
/// @return SCF result (energy [Ry]/[Ha], convergence flag, iterations).
PWResult run_restricted_scf(const cpu::Frame& frame, const PW& method,
                            PWProblem& prob);

} // namespace internal
} // namespace pw
} // namespace cpu
} // namespace ATV
