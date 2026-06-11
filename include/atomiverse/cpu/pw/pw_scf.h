/// @file cpu/pw/pw_scf.h
/// @brief pw::pw_scf() — thin functional shim over cpu::pw::PWSolver.
///
/// Mirrors `cpu::scf::scf_solve()`: constructs a transient `PWSolver`
/// and runs it. Prefer `cpu::Evaluate + cpu::pw::PWEvaluator` (or
/// `PWSolver` directly) in new code.
#pragma once

#include <string>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/pw/pw_solver.h"

namespace ATV {
namespace cpu {
namespace pw {

/// @brief Run a plane-wave LDA total energy calculation for a periodic system.
///
/// The Frame provides atom positions, species (type = Z), and simulation cell.
/// Electron count and valence charges are derived from atom types and UPF data.
///
/// @param frame      Input structure: positions [Å], atom types (Z), cell.
/// @param upf_path   Path to UPF pseudopotential file (empty = toy V_loc).
/// @param ecut_ry    Plane-wave energy cutoff [Ry].
/// @param k_mesh     Monkhorst-Pack mesh size per dimension.
/// @param tol        SCF convergence tolerance on density change.
/// @param max_iter   Maximum SCF iterations.
/// @return PWResult with energy [Ry] (+ energy_ha [Ha]), converged, iterations.
PWResult pw_scf(const cpu::Frame& frame, const std::string& upf_path = "",
                double ecut_ry = 30.0, int k_mesh = 4, double tol = 1e-6,
                int max_iter = 50);

} // namespace pw
} // namespace cpu
} // namespace ATV
