/// @file cpu/pw/solve_si.h
/// @brief pw::solve_si() — plane-wave DFT total energy for Si crystal.
#pragma once

#include "atomiverse/core/real.h"
#include <string>

namespace ATV {
namespace cpu {
namespace pw {

/// @brief Run plane-wave LDA total energy calculation for Si diamond crystal.
///
/// Uses 4x4x4 Monkhorst-Pack k-point mesh, E_cut from ecut_ry parameter.
/// Total energy returned in Rydberg.
///
/// @param upf_path Path to SG15 ONCV UPF v1 pseudopotential file (if empty, uses toy V_loc).
/// @param ecut_ry Plane-wave energy cutoff [Ry].
/// @param k_mesh Monkhorst-Pack mesh size per dimension.
/// @param tol SCF convergence tolerance.
/// @param max_iter Maximum SCF iterations.
/// @return Total energy [Ry].
double solve_si(const std::string& upf_path = "",
                double ecut_ry = 30.0,
                int k_mesh = 4,
                double tol = 1e-6,
                int max_iter = 50);

} // namespace pw
} // namespace cpu
} // namespace ATV
