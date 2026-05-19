/// @file cpu/pw/pw_scf.h
/// @brief pw::pw_scf() — generic plane-wave DFT total energy.
#pragma once

#include "atomiverse/core/real.h"
#include "atomiverse/cpu/frame.h"
#include <string>

namespace ATV {
namespace cpu {
namespace pw {

/// @brief Result of a plane-wave DFT total energy calculation.
struct PWResult {
    double energy = 0;       ///< Total energy [Ry].
    bool converged = false;  ///< Whether SCF converged within max_iter.
    int iterations = 0;      ///< Actual SCF iterations taken.
};

/// @brief Run plane-wave LDA total energy calculation for a periodic system.
///
/// The Frame provides atom positions, species (type = Z), and simulation cell.
/// Electron count and valence charges are derived from atom types and UPF data.
///
/// @param frame      Input structure: positions, atom types (Z), simulation cell.
/// @param upf_path   Path to UPF pseudopotential file (empty = toy V_loc).
/// @param ecut_ry    Plane-wave energy cutoff [Ry].
/// @param k_mesh     Monkhorst-Pack mesh size per dimension.
/// @param tol        SCF convergence tolerance on density change.
/// @param max_iter   Maximum SCF iterations.
/// @return PWResult with energy (Ry), converged, iterations.
PWResult pw_scf(const cpu::Frame& frame,
                const std::string& upf_path = "",
                double ecut_ry = 30.0,
                int k_mesh = 4,
                double tol = 1e-6,
                int max_iter = 50);

} // namespace pw
} // namespace cpu
} // namespace ATV
