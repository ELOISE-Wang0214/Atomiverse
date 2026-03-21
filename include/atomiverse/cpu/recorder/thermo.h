/// @file cpu/recorder/thermo.h
/// @brief CPU-side thermodynamic output (LAMMPS thermo_style one).
///
/// Called by cuda::Thermo after async D2H of pe/ke scalars.
/// The CUDA component owns no formatting or temperature logic — it is here.
#pragma once

#include <fstream>

#include "atomiverse/core/int.h"
#include "atomiverse/core/real.h"

namespace ATV {
namespace cpu {

/// @brief Print one thermodynamic output line.
///
/// Prints a column header on the first call (guarded by header_printed).
/// Writes to stdout and, if file is open, to the optional log file as well.
///
/// @param step           Current step number.
/// @param pe             Potential energy [eV].
/// @param ke             Kinetic energy [eV].
/// @param n_atoms        Number of atoms (used to compute degrees of freedom).
/// @param header_printed Guard flag; set to true after the first call.
/// @param file           Optional open output stream (may be a closed fstream).
void thermo_print(int step, Real pe, Real ke, Int n_atoms,
                  bool& header_printed, std::ofstream& file);

} // namespace cpu
} // namespace ATV
