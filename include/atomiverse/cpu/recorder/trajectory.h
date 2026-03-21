/// @file cpu/recorder/trajectory.h
/// @brief CPU-side plain-XYZ trajectory frame writer.
///
/// Called by cuda::Recorder after async D2H of positions.
/// The CUDA component owns no file-format logic — it is here.
#pragma once

#include <fstream>

#include "atomiverse/core/int.h"
#include "atomiverse/core/real.h"

namespace ATV {
namespace cpu {

/// @brief Append one plain-XYZ frame to an open file stream.
///
/// Format per frame:
/// @code
///   N
///   Step <step>
///   <type_id> <x> <y> <z>
///   ...
/// @endcode
///
/// @param os       Open output stream.
/// @param type_id  Atom type IDs [n_atoms] (written as integer species label).
/// @param x        x positions [n_atoms, Angstrom].
/// @param y        y positions [n_atoms, Angstrom].
/// @param z        z positions [n_atoms, Angstrom].
/// @param n_atoms  Number of atoms.
/// @param step     Current step number (written in comment line).
void trajectory_append(std::ofstream& os,
                       const Int* type_id,
                       const Real* x, const Real* y, const Real* z,
                       Int n_atoms, int step);

} // namespace cpu
} // namespace ATV
