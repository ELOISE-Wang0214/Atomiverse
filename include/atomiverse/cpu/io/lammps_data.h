/// @file cpu/io/lammps_data.h
/// @brief LAMMPS "molecular" atom_style data-file loader.
///
/// Wraps the molrs LAMMPS reader (Rust) and emits a fully populated
/// `cpu::Frame`: per-atom positions/masses/types, plus bonds, angles,
/// dihedrals, and the orthogonal simulation box. Indices in the
/// returned Frame are 0-indexed (LAMMPS files are 1-indexed; the bridge
/// performs the offset).
///
/// The Masses section is parsed here in C++ — molrs only puts per-atom
/// mass into the atoms block when atom_style="body". For
/// atom_style="molecular" (which the MSU PE deck uses) masses are
/// per-type and live in their own section.
#pragma once

#include <string>

namespace ATV {
namespace cpu {

class Frame;

/// @brief Read a LAMMPS atom_style=molecular data file and return a
///        `cpu::Frame` with positions, masses, types, bonds, angles,
///        dihedrals, and the orthogonal box populated.
///
/// @param path Filesystem path to the LAMMPS data file.
/// @return A populated `cpu::Frame`.
/// @throws std::runtime_error if the file cannot be opened or the
///         molrs parser rejects it.
Frame load_lammps_data(const std::string& path);

} // namespace cpu
} // namespace ATV
