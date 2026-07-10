/// @file cpu/frame/frame.h
/// @brief cpu::Frame -- thin holder of Atom, Box, and SystemScalars.
///
/// Frame is not a god struct. The per-atom container is `cpu::Atom`,
/// which carries positions, velocities, forces, type, tag, image,
/// per-atom bonded topology (`num_bond[i] + bond_atom[i][m]`, LAMMPS
/// style), and the `tag → local` map. Decomposition state
/// (`cpu::md::mpi::DomainDecomposition`) lives on `cpu::Backend`, not on
/// Frame.
///
/// The member `Frame::atoms` is named in the plural (the legacy idiom)
/// but its type is `cpu::Atom` (LAMMPS Atom-equivalent) — callsites
/// write `frame.atoms.x[i]`, `frame.atoms.tag[i]`, `frame.atoms.size()`
/// (returning `nlocal`).
///
/// Topology is per-atom only: populate with
/// `frame.atoms.add_bond / add_angle / add_dihedral` (tag map first —
/// `map_clear()` + `map_set(...)` — then `build_special_list()`), and
/// read with `for (i in [0, nlocal)) for (m in [0, atoms.num_bond[i]))`.
/// The historical flat `Bonds / Angles / Dihedrals` sub-structs are
/// deleted (arch-refactor-06-frame); the molrs boundary keeps the flat
/// Block schema, converted at `materialize()` / `write_back()` only.
///
/// @see cpu/frame/atom.h, cpu/frame/box.h, cpu/frame/materialize.h
/// @see core/system_scalars.h
#pragma once

#include "atomiverse/core/system_scalars.h"
#include "atomiverse/cpu/frame/atom.h"
#include "atomiverse/cpu/frame/box.h"

namespace ATV {
namespace cpu {

/// @brief CPU backend execution frame: thin holder of per-atom data,
/// cell, and scalar results.
class Frame {
  public:
    /// @brief Construct with an initial per-atom array capacity.
    ///
    /// `atoms` is allocated with capacity `n_atoms` and `nlocal = n_atoms`
    /// (single-rank semantics: every row is locally owned). Multi-rank
    /// drivers reduce `nlocal` via the first `cpu::comm::Exchange`.
    /// @param n_atoms Initial atom count / capacity.
    explicit Frame(int n_atoms) : atoms(n_atoms) { atoms.nlocal = n_atoms; }

    /// @brief Per-atom container. Positions, velocities, forces, type,
    /// tag, image, per-atom bonded topology, and the `tag → local` map.
    cpu::Atom atoms;

    /// @brief Per-frame system scalar results (energy, kinetic, conserved,
    /// stress, kinetic stress). See `core/system_scalars.h`.
    SystemScalars scalars{};

    /// @brief Simulation cell.
    cpu::Box box{};
};

} // namespace cpu
} // namespace ATV
