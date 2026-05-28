/// @file cpu/frame.h
/// @brief cpu::Frame -- thin holder of Atom, Box, and SystemScalars.
///
/// Post-refactor: Frame is no longer a god struct. The per-atom container
/// is `cpu::Atom`, which carries positions, velocities, forces, type, tag,
/// image, per-atom bonded topology (`num_bond[i] + bond_atom[i][m]`,
/// LAMMPS style), and the `tag → local` map. Decomposition state
/// (`cpu::Domain`) and comm state (`cpu::comm::Comm`) live on `cpu::Backend`,
/// not on Frame.
///
/// The member `Frame::atoms` is named in the plural (the legacy idiom)
/// but its type is `cpu::Atom` (LAMMPS Atom-equivalent) — callsites that
/// previously wrote `frame.atoms.x[i]`, `frame.atoms.tag[i]`,
/// `frame.atoms.size()` continue to work; the same fields and the
/// `size()` accessor (returning `nlocal`) are on `cpu::Atom`.
///
/// The old nested `Frame::Bonds / Angles / Dihedrals` flat-topology
/// sub-structs are deleted. Their data moves into `frame.atoms.num_bond /
/// bond_atom / bond_type / num_angle / ... / nspecial / special`.
/// Callers that read `frame.bonds.atomi[k]` must switch to the per-atom
/// iteration `for (i in [0, nlocal)) for (m in [0, atom.num_bond[i]))`.
///
/// @see cpu/atom.h, cpu/box.h, cpu/domain.h, cpu/comm/comm.h
/// @see core/system_scalars.h
#pragma once

#include "atomiverse/core/system_scalars.h"
#include "atomiverse/cpu/atom.h"
#include "atomiverse/cpu/box.h"

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
