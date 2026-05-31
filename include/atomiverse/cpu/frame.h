/// @file cpu/frame.h
/// @brief cpu::Frame -- thin holder of Atom, Box, and SystemScalars.
///
/// **Transitional state.** The nested `Frame::Bonds / Angles / Dihedrals`
/// sub-structs below are *dead storage* — no force kernel reads them
/// anymore (all bonded forces iterate `frame.atoms.num_bond[i]` etc. in
/// the LAMMPS per-atom layout). They remain as compile-shims so the ~30
/// tests + benchmarks that populate `frame.bonds.atomi = xt::zeros(...)`
/// still build until they're individually migrated to the
/// `frame.atoms.add_bond / add_angle / add_dihedral` setter API. Delete
/// in the follow-up "test fixture update" commit once all callers move.
///
/// Post-refactor: Frame is no longer a god struct. The per-atom container
/// is `cpu::Atom`, which carries positions, velocities, forces, type, tag,
/// image, per-atom bonded topology (`num_bond[i] + bond_atom[i][m]`,
/// LAMMPS style), and the `tag → local` map. Decomposition state
/// (`cpu::md::mpi::DomainDecomposition`) lives on `cpu::Backend`, not on
/// Frame.
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
/// @see cpu/atom.h, cpu/box.h, cpu/md/mpi/domain.h, cpu/comm/halo_exchange.h
/// @see core/system_scalars.h
#pragma once

#include <xtensor/containers/xtensor.hpp>

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

    // ── Transitional dead-storage shims ────────────────────────────────
    //
    // These nested types exist only so the ~30 tests + benchmarks that
    // pre-date the per-atom topology refactor still compile. They are
    // *write-only from the test side* and *read-only by nothing else* —
    // every force kernel and every halo / migration / NL path now reads
    // from `frame.atoms.num_bond / bond_atom / angle_atom / dihedral_atom`
    // (the LAMMPS per-atom layout). Materialize and the LAMMPS-data
    // reader populate `frame.atoms` directly; the legacy `Bonds / Angles
    // / Dihedrals` blocks they used to fill are no longer populated.
    //
    // Delete these sub-structs entirely once the test fixture cleanup
    // commit lands.

    struct Bonds {
        xt::xtensor<int, 1> atomi;
        xt::xtensor<int, 1> atomj;
        xt::xtensor<int, 1> type;
        int size() const { return static_cast<int>(atomi.size()); }
    } bonds;

    struct Angles {
        xt::xtensor<int, 1> atomi;
        xt::xtensor<int, 1> atomj;
        xt::xtensor<int, 1> atomk;
        xt::xtensor<int, 1> type;
        int size() const { return static_cast<int>(atomi.size()); }
    } angles;

    struct Dihedrals {
        xt::xtensor<int, 1> atomi;
        xt::xtensor<int, 1> atomj;
        xt::xtensor<int, 1> atomk;
        xt::xtensor<int, 1> atoml;
        xt::xtensor<int, 1> type;
        int size() const { return static_cast<int>(atomi.size()); }
    } dihedrals;
};

} // namespace cpu
} // namespace ATV
