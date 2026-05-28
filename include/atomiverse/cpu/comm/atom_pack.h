/// @file cpu/comm/atom_pack.h
/// @brief Pack / unpack helpers for cpu::Atom across MPI swaps.
///
/// Mirrors LAMMPS' `AtomVec::pack_comm / pack_border / pack_exchange`
/// trio. Free functions (not member methods on `cpu::Atom`) because we
/// have a single canonical atom layout and LAMMPS' AtomVec polymorphism
/// doesn't apply.
///
/// Pack layouts (per-atom doubles, dense, contiguous):
///
///   pack_comm     [x' y' z']                                       (3 dbl)
///   pack_border   [x' y' z' tag type img_x img_y img_z]            (8 dbl)
///   pack_exchange  variable-length; carries x/y/z/v/f/mass/charge/type/tag
///                  /image + per-atom topology + nspecial/special list.
///                  Returns the doubles-written count; receiver reads
///                  the same count.
///
/// `x'/y'/z'` denote positions with the per-face PBC shift already
/// applied (caller supplies the shift), so the receiver sees ghosts
/// in its own coordinate frame.
///
/// @see lammps/src/atom_vec.h for the upstream pack/unpack contracts
/// @see cpu/atom.h for the per-atom storage these functions operate on
#pragma once

#include <cstdint>

namespace ATV {
namespace cpu {

class Atom;
class Box;

namespace comm {

// ── Forward (steady-state ghost position refresh) ─────────────────────────

/// @brief Doubles written per atom by `pack_comm`.
inline constexpr int kCommDoublesPerAtom = 3;

/// @brief Pack positions of `n` atoms indexed by `sendlist` into `buf`.
/// Adds the per-face PBC shift to the packed coordinates.
/// @param atom Source atom container.
/// @param sendlist Local row indices to pack (length `n`).
/// @param n Number of atoms to pack.
/// @param sx X-axis PBC shift to add to packed `x` [Angstrom].
/// @param sy Y-axis PBC shift [Angstrom].
/// @param sz Z-axis PBC shift [Angstrom].
/// @param buf Output buffer; must hold at least `n * kCommDoublesPerAtom` doubles.
/// @return Number of doubles written.
int pack_comm(const Atom& atom, const int* sendlist, int n, double sx,
              double sy, double sz, double* buf);

/// @brief Unpack `n` positions from `buf` into ghost slots
/// `[first, first + n)` of `atom`.
/// @param atom Destination atom container; slot range must already be allocated.
/// @param n Number of atoms to unpack.
/// @param first First ghost row to write.
/// @param buf Source buffer.
void unpack_comm(Atom& atom, int n, int first, const double* buf);

// ── Border (full ghost rebuild: position + tag + type + image) ────────────

/// @brief Doubles written per atom by `pack_border`.
inline constexpr int kBorderDoublesPerAtom = 8;

/// @brief Pack position + tag + type + image of `n` atoms into `buf`.
/// Tag and type are cast to double for uniform packing (atom counts up
/// to 2^53 are exact; types up to 2^31 trivially fit). Image components
/// are individually cast.
/// @param atom Source atom container.
/// @param sendlist Local row indices.
/// @param n Number of atoms.
/// @param sx X-axis PBC shift [Angstrom].
/// @param sy Y-axis PBC shift [Angstrom].
/// @param sz Z-axis PBC shift [Angstrom].
/// @param buf Output buffer (≥ `n * kBorderDoublesPerAtom` doubles).
/// @return Number of doubles written.
int pack_border(const Atom& atom, const int* sendlist, int n, double sx,
                double sy, double sz, double* buf);

/// @brief Unpack position + tag + type + image into ghost slots
/// `[first, first + n)` of `atom`.
/// @param atom Destination atom container.
/// @param n Number of atoms.
/// @param first First ghost row to write.
/// @param buf Source buffer.
void unpack_border(Atom& atom, int n, int first, const double* buf);

// ── Exchange (atom migration: full per-atom state + per-atom topology) ────

/// @brief Pack the full per-atom payload for atom `i` into `buf`.
///
/// Payload includes:
///   - x, y, z, vx, vy, vz, fx, fy, fz, mass, charge   (11 doubles)
///   - type, tag                                       (2 ints, packed as doubles)
///   - image[3]                                        (3 ints, packed as doubles)
///   - num_bond + bond_atom[..] + bond_type[..]
///   - num_angle + angle_atom[..] + angle_type[..]
///   - num_dihedral + dihedral_atom[..] + dihedral_type[..]
///   - nspecial[3] + special[..]
///
/// Variable size; the receiver reads the exact same number of doubles.
/// Bonds-per-atom counts are written before their respective partner-tag
/// arrays so the receiver knows how many entries to consume.
///
/// @param atom Source atom container.
/// @param i Owner row index to pack.
/// @param buf Output buffer (must be sized large enough for the worst-case
///            payload — typically `kMaxExchangeDoubles` per atom; see
///            `exchange_doubles_for(atom, i)` for the exact count).
/// @return Number of doubles written.
int pack_exchange(const Atom& atom, int i, double* buf);

/// @brief Unpack one atom payload from `buf` into row `slot`.
/// @param atom Destination atom container; `slot < atom.nmax` must hold.
/// @param slot Local row to write.
/// @param buf Source buffer (must match what `pack_exchange` wrote).
/// @return Number of doubles read.
int unpack_exchange(Atom& atom, int slot, const double* buf);

/// @brief Worst-case doubles for a single `pack_exchange` payload. Used
/// to size send/recv buffers conservatively.
inline constexpr int kMaxExchangeDoubles =
    11               // positions / velocities / forces / mass / charge
    + 5              // type + tag + image[3]
    + 1 + 6 + 6      // num_bond + bond_atom (kMax=6) + bond_type
    + 1 + 12 + 6     // num_angle + angle_atom (2 * kMax=6) + angle_type
    + 1 + 24 + 8     // num_dihedral + dihedral_atom (3 * kMax=8) + dihedral_type
    + 3 + 16;        // nspecial[3] + special (kMax=16)

// ── Bonded skin extent ────────────────────────────────────────────────────

/// @brief Maximum intra-tuple atom-to-atom distance across all bonds,
/// angles, and dihedrals where every participant is locally resolvable.
///
/// Used by `cpu::comm::Borders` to pick a halo cutoff wide enough to
/// keep every bonded partner ghosted. Multi-rank callers reduce the
/// per-rank result with `MPI_Allreduce(MAX)` before applying.
///
/// @param atom Atom container with current positions and per-atom topology.
/// @param box Simulation cell, for minimum-image distance.
/// @return Max pairwise span in Angstroms; 0 if no bonded tuples present
///         or no tuple is fully locally resolved.
double compute_bonded_extension_local(const Atom& atom, const Box& box);

} // namespace comm
} // namespace cpu
} // namespace ATV
