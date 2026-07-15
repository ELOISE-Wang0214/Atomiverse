/// @file cpu/comm/bcast_frame.h
/// @brief `bcast_frame` -- replicate the root rank's whole Frame onto every
///        rank of `MPI_COMM_WORLD`.
///
/// The rank-0-reads / every-rank-computes primitive: the root packs each
/// owned atom row through the registry-driven `comm::pack_exchange` codec
/// (schema order is wire order; exchanged extension columns ride at the row
/// tail), broadcasts a fixed POD header (atom count, payload size, bonded
/// flag, `cpu::Box`, `SystemScalars`) followed by the variable-length row
/// payload, and every non-root rank reconstructs a full replica: all rows
/// owned (`nlocal == n`, `nghost == 0`), tag→local map rebuilt, box and
/// scalars overwritten. The root's own Frame is never modified.
///
/// The payload is pure byte transit -- no arithmetic touches any value in
/// flight, so every replicated field (positions [Angstrom], velocities
/// [Angstrom/fs], forces [eV/Angstrom], charges [e], scalars [eV], ...) is
/// bit-identical to the root's.
///
/// **Collective contract (lockstep):** under multi-rank every rank must call
/// this function in lockstep *before* any rank gate -- it enters MPI
/// collectives internally (same discipline as `gather_for_root_io`).
///
/// **Precondition:** the exchanged extension-column registration
/// (`Frame::insert(..., exchanged=true)`) must be identical on every rank.
/// A divergent registration changes the per-row payload width; the affected
/// non-root rank detects the consumed-doubles mismatch after decoding and
/// throws `std::runtime_error` (never a silent truncation). The collectives
/// complete before decoding, so the failure stays lockstep-safe.
///
/// Degrade path: on a non-MPI build, or when MPI is not initialized, or at
/// world size 1, the call returns immediately and the Frame is untouched
/// (`MPI_Initialized` probe, same as `scatter_atoms_from_root`).
///
/// @see cpu/comm/atom_pack.h    the row codec this primitive rides
/// @see cpu/comm/atom_migration.h  `scatter_atoms_from_root`, the
///      distribute-by-ownership sibling of this replicate-everywhere call
#pragma once

namespace ATV {
namespace cpu {

class Frame;

namespace comm {

/// @brief Broadcast the root rank's Frame to every rank of MPI_COMM_WORLD.
///
/// @param frame On the root rank: the source Frame (owned rows `[0, nlocal)`;
///              read-only, never modified). On every other rank: the
///              destination Frame, overwritten with a full replica (atom
///              rows, bonded topology, exchanged extras, box, scalars).
/// @param root  Rank holding the source Frame (default 0).
/// @throws std::runtime_error on every rank if the payload exceeds
///         `INT_MAX` doubles (no silent chunking), and on a non-root rank
///         whose decode consumes a different double count than the root
///         announced (exchanged-column registration divergence).
void bcast_frame(Frame& frame, int root = 0);

} // namespace comm
} // namespace cpu
} // namespace ATV
