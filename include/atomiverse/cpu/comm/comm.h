/// @file cpu/comm/comm.h
/// @brief cpu::comm::Comm -- LAMMPS-shaped CommBrick.
///
/// Single class (not a Component subclass) with four entry points
/// matching LAMMPS Comm::exchange / borders / forward_comm / reverse_comm.
/// The component-model adapters that bind these to the MD plan live in a
/// separate header; the logic itself is independent of the plan plumbing.
///
/// Owns a `CommState` PIMPL holding the per-swap state (`sendlist[6]`,
/// `recvnum[6]`, send/recv buffers) and the configured ghost cutoff.
///
/// @see lammps/src/comm_brick.h
/// @see cpu/atom.h, cpu/domain.h, cpu/comm/atom_pack.h
#pragma once

#include <array>
#include <memory>
#include <vector>

namespace ATV {
namespace cpu {

class Atom;
class Box;
class Domain;

namespace comm {

/// @brief Per-instance Comm state (PIMPL). Holds the per-swap sendlist,
/// recv counts, and reusable MPI send/recv buffers.
///
/// Definition lives in `comm.cpp` so the header stays free of MPI-adjacent
/// state. Same pattern as the deleted `halo_exchange.cpp::FaceBufs`.
struct CommState;

/// @brief LAMMPS-shaped CommBrick: borders, exchange, forward_comm, reverse_comm.
///
/// Lifecycle:
///   - construct with `cutghost = max(pair_cutoff + skin, bonded_extension)`
///   - call `build(atom, domain)` once when the driver enters a stage
///     to pre-allocate buffers against the rank's owner count
///   - call `exchange(atom, domain, box)` at every reneighbor event to
///     migrate owners that drifted out of `[sublo, subhi]`
///   - call `borders(atom, domain, box)` immediately after each
///     `exchange` to rebuild sendlists and ghost rows; finishes by
///     rebuilding `atom.map` over `[0, nlocal + nghost)`
///   - call `forward_comm(atom, domain, box)` on non-reneighbor steps to
///     refresh ghost positions cheaply
///   - call `reverse_comm(atom, domain)` after every force pass to gather
///     ghost-row forces back onto their owners
class Comm {
  public:
    /// @brief Construct with the configured ghost cutoff.
    /// @param cutghost Ghost cutoff [Angstrom]. Atoms within `cutghost`
    ///        of a slab face are shipped as ghosts to that face's neighbor.
    explicit Comm(double cutghost);

    ~Comm();
    Comm(Comm&&) noexcept;
    Comm& operator=(Comm&&) noexcept;
    Comm(const Comm&) = delete;
    Comm& operator=(const Comm&) = delete;

    /// @brief Configured ghost cutoff.
    double cutghost() const noexcept { return cutghost_; }

    /// @brief Override the ghost cutoff (e.g. after bonded_extension is
    /// known). Takes effect on the next `borders` call.
    /// @param c New cutoff [Angstrom].
    void set_cutghost(double c) noexcept { cutghost_ = c; }

    /// @brief One-time per-stage init: allocate per-swap buffers against
    /// the rank's current owner count.
    void build(Atom& atom, Domain& domain);

    /// @brief Migrate atoms that drifted out of this rank's slab.
    ///
    /// Mirrors LAMMPS `Comm::exchange()`. On entry expects `atom.nlocal`
    /// owners in `[0, nlocal)` (ghosts ignored). Classifies each by
    /// `[sublo, subhi]` (or `cuts` if RCB-balanced); packs migrants via
    /// `pack_exchange` (carries per-atom topology); single
    /// `MPI_Alltoallv` over byte buffers; gap-fill compaction of survivors
    /// to the front; unpacks received atoms into newly-grown rows. Grows
    /// `atom.nmax` if needed. Clears `atom.map`.
    /// @param atom Atom container.
    /// @param domain Cartesian topology + slab.
    /// @param box Simulation cell.
    void exchange(Atom& atom, Domain& domain, Box& box);

    /// @brief Rebuild per-swap sendlists and ghost rows.
    ///
    /// Mirrors LAMMPS `Comm::borders()`. 6 swaps (2 directions × 3 axes),
    /// cascading: axis-1's scan range includes ghosts received in axis-0,
    /// etc. After unpacking, rebuilds `atom.map` over `[0, nlocal + nghost)`.
    /// Grows `atom.nmax` whenever an incoming ghost batch would overflow.
    /// @param atom Atom container.
    /// @param domain Cartesian topology + slab.
    /// @param box Simulation cell.
    void borders(Atom& atom, Domain& domain, Box& box);

    /// @brief Refresh ghost positions only (no sendlist rebuild).
    ///
    /// Mirrors LAMMPS `Comm::forward_comm()`. Reuses `sendlist` from the
    /// last `borders`. 3 doubles per atom. Ghost slot indices stay stable,
    /// so the receiver writes back into the same ghost rows it filled at
    /// borders time.
    /// @param atom Atom container.
    /// @param domain Cartesian topology + slab.
    /// @param box Simulation cell.
    void forward_comm(Atom& atom, Domain& domain, Box& box);

    /// @brief Gather ghost-row forces back onto their owners.
    ///
    /// Mirrors LAMMPS `Comm::reverse_comm()`. Reverse order of the borders
    /// cascade (z → y → x). For each swap: pack the ghost rows currently
    /// at this rank for that swap; MPI send; receiver adds the forces
    /// into its owner rows at `sendlist` indices. Iterations cumulatively
    /// flow ghost contributions back to owners.
    /// @param atom Atom container.
    /// @param domain Cartesian topology.
    void reverse_comm(Atom& atom, Domain& domain);

    /// @brief Test-only: per-face send-buf capacity (doubles) after `build`.
    /// Used by regression tests to verify that two `Comm` instances do not
    /// share buffer storage.
    std::size_t send_buf_capacity_for_test() const;

  private:
    double cutghost_;
    std::unique_ptr<CommState> state_;
};

} // namespace comm
} // namespace cpu
} // namespace ATV
