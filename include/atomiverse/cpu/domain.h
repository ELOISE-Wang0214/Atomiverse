/// @file cpu/domain.h
/// @brief cpu::Domain -- LAMMPS-shaped spatial decomposition.
///
/// Replaces the deleted `cpu::md::mpi::DomainDecomposition`. Strictly the
/// LAMMPS `Domain` + `Comm::procgrid/myloc/procneigh` surface: this rank's
/// Cartesian-grid position, its owned spatial slab, and a single per-step
/// reneighbor flag set by the neighbor list and consumed by the comm
/// components. *No* atom counts, *no* tag→local map, *no* generation
/// counters — those move to `cpu::Atom`.
///
/// Lives at the backend root (sibling of `cpu::Box`, `cpu::Atom`), not
/// under `md::mpi`, because the cuts/sublo/subhi geometry is shared by
/// every spatial subsystem (md/, observer/, scf/ if it ever cares about
/// boxes).
///
/// @see cpu/atom.h for the per-atom counts and tag→local map
/// @see cpu/comm/comm.h for the components that consume reneigh_this_step
#pragma once

#include <array>
#include <vector>

#ifdef ATV_HAS_MPI
#include <mpi.h>
#endif

namespace ATV {
namespace cpu {

class Box;

/// @brief MPI Cartesian topology + this rank's spatial slab.
///
/// Construction defaults: dims `{1,1,1}`, periodic on all axes, cart_comm
/// null. Single-rank consumers can read every field without an
/// `mpi_size > 1` guard — the defaults are the single-rank semantics.
class Domain {
  public:
    // ── Cartesian topology (LAMMPS Comm::procgrid / myloc / procneigh) ─────

#ifdef ATV_HAS_MPI
    /// @brief Cartesian sub-communicator. `MPI_COMM_NULL` in single-rank mode.
    /// Set by `init_topology`; in this implementation it aliases the world
    /// communicator (see `init_topology` rationale on OpenMPI 5.0.x).
    MPI_Comm cart_comm = MPI_COMM_NULL;
#endif

    /// @brief Cartesian process grid dimensions. Single-rank: `{1, 1, 1}`.
    int dims[3] = {1, 1, 1};

    /// @brief This rank's coordinates in the Cartesian grid.
    int coords[3] = {0, 0, 0};

    /// @brief Periodic flags per axis (1 = periodic, 0 = non-periodic).
    int periods[3] = {1, 1, 1};

    /// @brief Face neighbor ranks in `{-x, +x, -y, +y, -z, +z}` order.
    /// Single-rank: all zeros (this rank is its own periodic neighbor).
    std::array<int, 6> face_neighbors{};

    // ── This rank's owned spatial slab (LAMMPS Domain::sublo/subhi) ────────

    /// @brief Low corner of this rank's slab [Angstrom].
    double sublo[3] = {};

    /// @brief High corner of this rank's slab [Angstrom].
    double subhi[3] = {};

    // ── Optional RCB-balanced cut planes (load balancer) ───────────────────

    /// @brief Per-axis cut planes as fractions of the box length, in `[0, 1]`.
    /// Empty: uniform brick decomposition. Size `dims[d] + 1` with
    /// `cuts[d].front() == 0` and `cuts[d].back() == 1`: non-uniform cuts
    /// set by RCBBalancer. Fractions survive NPT/Deform box rescaling.
    std::array<std::vector<double>, 3> cuts{};

    // ── Per-step reneighbor signal ─────────────────────────────────────────

    /// @brief Set by `NeighborList::launch` at the start of `pre_force`
    /// when reneighboring is needed (max displacement > skin/2, or
    /// scheduled-frequency boundary). Consumed by
    /// `cpu::comm::Exchange` / `Borders` / `Forward` to decide whether
    /// to run the full exchange + ghost rebuild or the cheap position
    /// refresh. Replaces both `neighbor_generation` and
    /// `halo_generation` counters from the old design.
    bool reneigh_this_step = false;

    // ── Lifecycle ──────────────────────────────────────────────────────────

#ifdef ATV_HAS_MPI
    /// @brief Build the Cartesian topology over `world`.
    ///
    /// Single-rank world (`world_size == 1`): no-op; defaults stand.
    /// Multi-rank: factorises `world_size` via `MPI_Dims_create` (or
    /// honours `ATV_FORCE_DIMS=Nx,Ny,Nz`), computes this rank's
    /// `coords`, and discovers its 6 face neighbors. `cart_comm`
    /// aliases `world` directly (no `MPI_Cart_create`, see note in
    /// the .cpp on OpenMPI 5.0.x macOS teardown).
    /// @param world Parent communicator.
    void init_topology(MPI_Comm world);

    /// @brief Release the Cartesian sub-communicator if it owns one.
    /// Safe to call when `cart_comm == MPI_COMM_NULL`.
    /// @param mpi_is_finalizing If true, skip MPI calls (the runtime is
    ///        shutting down).
    void release_comm(bool mpi_is_finalizing = false);
#endif

    /// @brief Refresh `sublo` / `subhi` against the current box.
    ///
    /// Single-rank: this rank's slab is the full box. Multi-rank: slab
    /// is either the `coords[d] / dims[d]` uniform brick or, when
    /// `cuts[d]` is populated, the corresponding fraction of the box
    /// length along axis `d`.
    /// @param box Simulation cell.
    void compute_sublo_subhi(const Box& box);
};

} // namespace cpu
} // namespace ATV
