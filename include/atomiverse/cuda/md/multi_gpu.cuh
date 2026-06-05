#pragma once

/// @file cuda/md/multi_gpu.cuh
/// @brief MultiGpu — single-process multi-GPU MD orchestrator. Owns the
///        DeviceManager (P2P over NVLink) and a 3D DomainDecomposition, splits
///        a global cpu::Frame into one device-resident subdomain frame per GPU,
///        and drives the full Verlet-cadence step loop over NVLink: partition /
///        gather, the 3-axis cascade halo (build_halo) + per-step position
///        forward_comm, atom migration across subdomains, the per-subdomain
///        cluster neighbor-list LJ force path, and velocity-Verlet NVE — all
///        without MPI (ghost/migration traffic is cudaMemcpyPeer).
///
/// gather()/gather_forces() merge owned atoms back into a global cpu::Frame by
/// persistent tag (so they survive migration). What is NOT yet here: device-
/// side cross-GPU reduction of energy/virial (gather is host-side) and bonded
/// interactions across subdomain boundaries.

#include <memory>
#include <vector>

#include "atomiverse/cuda/context.cuh"
#include "atomiverse/cuda/device_manager.cuh"
#include "atomiverse/cuda/domain_decomposition.cuh"
#include "atomiverse/cuda/frame.cuh"
#include "atomiverse/cuda/neighbor_list.cuh"
#include "atomiverse/cuda/potentials/pair_lj_cut.cuh"

namespace ATV {
namespace cpu {
class Frame;
}
namespace cuda {
namespace md {

class MultiGpu {
  public:
    /// @param n_gpus  Subdomains / GPUs to use (<=0 ⇒ all visible).
    explicit MultiGpu(int n_gpus = -1);
    ~MultiGpu();
    MultiGpu(const MultiGpu&) = delete;
    MultiGpu& operator=(const MultiGpu&) = delete;

    int count() const { return dm_.count(); }
    const DeviceManager& devices() const { return dm_; }
    const DomainDecomposition& domain() const { return *dd_; }

    /// Split `global` by owner slab and materialize each subdomain onto its
    /// GPU. `ghost_headroom` over-allocates each subdomain's per-atom arrays so
    /// halo ghosts can be appended later without reallocation (capacity =
    /// (n_local * (1 + ghost_headroom)) rounded up). Frees any prior frames.
    void partition(const cpu::Frame& global, double ghost_headroom = 0.5);

    /// Owned atom count on subdomain `slot`.
    int n_local(int slot) const;
    /// Ghost atom count currently appended on subdomain `slot` (after a halo).
    int n_ghost(int slot) const;
    /// The device frame for `slot` (owned in [0,n_local), ghosts in [n_local,
    /// n_local+n_ghost)).
    const cuda::Frame& frame(int slot) const { return frames_[slot]; }

    /// Establish the persistent ghost set (a "reneighbor"-cadence operation).
    /// For each face, each subdomain selects its owned atoms within `cutoff` of
    /// that face boundary and sends them to the face-neighbor GPU, appended as
    /// ghosts in [n_local, n_local+n_ghost). Positions are raw global
    /// coordinates (the neighbor list min-images against the global box, so no
    /// PBC shift is needed). Ghost x/y/z/tag/type are transferred over NVLink
    /// via cudaMemcpyPeer. Self-neighbors (dims along an axis == 1) are skipped.
    /// Unlike a per-step refresh, build_halo also RECORDS, per transfer, the
    /// source owned-atom indices and the destination ghost offset (the
    /// `HaloPlan`) so subsequent steps can refresh ghost POSITIONS cheaply via
    /// `forward_comm()` without re-selecting the ghost set. Call this every K
    /// steps (when the neighbor list is rebuilt); call `forward_comm()` on the
    /// intervening steps. Communication is a 3-axis CASCADE (x, then y, then z):
    /// the y-pass forwards owned+x-ghost atoms and the z-pass forwards
    /// owned+x+y-ghosts, so edge/corner ghosts reach diagonal-neighbor GPUs.
    /// Within an axis the sendable set is snapshotted at the axis start, so
    /// ghosts received during that axis are not immediately re-forwarded.
    void build_halo(double cutoff);

    /// Backward-compatible alias for `build_halo` (full ghost-set rebuild).
    void exchange_halo(double cutoff) { build_halo(cutoff); }

    /// Refresh ghost POSITIONS only, reusing the `HaloPlan` recorded by the
    /// last `build_halo`. For each recorded transfer, gather the current owned
    /// positions of the previously-selected source atoms and cudaMemcpyPeer
    /// them into the same destination ghost slots. The ghost SET (which atoms,
    /// their tags/types, and the per-GPU counts) is unchanged — only x/y/z move.
    /// This is the per-step "forward communicate" of standard domain-decomposed
    /// MD; it is valid as long as no atom has moved more than the neighbor-list
    /// skin since `build_halo` (the reneighbor cadence guarantees that).
    void forward_comm();

    /// Merge all subdomains' owned atoms back into one cpu::Frame (positions /
    /// velocities by tag order matching the original global frame size). Used
    /// to validate partition round-trips and, later, full multi-GPU runs.
    void gather(cpu::Frame& out) const;

    /// Compute LJ 12-6 forces on each subdomain's OWNED atoms from its
    /// owned+ghost set (brute force, global-box min-image). Run after
    /// exchange_halo: each owned atom's force is complete (ghosts supply the
    /// cross-boundary neighbors), so this needs no reverse halo. Validates the
    /// decomposition+halo+force path against a single-GPU run.
    void compute_lj_bruteforce(double cutoff, double epsilon, double sigma);

    /// Gather owned-atom forces back to global tag/index order (like gather()
    /// but for fx/fy/fz). `out` sized to the global atom count.
    void gather_forces(cpu::Frame& out) const;

    /// Run `steps` of velocity-Verlet NVE across the GPUs (brute-force LJ):
    /// each step refreshes the halo, recomputes forces (owned, from
    /// owned+ghost), and integrates owned atoms. Acceleration is f/m (consistent
    /// units; this validates parity, not a specific unit system). Assumes atoms
    /// stay in their subdomain over the run (no migration yet).
    void run_nve(int steps, double dt, double cutoff, double epsilon,
                 double sigma);

    // ── Cluster neighbor-list force path (production performance) ──────────
    // The real per-subdomain force engine: each GPU runs the optimized
    // atom-centric cluster LJ kernel over its owned+ghost set, instead of the
    // O(n²) brute force above. The neighbor list is rebuilt only on reneighbor
    // steps (with a fresh halo over cutoff+skin); intervening steps just
    // forward_comm ghost positions and let the Verlet skin-gating reuse the CSR.

    /// Construct the per-subdomain force engine (Context + NeighborList +
    /// PairLJCut). Call once after `partition`. `skin` < 0 ⇒ auto (0.3*cutoff).
    /// The halo is built over `cutoff + skin` so ghosts cover the skin region
    /// and no cross-boundary neighbor appears mid-interval.
    void setup_force_engine(double cutoff, double epsilon, double sigma,
                            double skin = -1.0);

    /// Reneighbor step: migrate atoms that crossed subdomain boundaries to
    /// their new owner GPU, rebuild the halo (cutoff+skin), and rebuild the
    /// per-subdomain cluster neighbor list over the new owned+ghost set. Call
    /// every K steps.
    void reneighbor();

    /// Re-bucket owned atoms by their CURRENT owner subdomain: each atom that
    /// drifted out of its subdomain is gathered with full per-atom state on its
    /// old owner GPU and cudaMemcpyPeer'd to its new owner (single-process
    /// equivalent of MPI_Alltoallv). Ghosts are dropped (rebuilt by build_halo).
    /// Owned-atom global identity is carried by the persistent tag; gather()
    /// resolves by tag afterwards. No-op on a single GPU.
    void migrate();

    /// Compute LJ forces on every subdomain's owned atoms via the cluster
    /// kernel (zero forces over owned+ghost, refresh/launch the neighbor list,
    /// run PairLJCut). Assumes the halo/forward_comm for this step already ran.
    void compute_lj_nlist();

    /// Velocity-Verlet NVE using the cluster neighbor-list force path.
    /// `reneigh_every` steps between full halo+neighbor rebuilds; positions are
    /// forward-communicated on the intervening steps.
    void run_nve_nlist(int steps, double dt, double cutoff, double epsilon,
                       double sigma, double skin = -1.0, int reneigh_every = 10);

    /// One recorded P2P ghost transfer (one src subdomain → one dst face
    /// neighbor). `d_src_idx` (on the src device) lists the owned-atom indices
    /// selected at build_halo; their current positions are gathered into
    /// `d_sx/sy/sz` (also src-device) and copied into dst's ghost region at
    /// [base, base+cnt). Stable for the life of a reneighbor interval.
    struct HaloPlan {
        int src = 0;   // source subdomain slot
        int dst = 0;   // destination subdomain slot
        int axis = 0;  // cascade stage (0=x,1=y,2=z); replayed in order
        int cnt = 0;   // atoms transferred
        int base = 0;  // dst ghost offset (n_local + accumulated)
        int* d_src_idx = nullptr;          // [cnt] src-array indices (owned OR
                                           // prior-axis ghost) on src device
        real_t *d_sx = nullptr, *d_sy = nullptr, *d_sz = nullptr; // [cnt] send bufs
    };

  private:
    void destroy_frames();
    void destroy_halo_plan();
    void destroy_force_engine();

    DeviceManager dm_;
    std::unique_ptr<DomainDecomposition> dd_;
    std::vector<cuda::Frame> frames_; // per-slot device frame
    std::vector<int> n_local_;        // per-slot owned count
    std::vector<HaloPlan> halo_plan_; // recorded by build_halo, replayed by forward_comm
    int n_global_ = 0;
    double box_L_[3] = {0, 0, 0}; // global orthorhombic box edge lengths

    // Cluster force engine (one per subdomain); built by setup_force_engine.
    std::vector<cuda::Context> ctx_;
    std::vector<std::unique_ptr<cuda::NeighborList>> nlist_;
    std::vector<std::unique_ptr<cuda::PairLJCut>> pair_;
    std::vector<char> nlist_built_; // per-slot: NeighborList::build called (vector<bool> avoided)
    std::vector<int> nlist_cap_;    // per-slot: atom count the nlist buffers were built for
    bool engine_ready_ = false;
    double force_cutoff_ = 0, force_skin_ = 0;
};

} // namespace md
} // namespace cuda
} // namespace ATV
