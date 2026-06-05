#pragma once

/// @file cuda/md/multi_gpu.cuh
/// @brief MultiGpu — single-process multi-GPU MD orchestrator. Owns the
///        DeviceManager (P2P over NVLink) and a 3D DomainDecomposition, splits
///        a global cpu::Frame into one device-resident subdomain frame per GPU,
///        and (in later phases) drives the per-step P2P halo exchange,
///        migration, and cross-GPU reduction.
///
/// This first cut implements partition (global frame -> N owned subdomain
/// frames, materialized on their devices) and gather (merge the subdomain
/// frames back into a global cpu::Frame, by tag) so the decomposition can be
/// validated against the single-GPU result before the comm layer lands.

#include <memory>
#include <vector>

#include "atomiverse/cuda/device_manager.cuh"
#include "atomiverse/cuda/domain_decomposition.cuh"
#include "atomiverse/cuda/frame.cuh"

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
    /// intervening steps. This single pass is complete for a 1D decomposition;
    /// the 3-axis cascade for edge/corner ghosts is `build_halo`'s next step.
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

    /// One recorded P2P ghost transfer (one src subdomain → one dst face
    /// neighbor). `d_src_idx` (on the src device) lists the owned-atom indices
    /// selected at build_halo; their current positions are gathered into
    /// `d_sx/sy/sz` (also src-device) and copied into dst's ghost region at
    /// [base, base+cnt). Stable for the life of a reneighbor interval.
    struct HaloPlan {
        int src = 0;   // source subdomain slot
        int dst = 0;   // destination subdomain slot
        int cnt = 0;   // atoms transferred
        int base = 0;  // dst ghost offset (n_local + accumulated)
        int* d_src_idx = nullptr;          // [cnt] owned indices on src device
        real_t *d_sx = nullptr, *d_sy = nullptr, *d_sz = nullptr; // [cnt] send bufs
    };

  private:
    void destroy_frames();
    void destroy_halo_plan();

    DeviceManager dm_;
    std::unique_ptr<DomainDecomposition> dd_;
    std::vector<cuda::Frame> frames_; // per-slot device frame
    std::vector<int> n_local_;        // per-slot owned count
    std::vector<std::vector<int>>
        owned_idx_;   // per-slot: global index of each owned atom (partition order)
    std::vector<HaloPlan> halo_plan_; // recorded by build_halo, replayed by forward_comm
    int n_global_ = 0;
    double box_L_[3] = {0, 0, 0}; // global orthorhombic box edge lengths
};

} // namespace md
} // namespace cuda
} // namespace ATV
