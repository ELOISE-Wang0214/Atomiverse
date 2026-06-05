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
    /// The device frame for `slot` (owned in [0,n_local), ghosts later).
    const cuda::Frame& frame(int slot) const { return frames_[slot]; }

    /// Merge all subdomains' owned atoms back into one cpu::Frame (positions /
    /// velocities by tag order matching the original global frame size). Used
    /// to validate partition round-trips and, later, full multi-GPU runs.
    void gather(cpu::Frame& out) const;

  private:
    void destroy_frames();

    DeviceManager dm_;
    std::unique_ptr<DomainDecomposition> dd_;
    std::vector<cuda::Frame> frames_; // per-slot device frame
    std::vector<int> n_local_;        // per-slot owned count
    std::vector<std::vector<int>>
        owned_idx_;   // per-slot: global index of each owned atom (partition order)
    int n_global_ = 0;
};

} // namespace md
} // namespace cuda
} // namespace ATV
