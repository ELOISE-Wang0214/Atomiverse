/// @file cpu/md/mpi/halo_exchange.h
/// @brief Ghost-atom halo exchange components for MPI domain decomposition.
///
/// Phase 1: both classes are structural no-ops that register default bindings.
/// Phase 2 implements the real non-blocking comm patterns:
///   - `HaloExchangeForward::launch` issues `MPI_Isend/Irecv` on the 6 face
///     neighbors carrying owned positions into ghost rows.
///   - `HaloExchangeReverse::launch` reduces ghost-row force contributions
///     back onto the owner.
///
/// Both share a `HaloExchangeState` record (buffers, neighbor topology,
/// persistent request handles) held inside `DomainDecomposition`; the two
/// components are thin dispatchers over it.
///
/// @see .claude/specs/parallel-md-engine.md  § A.4
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cpu/component.h"
#include "atomiverse/core/binding.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace md {
namespace mpi {

/// @brief Forward pass: owned positions → ghost rows on face neighbors.
///
/// Phase 1: no-op. Phase 2: `MPI_Isend/Irecv` the six faces. In single-rank
/// mode (`runtime.mpi_size == 1`) the Phase 2 implementation will take an
/// early return, preserving single-rank throughput.
class HaloExchangeForward : public cpu::ComponentBase {
public:
    std::string name() const override { return "halo_forward"; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::pre_force), -20}};
    }
    void launch(const cpu::LaunchContext& ctx) override;
};

/// @brief Reverse pass: ghost-row forces → owner's force arrays.
///
/// Phase 1: no-op. Phase 2: reduce `fx/fy/fz` contributions on ghost rows
/// back onto the owning rank with `MPI_Isend/Irecv` plus local atomic adds.
class HaloExchangeReverse : public cpu::ComponentBase {
public:
    std::string name() const override { return "halo_reverse"; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::post_force), -1000}};
    }
    void launch(const cpu::LaunchContext& ctx) override;
};

} // namespace mpi
} // namespace md
} // namespace cpu
} // namespace ATV
