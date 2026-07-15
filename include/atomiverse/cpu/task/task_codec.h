/// @file cpu/task/task_codec.h
/// @brief MPI-free byte-buffer codec for TaskSpec
///        (mpi-launcher-embedded-python-02).
///
/// Wire format (native endianness -- every rank of one MPI job runs the same
/// binary, so cross-endian portability is a non-goal by contract):
///
///   magic 'ATVT' (4 bytes)
///   u32  version            (== TaskSpec::kWireVersion)
///   u8   task_kind          (TaskKind)
///   u32  n_entries
///   entry * n_entries:
///     u16  component_id     (ComponentParams variant index)
///     u8   has_binding
///     [u32 phase_slot, i32 order]   iff has_binding
///     per-field payload     (vectors u64-length-prefixed)
///   i64  n_steps
///
/// Phase slots cross the wire via `ATV::slot_of(phase)` -- the one legal
/// enum->slot conversion -- and are range-checked against
/// `ATV::slot_count<md::Phase>()` on decode before the enum is rebuilt.
///
/// Decode rebuilds the TaskSpec through its validating constructors, so
/// task-kind invariant violations on the wire (a SinglePoint buffer carrying
/// more than one entry, a non-AM1 entry, or an MD buffer carrying an AM1
/// entry) throw rather than producing a half-valid spec. Unknown version /
/// task_kind / component_id and truncated buffers throw
/// `std::runtime_error` naming the offending value.
///
/// Free-function pair by design (the `comm::pack_exchange` precedent):
/// serialization is a boundary operation, not a component factory.
///
/// @see cpu/task/task_spec.h
#pragma once

#include <cstddef>
#include <vector>

#include "atomiverse/cpu/task/task_spec.h"

namespace ATV {
namespace cpu {
namespace task {

/// @brief Serialize a TaskSpec into a self-contained byte buffer.
/// @param spec The validated task description.
/// @return Wire bytes (see the format above).
std::vector<std::byte> pack(const TaskSpec& spec);

/// @brief Reconstruct a TaskSpec from wire bytes.
/// @param data Buffer produced by `pack` (same binary, same registry).
/// @param size Buffer length in bytes.
/// @return The decoded, re-validated TaskSpec.
/// @throws std::runtime_error on bad magic, version mismatch, unknown
///         task_kind / component_id, out-of-range phase slot, truncation,
///         or a wire-level task-kind invariant violation.
/// @throws std::invalid_argument when the re-validating TaskSpec
///         constructor rejects the decoded content.
TaskSpec unpack(const std::byte* data, std::size_t size);

} // namespace task
} // namespace cpu
} // namespace ATV
