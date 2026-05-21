/// @file core/types.h
/// @brief Core vocabulary types for the Atomiverse execution model.
/// CUDA-free header.
#pragma once

namespace ATV {

/// @brief Execution scheduling trait for a stage binding.
///
/// The Backend maps this trait to execution streams:
/// - **blocking** -- main stream (serialized with other blocking ops)
/// - **deferrable** -- secondary stream (may overlap with main stream work)
enum class Scheduling {
    blocking,  ///< Execute on main stream (serialized).
    deferrable ///< Execute on secondary stream (may overlap).
};

/// @brief Backpressure policy for the CPU analysis executor's snapshot ring.
///
/// Applied when a new snapshot is acquired while the ring is saturated:
/// - **drop_oldest** -- overwrite the oldest slot, increment drop counter.
/// Allocation-free hot path; recommended default for production throughput.
/// - **block** -- hot-path submit blocks until a worker frees a slot.
/// Allocation-free hot path; use when every snapshot must be delivered.
/// - **grow** -- **opt-in** policy that expands ring capacity on saturation.
/// Growth path allocates inside `AnalysisExecutor::acquire_slot()` and is
/// therefore a documented exception to the "no allocation in launch()" rule.
/// Use only in diagnostic / benchmark runs where neither drops nor stalls
/// are acceptable. Not recommended for production hot paths.
enum class BackpressurePolicy { drop_oldest, block, grow };

// Forward type declarations (actual types defined in their respective headers)
class Mat3;
class Vec3;

using Mat33 = Mat3;

} // namespace ATV
