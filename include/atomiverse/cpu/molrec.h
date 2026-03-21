/// @file cpu/molrec.h
/// @brief Trajectory container backed by molrs (PIMPL).
///
/// No MolFrame, no C++ mirror — all I/O state lives in molrs.
/// CUDA-safe: no rust/cxx.h in this header.
///
/// @see cpu::XYZ (recorder/xyz.h)
#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atomiverse/core/int.h"
#include "atomiverse/core/real.h"

namespace ATV {
namespace cpu {

class Frame; // forward decl, defined in cpu/frame.h

/// @brief Molecular trajectory container backed by molrs (PIMPL).
///
/// All frame data lives in molrs::MolRec on the Rust side.
/// C++ only stages data and commits — no local frame storage.
class MolRec {
public:
    MolRec();
    ~MolRec();
    MolRec(MolRec&&) noexcept;
    MolRec& operator=(MolRec&&) noexcept;
    MolRec(const MolRec&) = delete;
    MolRec& operator=(const MolRec&) = delete;

    // ── Staging API ───────────────────────────────────────────────────────────

    /// @brief Copy atom geometry from a cpu::Frame into the staging area.
    void set_geometry(const Frame& frame);

    /// @brief Attach a per-atom Real field.
    void add_field(const std::string& name, std::vector<Real> values);

    /// @brief Attach a frame-level scalar (always f64 across FFI).
    void add_scalar(const std::string& name, Real value);

    /// @brief Attach a frame-level string.
    void add_string(const std::string& name, std::string value);

    /// @brief Commit the staging frame into the trajectory.
    void commit_frame();

    // ── Access ────────────────────────────────────────────────────────────────

    /// @brief Number of committed frames.
    Int n_frames() const;

    /// @brief Discard all frames and reset.
    void clear();

    /// @brief Opaque FFI handle (void* to avoid CXX leak).
    void* ffi_handle();
    const void* ffi_handle() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// ── Free functions ────────────────────────────────────────────────────────────

/// @brief Print a human-readable summary of the latest frame to stdout.
void print_summary(const MolRec& rec);

} // namespace cpu
} // namespace ATV
