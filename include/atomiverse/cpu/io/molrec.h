/// @file cpu/io/molrec.h
/// @brief Trajectory container backed by molrs (PIMPL).
///
/// No MolFrame, no C++ mirror — all I/O state lives in molrs.
/// CUDA-safe: no rust/cxx.h in this header.
///
/// @see cpu::XYZ (recorder/xyz.h)
#pragma once

#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

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

    /// @brief Stage a frame directly from raw arrays (no cpu::Frame needed).
    ///
    /// For async analysis observers that consume a `FrameSnapshot`. `type_id`
    /// may be dummy (e.g. all 1) for element-agnostic analyses such as MSD.
    /// @param type_id per-atom atomic number Z (length N).
    /// @param x,y,z   per-atom coordinates [Angstrom] (length N).
    /// @param box_mat row-major 3x3 cell matrix [Angstrom].
    void stage_raw(const std::vector<std::int32_t>& type_id,
                   const std::vector<double>& x, const std::vector<double>& y,
                   const std::vector<double>& z,
                   const std::array<double, 9>& box_mat);

    /// @brief Attach a per-atom double field.
    void add_field(const std::string& name, std::vector<double> values);

    /// @brief Attach a frame-level scalar (always f64 across FFI).
    void add_scalar(const std::string& name, double value);

    /// @brief Attach a frame-level string.
    void add_string(const std::string& name, std::string value);

    /// @brief Commit the staging frame into the trajectory.
    void commit_frame();

    // ── Access ────────────────────────────────────────────────────────────────

    /// @brief Number of committed frames.
    int n_frames() const;

    /// @brief Discard all frames and reset.
    void clear();

    // ── Trajectory analysis (delegates to molrs `compute` kernels) ────────────

    /// @brief Mean squared displacement over the committed trajectory.
    ///
    /// Delegates to molrs `compute::MSD` (Direct mode: `⟨|r(t) − r(0)|²⟩` with
    /// the first frame as reference — LAMMPS `compute msd` semantics). No
    /// analysis math runs in C++. Returns the mean-MSD curve, one value per
    /// committed frame (index 0 is the reference, always 0); empty when fewer
    /// than two frames are committed.
    std::vector<double> compute_msd() const;

    /// @brief Einstein self-diffusion coefficient D [length²/time].
    ///
    /// Delegates to molrs `EinsteinDiffusion` (windowed MSD) + `LinearFit`:
    /// `D = slope / (2·dims)`, slope from an OLS fit of MSD(t) over the lag
    /// window `[fit_lo, fit_hi]` (fractions of the last lag). No analysis math
    /// in C++.
    /// @param dt     time between committed frames.
    /// @param dims   spatial dimensions (default 3).
    /// @param fit_lo,fit_hi  fit window fractions, `0 ≤ lo < hi ≤ 1`.
    /// @return D, or NaN on fewer than two frames / a fit error.
    double compute_diffusion(double dt, int dims = 3, double fit_lo = 0.2,
                             double fit_hi = 0.8) const;

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
