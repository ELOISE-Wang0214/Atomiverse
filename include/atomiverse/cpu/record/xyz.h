/// @file cpu/record/xyz.h
/// @brief XYZ trajectory writer — backed by molrs.
///
/// The writer takes a `cpu::Frame` and hands its raw per-atom arrays and
/// metadata to molrs one frame at a time (molrs `write_xyz`); there is
/// no C++ record type.
///
/// @see cpu::write_xyz_frame / cpu::write_mrec_frame (cpu/io/frame_writer.h)
#pragma once

#include <string>
#include <vector>

namespace ATV {
namespace cpu {

class Frame; // fwd — CUDA-safe (this header is pulled by 02_rhf.cu)

/// @brief XYZ trajectory writer (standard element + coordinates).
///
/// `write`/`append` write one frame's geometry through molrs. Extraction
/// scratch is reused across calls, so per-frame appends on a long-lived
/// writer do not reallocate on the steady-state path.
class XYZ {
  public:
    /// @brief Construct.
    /// @param filename Output XYZ path.
    explicit XYZ(const std::string& filename);

    void write(const Frame& frame);  ///< truncate/create, write one frame
    void append(const Frame& frame); ///< append one frame

  private:
    void write_impl(const Frame& frame, bool append);

    std::string filename_;
    std::vector<int> type_id_; ///< per-frame extraction scratch (reused)
    std::vector<double> x_, y_, z_;
};

} // namespace cpu
} // namespace ATV
