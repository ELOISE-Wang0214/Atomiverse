/// @file cpu/io/frame_writer.h
/// @brief Per-frame trajectory writers — raw arrays in, delegating to molrs.
///
/// The engine hands molrs one frame's worth of raw per-atom arrays; molrs builds
/// the `Frame` and writes it (XYZ or `*.mrec`). No C++ record/accumulator type.
///
/// CUDA-safe: only `std::vector`/`std::array`/`std::string` (no `rust/cxx.h`).
#pragma once

#include <array>
#include <string>
#include <vector>

namespace ATV {
namespace cpu {

/// @brief Write one frame to a standard-XYZ file (element + coordinates).
///
/// Delegates to molrs `write_xyz`. @param atomic_number per-atom
/// atomic number Z (element symbol derived in molrs). @param x,y,z per-atom
/// coordinates [Å]. @param box_h row-major 3×3 cell matrix H [Å]; all zeros
/// means no cell, any other singular H is an error. @param append `false` truncates/creates the file, `true`
/// appends the frame.
void write_xyz_frame(const std::string& path,
                     const std::vector<int>& atomic_number,
                     const std::vector<double>& x, const std::vector<double>& y,
                     const std::vector<double>& z,
                     const std::array<double, 9>& box_h, bool append);

/// @brief Write one frame plus named per-atom fields as a `*.mrec` record
///        whose `frame` section is that frame (molrs' Zarr-v3 scientific
///        record).
///
/// Delegates to molrs `write_mrec_frame`; read it back with molrs
/// `read_mrec_frame`. The checkpoint door: one frame, whole, replaced on every
/// call — for a streaming trajectory use cpu::MrecRecorder.
/// @param path store path; must not use the retired `.zarr` suffix (molrs
/// refuses it), conventionally `*.mrec`. @param field_names one name per
/// extra per-atom column. @param field_data row-major `[n_fields, n_atoms]`,
/// matching `field_names` order (a length other than n_fields × n_atoms is
/// an error). @param box_h row-major 3×3 cell matrix H [Å]; all zeros means
/// no cell, any other singular H is an error.
void write_mrec_frame(const std::string& path,
                      const std::vector<int>& atomic_number,
                      const std::vector<double>& x,
                      const std::vector<double>& y,
                      const std::vector<double>& z,
                      const std::array<double, 9>& box_h,
                      const std::vector<std::string>& field_names,
                      const std::vector<double>& field_data);

} // namespace cpu
} // namespace ATV
