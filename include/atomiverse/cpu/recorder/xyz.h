/// @file cpu/recorder/xyz.h
/// @brief XYZ, extxyz, and Zarr trajectory writers — all backed by molrs.
///
/// @see cpu::MolRec (cpu/molrec.h)
#pragma once

#include <string>

#include "atomiverse/cpu/molrec.h"

namespace ATV {
namespace cpu {

/// @brief XYZ and extxyz trajectory writer.
class XYZ {
  public:
    explicit XYZ(const std::string& filename);
    ~XYZ();
    XYZ(const XYZ&) = delete;
    XYZ& operator=(const XYZ&) = delete;

    void write(const MolRec& rec);
    void write_ext(const MolRec& rec);
    void append(const MolRec& rec);
    void append_ext(const MolRec& rec);

  private:
    std::string filename_;
};

/// @brief Zarr V3 trajectory writer backed by molrs::MolRec::write_zarr.
class ZarrWriter {
  public:
    explicit ZarrWriter(const std::string& path);
    ~ZarrWriter();
    ZarrWriter(const ZarrWriter&) = delete;
    ZarrWriter& operator=(const ZarrWriter&) = delete;

    /// @brief Write all committed frames to the Zarr store.
    void write(const MolRec& rec);

  private:
    std::string path_;
};

class Frame; // fwd

/// @brief Zarr V3 reader — loads the first frame of a Zarr-stored
///        trajectory into a cpu::Frame.
///
/// Companion to `ZarrWriter`. Used by checkpoint reload: stage 1 dumps
/// its end-state via `ZarrWriter`, subsequent debug iterations call
/// `ZarrReader::read_into(frame)` to skip stage 1.
///
/// The reader expects per-atom columns `x, y, z, vx, vy, vz, mass, tag,
/// image_x, image_y, image_z` (all f64) and an integer `type` column
/// in the `atoms` block, plus a 9-element `simbox` (row-major
/// h-matrix). These are exactly the fields `dump_frame_to_zarr()`
/// writes — see benchmarks/md/polyethylene_tension.cpp.
class ZarrReader {
  public:
    explicit ZarrReader(const std::string& path);
    ~ZarrReader();
    ZarrReader(const ZarrReader&) = delete;
    ZarrReader& operator=(const ZarrReader&) = delete;

    /// @brief Read the first frame of the Zarr store into `frame`.
    /// Overwrites `frame.atoms[0..n_global)` and `frame.box`. Throws on
    /// I/O error or if `n_global > frame.atoms.size()`.
    /// @return Number of atoms loaded (== `n_global` from the file).
    int read_into(Frame& frame);

  private:
    std::string path_;
};

} // namespace cpu
} // namespace ATV
