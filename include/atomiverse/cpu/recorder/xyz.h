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

} // namespace cpu
} // namespace ATV
