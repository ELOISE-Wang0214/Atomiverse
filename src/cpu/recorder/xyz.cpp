/// @file cpu/recorder/xyz.cpp
/// @brief XYZ, extxyz, and Zarr writers — all delegating to molrs FFI.

#include "atomiverse/cpu/recorder/xyz.h"
#include "atomiverse/molrs.hpp"

namespace ATV {
namespace cpu {

// ── XYZ ──────────────────────────────────────────────────────────────────────

XYZ::XYZ(const std::string& filename) : filename_(filename) {}
XYZ::~XYZ() = default;

void XYZ::write(const MolRec& rec) {
    auto& h =
        *static_cast<molrs::AtvMolRec*>(const_cast<MolRec&>(rec).ffi_handle());
    molrs::xyz_write({filename_.data(), filename_.size()}, h);
}
void XYZ::write_ext(const MolRec& rec) {
    auto& h =
        *static_cast<molrs::AtvMolRec*>(const_cast<MolRec&>(rec).ffi_handle());
    molrs::xyz_write_ext({filename_.data(), filename_.size()}, h);
}
void XYZ::append(const MolRec& rec) {
    auto& h =
        *static_cast<molrs::AtvMolRec*>(const_cast<MolRec&>(rec).ffi_handle());
    molrs::xyz_append({filename_.data(), filename_.size()}, h);
}
void XYZ::append_ext(const MolRec& rec) {
    auto& h =
        *static_cast<molrs::AtvMolRec*>(const_cast<MolRec&>(rec).ffi_handle());
    molrs::xyz_append_ext({filename_.data(), filename_.size()}, h);
}

// ── ZarrWriter ───────────────────────────────────────────────────────────────

ZarrWriter::ZarrWriter(const std::string& path) : path_(path) {}
ZarrWriter::~ZarrWriter() = default;

void ZarrWriter::write(const MolRec& rec) {
    auto& h =
        *static_cast<molrs::AtvMolRec*>(const_cast<MolRec&>(rec).ffi_handle());
    molrs::molrec_write_zarr({path_.data(), path_.size()}, h);
}

} // namespace cpu
} // namespace ATV
