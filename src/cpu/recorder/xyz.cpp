/// @file cpu/recorder/xyz.cpp
/// @brief XYZ, extxyz, and Zarr writers — all delegating to molrs FFI.

#include "atomiverse/cpu/recorder/xyz.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/molrs.hpp"

#include <stdexcept>

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

// ── ZarrReader ───────────────────────────────────────────────────────────────

ZarrReader::ZarrReader(const std::string& path) : path_(path) {}
ZarrReader::~ZarrReader() = default;

int ZarrReader::read_into(Frame& frame) {
    // Get the first frame of the Zarr store as a molrs FrameRef. Column
    // names follow the convention established by `dump_frame_to_zarr()`
    // in the polyethylene_tension bench: f64 columns for positions /
    // velocities / mass / tag / image components, i32 column for type.
    auto fref =
        molrs::molrec_read_zarr_first_frame({path_.data(), path_.size()});
    auto& fr = *fref;

    auto x = molrs::frame_column_f64(fr, "atoms", "x");
    auto y = molrs::frame_column_f64(fr, "atoms", "y");
    auto z = molrs::frame_column_f64(fr, "atoms", "z");
    auto vx = molrs::frame_column_f64(fr, "atoms", "vx");
    auto vy = molrs::frame_column_f64(fr, "atoms", "vy");
    auto vz = molrs::frame_column_f64(fr, "atoms", "vz");
    auto mass = molrs::frame_column_f64(fr, "atoms", "mass");
    auto tag = molrs::frame_column_f64(fr, "atoms", "tag");
    auto img_x = molrs::frame_column_f64(fr, "atoms", "image_x");
    auto img_y = molrs::frame_column_f64(fr, "atoms", "image_y");
    auto img_z = molrs::frame_column_f64(fr, "atoms", "image_z");
    // `type` is stored as an f64 add_field (matches dump_frame_to_zarr).
    // Set_geometry's "element" string column is also present but lossy
    // for arbitrary integer types; the explicit "type" field is the
    // authoritative source.
    auto type = molrs::frame_column_f64(fr, "atoms", "type");
    auto box_h = molrs::frame_simbox(fr);

    const int n = static_cast<int>(x.size());
    if (n == 0)
        throw std::runtime_error("ZarrReader: empty Zarr at " + path_);
    if (n > frame.atoms.size())
        throw std::runtime_error(
            "ZarrReader: n_global=" + std::to_string(n)
            + " > frame.atoms.size()=" + std::to_string(frame.atoms.size())
            + " at " + path_);
    if (box_h.size() != 9)
        throw std::runtime_error("ZarrReader: simbox.size != 9 at " + path_);

    ATV::Mat3 cell;
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c)
            cell(r, c) = box_h[r * 3 + c];
    frame.box.set_cell(cell);

    for (int i = 0; i < n; ++i) {
        frame.atoms.x[i] = x[i];
        frame.atoms.y[i] = y[i];
        frame.atoms.z[i] = z[i];
        frame.atoms.vx[i] = vx[i];
        frame.atoms.vy[i] = vy[i];
        frame.atoms.vz[i] = vz[i];
        frame.atoms.mass[i] = mass[i];
        frame.atoms.type[i] = static_cast<int>(type[i]);
        frame.atoms.tag[i] = static_cast<std::int64_t>(tag[i]);
        frame.atoms.image(i, 0) = static_cast<std::int8_t>(img_x[i]);
        frame.atoms.image(i, 1) = static_cast<std::int8_t>(img_y[i]);
        frame.atoms.image(i, 2) = static_cast<std::int8_t>(img_z[i]);
    }
    return n;
}

} // namespace cpu
} // namespace ATV
