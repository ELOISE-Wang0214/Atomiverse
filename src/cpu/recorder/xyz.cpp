/// @file cpu/recorder/xyz.cpp
/// @brief XYZ writer + Zarr reader — delegating to molrs FFI, one frame at a time.

#include "atomiverse/cpu/recorder/xyz.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/io/frame_writer.h"
#include "atomiverse/molrs.hpp"

#include <array>
#include <stdexcept>
#include <string>
#include <vector>

namespace ATV {
namespace cpu {

namespace {
/// Extract one frame's raw arrays (type_id, x, y, z, box) for the molrs writer.
void extract(const Frame& frame, std::vector<int>& type_id,
             std::vector<double>& x, std::vector<double>& y,
             std::vector<double>& z, std::array<double, 9>& box) {
    const int n = static_cast<int>(frame.atoms.size());
    type_id.resize(n);
    x.resize(n);
    y.resize(n);
    z.resize(n);
    for (int i = 0; i < n; ++i) {
        type_id[i] = static_cast<int>(frame.atoms.type[i]);
        x[i] = static_cast<double>(frame.atoms.x[i]);
        y[i] = static_cast<double>(frame.atoms.y[i]);
        z[i] = static_cast<double>(frame.atoms.z[i]);
    }
    const auto& cell = frame.box.cell();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            box[static_cast<std::size_t>(i * 3 + j)] =
                static_cast<double>(cell(i, j));
}
} // namespace

// ── XYZ ──────────────────────────────────────────────────────────────────────

XYZ::XYZ(const std::string& filename) : filename_(filename) {}

void XYZ::write(const Frame& frame) { write_impl(frame, /*append=*/false, {}); }
void XYZ::append(const Frame& frame) { write_impl(frame, /*append=*/true, {}); }

void XYZ::write(const Frame& frame, const XYZMeta& meta) {
    write_impl(frame, /*append=*/false, meta);
}
void XYZ::append(const Frame& frame, const XYZMeta& meta) {
    write_impl(frame, /*append=*/true, meta);
}

void XYZ::write_impl(const Frame& frame, bool append, const XYZMeta& meta) {
    std::array<double, 9> box{};
    extract(frame, type_id_, x_, y_, z_, box);
    if (meta.empty()) {
        write_frame_xyz(filename_, type_id_, x_, y_, z_, box, append);
    } else {
        write_frame_xyz(filename_, type_id_, x_, y_, z_, box, append, meta);
    }
}

// ── ZarrReader ───────────────────────────────────────────────────────────────

ZarrReader::ZarrReader(const std::string& path) : path_(path) {}
ZarrReader::~ZarrReader() = default;

int ZarrReader::read_into(Frame& frame) {
    // Get the first frame of the Zarr store as a molrs FrameRef. Column
    // names follow the convention established by `dump_frame_to_zarr()`
    // in the polyethylene_tension bench: f64 columns for positions /
    // velocities / mass / tag / image components, f64 column for type.
    auto fref = molrs::read_frame_zarr_first({path_.data(), path_.size()});
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
    auto type = molrs::frame_column_f64(fr, "atoms", "type");
    auto box_h = molrs::frame_simbox(fr);

    const int n = static_cast<int>(x.size());
    if (n == 0)
        throw std::runtime_error("ZarrReader: empty Zarr at " + path_);
    if (n > frame.atoms.size())
        throw std::runtime_error("ZarrReader: n_global=" + std::to_string(n)
                                 + " > frame.atoms.size()="
                                 + std::to_string(frame.atoms.size()) + " at "
                                 + path_);
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
