/// @file cpu/io/frame_writer.cpp
/// @brief Per-frame trajectory writers — thin marshaling to molrs FFI.

#include "atomiverse/cpu/io/frame_writer.h"
#include "atomiverse/molrs.hpp"

#include <cstdint>

#include "readable_path.h"

namespace ATV {
namespace cpu {

namespace {
/// The cell matrix H as molrs takes it: all nine values, or none when the
/// engine has no cell (an all-zero H, a default-constructed Box). A non-zero
/// singular H is still sent — molrs refuses it rather than dropping the box.
::rust::Slice<const double> molrs_h(const std::array<double, 9>& box_h) {
    for (const double v : box_h)
        if (v != 0.0)
            return {box_h.data(), box_h.size()};
    return {box_h.data(), 0};
}
} // namespace

void write_xyz_frame(const std::string& path,
                     const std::vector<int>& atomic_number,
                     const std::vector<double>& x, const std::vector<double>& y,
                     const std::vector<double>& z,
                     const std::array<double, 9>& box_h, bool append) {
    // The bridge message does not name the file; with_molrs_context does.
    internal::with_molrs_context("cpu::write_xyz_frame", path, [&] {
        molrs::write_xyz(
            {path.data(), path.size()},
            {reinterpret_cast<const std::int32_t*>(atomic_number.data()),
             atomic_number.size()},
            {x.data(), x.size()}, {y.data(), y.size()}, {z.data(), z.size()},
            molrs_h(box_h), rust::Vec<molrs::KeyedMetaValue>{}, append);
    });
}

void write_mrec_frame(const std::string& path,
                      const std::vector<int>& atomic_number,
                      const std::vector<double>& x,
                      const std::vector<double>& y,
                      const std::vector<double>& z,
                      const std::array<double, 9>& box_h,
                      const std::vector<std::string>& field_names,
                      const std::vector<double>& field_data) {
    rust::Vec<rust::String> names;
    for (const auto& n : field_names)
        names.push_back(rust::String(n));
    internal::with_molrs_context("cpu::write_mrec_frame", path, [&] {
        molrs::write_mrec_frame(
            {path.data(), path.size()},
            {reinterpret_cast<const std::int32_t*>(atomic_number.data()),
             atomic_number.size()},
            {x.data(), x.size()}, {y.data(), y.size()}, {z.data(), z.size()},
            molrs_h(box_h), std::move(names),
            {field_data.data(), field_data.size()});
    });
}

} // namespace cpu
} // namespace ATV
