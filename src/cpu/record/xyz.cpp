/// @file cpu/record/xyz.cpp
/// @brief XYZ writer — delegating to molrs FFI, one frame at a time.

#include "atomiverse/cpu/record/xyz.h"
#include "../io/molrs_meta.h"
#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/molrs.hpp"

#include <array>
#include <string>
#include <vector>

namespace ATV {
namespace cpu {

namespace {
/// @brief Stopgap launch-boundary span resolution for a `const cpu::Frame&`.
///
/// Extract one frame's raw arrays (atomic number, x, y, z, cell matrix H) for
/// the molrs writer.
void extract(const Frame& frame, std::vector<int>& type_id,
             std::vector<double>& x, std::vector<double>& y,
             std::vector<double>& z, std::array<double, 9>& box) {
    const int n = static_cast<int>(frame.nlocal());
    type_id.resize(n);
    x.resize(n);
    y.resize(n);
    z.resize(n);
    // Hoisted once before the loop -- the 13.7x rule.
    const int* __restrict sp_type = frame.span(kAtom_type);
    const StorageReal* __restrict sp_x = frame.span(kAtom_x);
    const StorageReal* __restrict sp_y = frame.span(kAtom_y);
    const StorageReal* __restrict sp_z = frame.span(kAtom_z);
    for (int i = 0; i < n; ++i) {
        type_id[i] = sp_type[i];
        x[i] = static_cast<double>(sp_x[i]);
        y[i] = static_cast<double>(sp_y[i]);
        z[i] = static_cast<double>(sp_z[i]);
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

void XYZ::write(const Frame& frame) { write_impl(frame, /*append=*/false); }
void XYZ::append(const Frame& frame) { write_impl(frame, /*append=*/true); }

void XYZ::write_impl(const Frame& frame, bool append) {
    std::array<double, 9> box{};
    extract(frame, type_id_, x_, y_, z_, box);
    molrs::write_xyz(
        {filename_.data(), filename_.size()},
        {reinterpret_cast<const std::int32_t*>(type_id_.data()),
         type_id_.size()},
        {x_.data(), x_.size()}, {y_.data(), y_.size()}, {z_.data(), z_.size()},
        // No periodic cell (det == 0) crosses as no box, as in to_molrs.
        {box.data(), frame.box.volume() > 0.0 ? box.size() : 0},
        internal::encode_molrs_meta(frame.meta), append);
}

} // namespace cpu
} // namespace ATV
