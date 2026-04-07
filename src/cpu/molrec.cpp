/// @file cpu/molrec.cpp
/// @brief MolRec — pure PIMPL forwarding to molrs.
///
/// Bridge always uses f64 (molrs-core hardcodes F = f64).
/// Zero-copy when ATV_REAL=double; promotion when ATV_REAL=float.

#include "atomiverse/cpu/molrec.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/molrs.hpp"

#include <cstdint>
#include <type_traits>
#include <vector>

namespace ATV {
namespace cpu {

// ── Promotion helper ─────────────────────────────────────────────────────────
// Zero-copy when T == double; promotes to double otherwise.

template<typename T>
static rust::Slice<const double> as_f64(const T* data, size_t n,
                                         std::vector<double>& buf)
{
    if constexpr (std::is_same_v<T, double>) {
        return {data, n};
    } else {
        buf.resize(n);
        for (size_t i = 0; i < n; ++i) buf[i] = static_cast<double>(data[i]);
        return {buf.data(), n};
    }
}

// ── PIMPL ────────────────────────────────────────────────────────────────────

struct MolRec::Impl {
    rust::Box<molrs::AtvMolRec> inner;
    Impl() : inner(molrs::molrec_new()) {}
};

// ── Lifecycle ────────────────────────────────────────────────────────────────

MolRec::MolRec() : impl_(std::make_unique<Impl>()) {}
MolRec::~MolRec() = default;
MolRec::MolRec(MolRec&&) noexcept = default;
MolRec& MolRec::operator=(MolRec&&) noexcept = default;

void* MolRec::ffi_handle() { return &*impl_->inner; }
const void* MolRec::ffi_handle() const { return &*impl_->inner; }

// ── Staging API ──────────────────────────────────────────────────────────────

void MolRec::set_geometry(const Frame& frame) {
    const int n = static_cast<int>(frame.atoms.size());

    std::vector<int32_t> tid(n);
    for (int i = 0; i < n; ++i)
        tid[i] = static_cast<int32_t>(frame.atoms.type[i]);

    double box_mat[9];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            box_mat[i * 3 + j] = static_cast<double>(frame.box.data[i][j]);

    std::vector<double> xb, yb, zb; // unused when double == double
    molrs::molrec_set_geometry(
        *impl_->inner,
        {tid.data(), tid.size()},
        as_f64(frame.atoms.x.data(), n, xb),
        as_f64(frame.atoms.y.data(), n, yb),
        as_f64(frame.atoms.z.data(), n, zb),
        {box_mat, 9}
    );
}

void MolRec::add_field(const std::string& name, std::vector<double> values) {
    std::vector<double> buf;
    molrs::molrec_add_field(
        *impl_->inner, {name.data(), name.size()},
        as_f64(values.data(), values.size(), buf)
    );
}

void MolRec::add_scalar(const std::string& name, double value) {
    molrs::molrec_add_scalar(*impl_->inner, {name.data(), name.size()},
                             static_cast<double>(value));
}

void MolRec::add_string(const std::string& name, std::string value) {
    molrs::molrec_add_string(
        *impl_->inner, {name.data(), name.size()}, {value.data(), value.size()}
    );
}

void MolRec::commit_frame() {
    molrs::molrec_commit_frame(*impl_->inner);
}

int MolRec::n_frames() const {
    return static_cast<int>(molrs::molrec_n_frames(
        *static_cast<const molrs::AtvMolRec*>(ffi_handle())));
}

void MolRec::clear() {
    molrs::molrec_clear(*impl_->inner);
}

// ── print_summary ────────────────────────────────────────────────────────────

void print_summary(const MolRec& rec) {
    molrs::molrec_print_summary(
        *static_cast<const molrs::AtvMolRec*>(rec.ffi_handle()));
}

} // namespace cpu
} // namespace ATV
