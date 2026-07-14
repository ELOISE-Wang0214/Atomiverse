/// @file cpu/frame/column.h
/// @brief CPU column machinery for the frame-block registry: typed spans
///        (launch-boundary resolution), the type-erased ColumnRef the
///        structural encoders iterate, and registry-owned extension columns.
#pragma once

#include <cstdint>
#include <stdexcept>
#include <string>
#include <variant>

#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/precision.h"

namespace ATV {
namespace cpu {

/// @brief A resolved column: raw pointer + column count. The return type of
///        `Frame::span` — resolve ONCE at the launch boundary, loop over the
///        raw pointer (the 13.7× rule: never look a key up inside a loop).
/// @tparam T Element type.
template <typename T>
struct Span {
    T* data = nullptr; ///< Flat storage, row-major `[nmax, cols]`.
    int cols = 1;      ///< Columns per row (1 for plain 1D columns).

    /// @brief 1D convenience: decay to the raw pointer.
    operator T*() const noexcept { return data; }

    /// @brief 2D element access: row `i`, column `c`.
    T& operator()(int i, int c) const { return data[i * cols + c]; }
};

/// @brief Type-erased view of one registry column — what the structural
///        encoders (copy_row, pack/unpack_exchange) iterate. Carries the
///        schema facts it needs by value (no pointer into the schema table,
///        so registry-owned extras views are lifetime-safe). Values cross
///        through `double`, matching the exchange payload's cast-to-double
///        wire convention.
struct ColumnRef {
    const char* key = nullptr; ///< Column name (diagnostics).
    ColDType dtype = ColDType::Storage; ///< Element tier.
    ColGroup group = ColGroup::Core;    ///< Exchange payload section.
    int epi = 0;                        ///< Entries per group item.
    unsigned flags = kColNone;          ///< kCol* bits.
    void* data = nullptr;               ///< Flat storage `[nmax, cols]`.
    int cols = 0;                       ///< Runtime column count.

    /// @brief Read element (row, c) widened to double.
    double get(int row, int c) const {
        const std::size_t k =
            static_cast<std::size_t>(row) * cols + static_cast<std::size_t>(c);
        switch (dtype) {
        case ColDType::Storage: return static_cast<double>(static_cast<StorageReal*>(data)[k]);
        case ColDType::Acc:     return static_cast<double>(static_cast<AccReal*>(data)[k]);
        case ColDType::Tag:     return static_cast<double>(static_cast<TagInt*>(data)[k]);
        case ColDType::Int:     return static_cast<double>(static_cast<int*>(data)[k]);
        case ColDType::Int8:    return static_cast<double>(static_cast<std::int8_t*>(data)[k]);
        }
        return 0.0; // unreachable
    }

    /// @brief Write element (row, c) narrowed from double.
    void set(int row, int c, double v) const {
        const std::size_t k =
            static_cast<std::size_t>(row) * cols + static_cast<std::size_t>(c);
        switch (dtype) {
        case ColDType::Storage: static_cast<StorageReal*>(data)[k] = static_cast<StorageReal>(v); break;
        case ColDType::Acc:     static_cast<AccReal*>(data)[k] = static_cast<AccReal>(v); break;
        case ColDType::Tag:     static_cast<TagInt*>(data)[k] = static_cast<TagInt>(v); break;
        case ColDType::Int:     static_cast<int*>(data)[k] = static_cast<int>(v); break;
        case ColDType::Int8:    static_cast<std::int8_t*>(data)[k] = static_cast<std::int8_t>(v); break;
        }
    }

    /// @brief Copy element (row `src`, c) to (row `dst`, c) without a
    ///        double round-trip (bitwise row moves).
    void copy_element(int dst, int src, int c) const {
        const std::size_t kd =
            static_cast<std::size_t>(dst) * cols + static_cast<std::size_t>(c);
        const std::size_t ks =
            static_cast<std::size_t>(src) * cols + static_cast<std::size_t>(c);
        switch (dtype) {
        case ColDType::Storage: { auto* p = static_cast<StorageReal*>(data); p[kd] = p[ks]; break; }
        case ColDType::Acc:     { auto* p = static_cast<AccReal*>(data); p[kd] = p[ks]; break; }
        case ColDType::Tag:     { auto* p = static_cast<TagInt*>(data); p[kd] = p[ks]; break; }
        case ColDType::Int:     { auto* p = static_cast<int*>(data); p[kd] = p[ks]; break; }
        case ColDType::Int8:    { auto* p = static_cast<std::int8_t*>(data); p[kd] = p[ks]; break; }
        }
    }
};

/// @brief A registry-OWNED extension column (`Frame::insert`): a per-atom
///        payload some component or workflow attaches at runtime. Stored
///        rank-2 `[nmax, cols]` regardless of cols so growth is uniform.
struct ExtraColumn {
    std::string key;        ///< Column name (block is "atoms" in v1).
    ColDType dtype;         ///< Element tier.
    int cols = 1;           ///< Columns per row.
    bool migrates = true;   ///< Rides copy_row (row identity follows the atom).
    bool exchanged = false; ///< Rides pack/unpack_exchange (MPI + MC patch).

    /// One alternative per dtype tier; always `[nmax, cols]`.
    std::variant<xt::xtensor<StorageReal, 2>, xt::xtensor<AccReal, 2>,
                 xt::xtensor<TagInt, 2>, xt::xtensor<int, 2>,
                 xt::xtensor<std::int8_t, 2>>
        data;
};

/// @brief Capacity of the fixed ColumnRef scratch arrays the structural
///        encoders use (schema columns + extension-column headroom).
inline constexpr int kMaxExchangeColumns = 64;

/// @brief Thrown by `Frame::span` on a dtype-mismatched or unknown key.
class column_error : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

} // namespace cpu
} // namespace ATV
