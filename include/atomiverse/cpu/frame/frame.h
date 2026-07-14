/// @file cpu/frame/frame.h
/// @brief cpu::Frame -- Atom + Box + SystemScalars behind the molrs-shaped
///        mapping facade (frame-block-01).
///
/// The signed Frame verdict (docs/developer/frame-block-redesign.md): Frame
/// aligns with molrs — `{ blocks: string→Block, meta: string→MetaValue,
/// simbox }` — as an independent per-backend implementation. Step 01 keeps
/// the physical storage untouched (named members on `cpu::Atom`, plus `Box`
/// and `SystemScalars`) and fronts it with the mapping API:
///
///   - `span(kAtom_x)`      — typed column resolution, ONCE per launch
///   - `has("atoms","rho")` — presence check (schema, topology group, extras)
///   - `insert(...)`        — registry-owned extension columns
///   - `meta("nlocal")`     — the ONE typed-meta bookkeeping region (views
///                            over counts + scalars in this step; physical
///                            ownership moves in the terminal step)
///
/// `box` is the molrs `simbox` slot — outside the mapping on both sides of
/// the FFI, by molrs's own design. The `tag → row` map stays an engine-local
/// derived index (never serialized; not a meta value).
///
/// ## The 13.7× rule
///
/// Resolution stops at the launch boundary: resolve spans/meta once at the
/// top of `launch()`, loop over raw pointers. A key lookup inside a loop is
/// 13.7× slower (measured, docs/developer/frame-block-redesign.md §4). Debug
/// The driver opens a resolution window around every component launch and
/// throws past the budget — the discipline is machine-checked, not reviewed.
/// Counting is constant-time and sits only at resolution sites (which are
/// launch-boundary-rare by construction); the LJ ⚡ gate confirms zero cost.
///
/// @see cpu/frame/atom.h (the registry binding), core/frame_schema.h (THE
///      field list), cpu/frame/column.h (Span/ColumnRef/ExtraColumn)
#pragma once

#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <variant>

#include "atomiverse/core/system_scalars.h"
#include "atomiverse/cpu/frame/atom.h"
#include "atomiverse/cpu/frame/box.h"
#include "atomiverse/cpu/frame/column.h"

namespace ATV {
namespace cpu {

/// @brief A typed meta value — the extension of molrs's `string→string`
///        meta agreed in the signed verdict (revision 2026-07-14).
using MetaValue = std::variant<std::int64_t, double, std::string>;

// The bookkeeping VIEWS the meta region exposes in step 01. Ownership stays
// on the members (hot reads like `nlocal` keep their direct paths until the
// terminal step); meta is the one region that names them.
#define ATV_FRAME_META_INT_VIEWS(X)                                            \
    X(nlocal, atoms.nlocal)                                                    \
    X(nghost, atoms.nghost)                                                    \
    X(nmax, atoms.nmax)
#define ATV_FRAME_META_REAL_VIEWS(X)                                           \
    X(energy, scalars.energy)                                                  \
    X(energy_pair, scalars.energy_pair)                                        \
    X(energy_bond, scalars.energy_bond)                                        \
    X(energy_angle, scalars.energy_angle)                                      \
    X(energy_dihedral, scalars.energy_dihedral)                                \
    X(kinetic, scalars.kinetic)                                                \
    X(conserved, scalars.conserved)

/// @brief Default per-launch resolution budget the driver arms: generous
///        against legitimate launch-boundary resolution (every schema column
///        plus meta twice over), tiny against an in-loop lookup (N atoms).
inline constexpr int kResolveBudgetPerLaunch = 64;

/// @brief CPU backend execution frame: per-atom data behind the mapping
///        facade, cell (simbox), scalar results, and the typed meta region.
class Frame {
  public:
    /// @brief Construct with an initial per-atom array capacity.
    ///
    /// `atoms` is allocated with capacity `n_atoms` and `nlocal = n_atoms`
    /// (single-rank semantics: every row is locally owned). Multi-rank
    /// drivers reduce `nlocal` via the first `cpu::comm::Exchange`.
    /// @param n_atoms Initial atom count / capacity.
    explicit Frame(int n_atoms) : atoms(n_atoms) { atoms.nlocal = n_atoms; }

    /// @brief Per-atom container. Positions, velocities, forces, type,
    /// tag, image, per-atom bonded topology, and the `tag → local` map.
    cpu::Atom atoms;

    /// @brief Per-frame system scalar results (energy, kinetic, conserved,
    /// stress, kinetic stress). See `core/system_scalars.h`.
    SystemScalars scalars{};

    /// @brief Simulation cell (the molrs `simbox` slot).
    cpu::Box box{};

    // ── Mapping API (launch-boundary resolution) ────────────────────────────

    /// @brief Resolve a typed column to a raw span. Call ONCE at the top of
    ///        `launch()`; loop over the span, never over this call.
    /// @tparam T Element type carried by the key (checked against the bound
    ///         member's element type).
    /// @param key A schema key constant (`kAtom_x`, …) or a hand-built key
    ///        for an extension column.
    /// @return `{data, cols}` view of the live storage (invalidated by
    ///         grow/set_capacity, like every raw pointer).
    /// @throws column_error on an unknown key or an element-type mismatch.
    template <typename T>
    Span<T> span(const ColKey<T>& key) {
        note_resolution_();
        if (std::strcmp(key.block, "atoms") != 0)
            throw column_error(std::string("Frame::span: unknown block '")
                               + key.block + "'");
        Span<T> out{};
        bool found = false;
        atoms.for_each_column([&](const FieldSpec& spec, auto& arr) {
            using value_t = typename std::decay_t<decltype(arr)>::value_type;
            if (found || std::strcmp(spec.key, key.key) != 0)
                return;
            if constexpr (std::is_same_v<value_t, T>) {
                out = Span<T>{arr.data(), atoms.cols_of(spec.stride)};
                found = true;
            } else {
                throw column_error(
                    std::string("Frame::span: element-type mismatch on '")
                    + key.key + "'");
            }
        });
        if (!found) {
            for (auto& e : atoms.extras) {
                if (e.key != key.key)
                    continue;
                std::visit(
                    [&](auto& arr) {
                        using value_t =
                            typename std::decay_t<decltype(arr)>::value_type;
                        if constexpr (std::is_same_v<value_t, T>) {
                            out = Span<T>{arr.data(), e.cols};
                            found = true;
                        } else {
                            throw column_error(
                                std::string("Frame::span: element-type "
                                            "mismatch on extension column '")
                                + key.key + "'");
                        }
                    },
                    e.data);
                break;
            }
        }
        if (!found)
            throw column_error(std::string("Frame::span: no column '")
                               + key.key + "' in block 'atoms'");
        return out;
    }

    /// @brief Presence check. Schema columns are always present except the
    ///        lazy wide-topology group, which exists once the system carries
    ///        bonded topology; extension columns exist once inserted.
    /// @param block Block name ("atoms" in step 01).
    /// @param key Column name.
    bool has(const char* block, const char* key) const {
        note_resolution_();
        if (std::strcmp(block, "atoms") != 0)
            return false;
        bool present = false;
        bool seen = false;
        atoms.for_each_column([&](const FieldSpec& spec, auto& /*arr*/) {
            if (seen || std::strcmp(spec.key, key) != 0)
                return;
            seen = true;
            present = !(spec.flags & kColTopology) || atoms.has_bonded();
        });
        if (!seen)
            for (const auto& e : atoms.extras)
                if (e.key == key)
                    return true;
        return present;
    }

    /// @brief Insert a registry-owned per-atom extension column (block
    ///        "atoms"), zero-filled at the current capacity. It is sized
    ///        with the frame from then on and rides `copy_row` / the
    ///        exchange payload per its tags.
    /// @param key Column name; must not collide with schema or extras.
    /// @param dtype Element tier.
    /// @param cols Columns per row (>= 1).
    /// @param migrates Rides copy_row (default true).
    /// @param exchanged Rides pack/unpack_exchange (default false).
    /// @return The owned column (reference valid until the next insert).
    ExtraColumn& insert(const std::string& key, ColDType dtype, int cols = 1,
                        bool migrates = true, bool exchanged = false) {
        return atoms.add_column(key, dtype, cols, migrates, exchanged);
    }

    // ── Typed meta (the ONE bookkeeping region) ─────────────────────────────

    /// @brief Read a meta entry: a bookkeeping view (`nlocal`, `energy`, …)
    ///        or a user-owned entry.
    /// @throws std::out_of_range on an unknown key.
    MetaValue meta(const std::string& key) const {
        note_resolution_();
#define ATV_META_INT_GET(name, member)                                         \
    if (key == #name)                                                          \
        return MetaValue{static_cast<std::int64_t>(member)};
        ATV_FRAME_META_INT_VIEWS(ATV_META_INT_GET)
#undef ATV_META_INT_GET
#define ATV_META_REAL_GET(name, member)                                        \
    if (key == #name)                                                          \
        return MetaValue{static_cast<double>(member)};
        ATV_FRAME_META_REAL_VIEWS(ATV_META_REAL_GET)
#undef ATV_META_REAL_GET
        return owned_meta_.at(key);
    }

    /// @brief Write a meta entry. View keys route to the live member;
    ///        anything else lands in the owned map.
    /// @throws std::bad_variant_access if a view key is written with the
    ///         wrong alternative.
    void set_meta(const std::string& key, MetaValue v) {
#define ATV_META_INT_SET(name, member)                                         \
    if (key == #name) {                                                        \
        member = static_cast<int>(std::get<std::int64_t>(v));                  \
        return;                                                                \
    }
        ATV_FRAME_META_INT_VIEWS(ATV_META_INT_SET)
#undef ATV_META_INT_SET
#define ATV_META_REAL_SET(name, member)                                        \
    if (key == #name) {                                                        \
        member = std::get<double>(v);                                          \
        return;                                                                \
    }
        ATV_FRAME_META_REAL_VIEWS(ATV_META_REAL_SET)
#undef ATV_META_REAL_SET
        owned_meta_[key] = std::move(v);
    }

    /// @brief Whether a meta key exists (view or owned).
    bool has_meta(const std::string& key) const {
#define ATV_META_HAS(name, member)                                             \
    if (key == #name)                                                          \
        return true;
        ATV_FRAME_META_INT_VIEWS(ATV_META_HAS)
        ATV_FRAME_META_REAL_VIEWS(ATV_META_HAS)
#undef ATV_META_HAS
        return owned_meta_.count(key) != 0;
    }

    // ── Debug resolution accounting (the machine-checked discipline) ────────

    /// @brief Open a resolution window with a budget (driver-side, once per
    ///        component launch).
    /// @param budget Max key resolutions allowed inside the window.
    void begin_resolution_window(int budget) const {
        resolve_budget_ = budget;
        resolve_count_ = 0;
    }

    /// @brief Close the window; throws if the budget was exceeded — i.e.
    ///        some code resolved keys inside a loop.
    /// @param who Diagnostic label (component name).
    /// @throws std::logic_error on a budget violation.
    void end_resolution_window(const char* who) const {
        const int count = resolve_count_;
        const int budget = resolve_budget_;
        resolve_budget_ = 0;
        if (budget > 0 && count > budget)
            throw std::logic_error(
                std::string("Frame: '") + who + "' resolved "
                + std::to_string(count) + " keys in one launch (budget "
                + std::to_string(budget)
                + ") — key lookups belong at the launch boundary, never "
                  "inside a loop (the 13.7x rule)");
    }

    /// @brief Resolutions since the window opened (tests / diagnostics).
    int resolution_count() const { return resolve_count_; }

  private:
    /// @brief Count one key resolution.
    void note_resolution_() const { ++resolve_count_; }

    /// User-owned meta entries (non-view keys).
    std::unordered_map<std::string, MetaValue> owned_meta_;
    mutable int resolve_count_ = 0;  ///< Resolutions in the open window.
    mutable int resolve_budget_ = 0; ///< 0 = no window armed.
};

} // namespace cpu
} // namespace ATV
