/// @file core/frame_schema.h
/// @brief THE single enumeration of per-atom fields — the frame-block
///        registry's schema (backend-agnostic).
///
/// This is step 1 of the signed Frame verdict (frame-block-01): one X-macro
/// list is the *only* place the per-atom field set is written down. Every
/// consumer that used to hand-copy the list — `Atom::grow`, `set_capacity`,
/// `copy_row`, `comm::pack_exchange` / `unpack_exchange` — derives from this
/// table instead. The grow data-destruction bug of probe 4 was a direct
/// product of the old "copy the list at every site" pattern; with one list
/// there is no second copy to fall out of sync.
///
/// ## The macro row
///
/// `X(member, key, Elem, DTag, Stride, Group, epi, flags)`
///
/// - `member` — the `cpu::Atom` data member the CPU backend binds.
/// - `key`    — the molrs exchange-schema column name (block `"atoms"`),
///              1:1 with the member name by the signed verdict.
/// - `Elem`   — element type from the ATV precision tiers
///              (`StorageReal` / `AccReal` / `TagInt` / `int` / `std::int8_t`).
/// - `DTag`   — the matching `ColDType` enumerator (runtime dispatch tag).
/// - `Stride` — `ColStride`: how many columns per row, and who owns the
///              runtime stride (`One`, `Fixed3`, or a per-atom topology
///              stride like `BondPA`).
/// - `Group`  — `ColGroup`: which counted section of the exchange payload
///              the column belongs to (`Core` is packed flat; `Bond` /
///              `Angle` / `Dihedral` items are length-prefixed by their
///              count column; `Special` is the nspecial-prefixed list).
/// - `epi`    — entries per group item (e.g. `angle_atom` stores 2 partner
///              tags per angle; 0 for non-group columns).
/// - `flags`  — `kColCount` (this column IS its group's length prefix),
///              `kColTopology` (lazy wide storage, only allocated once the
///              system has bonded topology), `kColIdentityTail` (new tail
///              rows are seeded with the row index — the tag invariant).
///
/// The order of rows is load-bearing: it IS the exchange wire order
/// (Core scalars, then bond / angle / dihedral counted groups, then the
/// special-neighbor section) — pinned by the MPI migration suite and the
/// MC transaction-rollback tests.
#pragma once

#include <cstdint>

#include "atomiverse/core/precision.h"

namespace ATV {

/// @brief Runtime dtype tag of a registry column (ATV precision tiers).
enum class ColDType : std::uint8_t { Storage, Acc, Tag, Int, Int8 };

/// @brief Column width semantics: fixed widths, or which runtime per-atom
///        topology stride sizes the column.
enum class ColStride : std::uint8_t {
    One,         ///< plain 1D column
    Fixed3,      ///< [nmax, 3] (image, nspecial)
    BondPA,      ///< [nmax, bond_per_atom]
    AnglePA2,    ///< [nmax, 2 * angle_per_atom]
    AnglePA,     ///< [nmax, angle_per_atom]
    DihedralPA3, ///< [nmax, 3 * dihedral_per_atom]
    DihedralPA,  ///< [nmax, dihedral_per_atom]
    MaxSpecial,  ///< [nmax, maxspecial]
};

/// @brief Which counted section of the exchange payload a column rides.
enum class ColGroup : std::uint8_t { Core, Bond, Angle, Dihedral, Special };

inline constexpr unsigned kColNone = 0;
inline constexpr unsigned kColCount = 1u << 0;    ///< group length prefix
inline constexpr unsigned kColTopology = 1u << 1; ///< lazy wide storage
inline constexpr unsigned kColIdentityTail = 1u << 2; ///< tag tail seeding
inline constexpr unsigned kColWriteBack = 1u << 3; ///< in the write_back contract

/// @brief One schema row — everything a registry consumer needs to know
///        about a column except the backend's storage object.
struct FieldSpec {
    const char* key;    ///< molrs exchange-schema column name
    ColDType dtype;     ///< element tier
    ColStride stride;   ///< width semantics
    ColGroup group;     ///< exchange payload section
    int epi;            ///< entries per group item (0 = not a group payload)
    unsigned flags;     ///< kCol* bits
};

// ── THE field list ──────────────────────────────────────────────────────────
//
//                member         key             Elem         DTag     Stride       Group     epi flags
#define ATV_ATOM_FIELDS(X)                                                                              \
    X(x,             "x",            StorageReal, Storage, One,         Core,     0, kColWriteBack)          \
    X(y,             "y",            StorageReal, Storage, One,         Core,     0, kColWriteBack)          \
    X(z,             "z",            StorageReal, Storage, One,         Core,     0, kColWriteBack)          \
    X(vx,            "vx",           StorageReal, Storage, One,         Core,     0, kColWriteBack)          \
    X(vy,            "vy",           StorageReal, Storage, One,         Core,     0, kColWriteBack)          \
    X(vz,            "vz",           StorageReal, Storage, One,         Core,     0, kColWriteBack)          \
    X(fx,            "fx",           AccReal,     Acc,     One,         Core,     0, kColWriteBack)          \
    X(fy,            "fy",           AccReal,     Acc,     One,         Core,     0, kColWriteBack)          \
    X(fz,            "fz",           AccReal,     Acc,     One,         Core,     0, kColWriteBack)          \
    X(mass,          "mass",         StorageReal, Storage, One,         Core,     0, kColNone)          \
    X(charge,        "charge",       StorageReal, Storage, One,         Core,     0, kColWriteBack)          \
    X(type,          "type",         int,         Int,     One,         Core,     0, kColNone)          \
    X(tag,           "tag",          TagInt,      Tag,     One,         Core,     0, kColIdentityTail)  \
    X(image,         "image",        std::int8_t, Int8,    Fixed3,      Core,     0, kColNone)          \
    X(num_bond,      "num_bond",     int,         Int,     One,         Bond,     0, kColCount)         \
    X(bond_atom,     "bond_atom",    TagInt,      Tag,     BondPA,      Bond,     1, kColTopology)      \
    X(bond_type,     "bond_type",    int,         Int,     BondPA,      Bond,     1, kColTopology)      \
    X(num_angle,     "num_angle",    int,         Int,     One,         Angle,    0, kColCount)         \
    X(angle_atom,    "angle_atom",   TagInt,      Tag,     AnglePA2,    Angle,    2, kColTopology)      \
    X(angle_type,    "angle_type",   int,         Int,     AnglePA,     Angle,    1, kColTopology)      \
    X(num_dihedral,  "num_dihedral", int,         Int,     One,         Dihedral, 0, kColCount)         \
    X(dihedral_atom, "dihedral_atom",TagInt,      Tag,     DihedralPA3, Dihedral, 3, kColTopology)      \
    X(dihedral_type, "dihedral_type",int,         Int,     DihedralPA,  Dihedral, 1, kColTopology)      \
    X(nspecial,      "nspecial",     int,         Int,     Fixed3,      Special,  0, kColCount)         \
    X(special,       "special",      TagInt,      Tag,     MaxSpecial,  Special,  1, kColTopology)

/// @brief The schema table, generated from the list above.
inline constexpr FieldSpec kAtomFieldSpecs[] = {
#define ATV_FIELD_SPEC_ROW(member, key, Elem, DTag, Stride, Group, epi, flags)                          \
    FieldSpec{key, ColDType::DTag, ColStride::Stride, ColGroup::Group, epi, flags},
    ATV_ATOM_FIELDS(ATV_FIELD_SPEC_ROW)
#undef ATV_FIELD_SPEC_ROW
};

/// @brief Number of per-atom schema columns.
inline constexpr int kAtomFieldCount =
    static_cast<int>(sizeof(kAtomFieldSpecs) / sizeof(FieldSpec));

/// @brief Compile-time typed column key: the unit of launch-boundary
///        resolution (`frame.span(kAtom_x)`), aligned with the molrs
///        two-level naming (`block` / `key`).
///
/// Keys are constants, not bare strings, so a dtype mismatch is a compile
/// error rather than a runtime surprise.
/// @tparam T Element type the column must carry.
template <typename T>
struct ColKey {
    const char* block; ///< molrs block name (per-atom columns: "atoms")
    const char* key;   ///< molrs column name
};

/// @brief Typed keys for every schema column (kAtom_x, kAtom_vy, …).
#define ATV_FIELD_KEY(member, key, Elem, DTag, Stride, Group, epi, flags)                               \
    inline constexpr ColKey<Elem> kAtom_##member{"atoms", key};
ATV_ATOM_FIELDS(ATV_FIELD_KEY)
#undef ATV_FIELD_KEY

} // namespace ATV
