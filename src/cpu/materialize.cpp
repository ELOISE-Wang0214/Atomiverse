/// @file src/cpu/materialize.cpp
/// @brief Implementation of cpu::materialize / cpu::write_back.

#include "atomiverse/cpu/materialize.h"

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "atomiverse/core/mat3.h"

namespace ATV {
namespace cpu {

namespace {

// ── Required-field readers (fail loud, naming the field) ────────────────

/// @brief Read a required f64 column or throw with the field name.
::rust::Vec<double> require_f64(const molrs::FrameRef& src, const char* block,
                                const char* col) {
    auto v = molrs::frame_column_f64(src, block, col);
    if (v.empty()) {
        throw std::runtime_error(std::string("cpu::materialize: required "
                                             "field '")
                                 + block + "." + col
                                 + "' is missing or empty in the molrs.Frame");
    }
    return v;
}

/// @brief Read a required i32 column or throw with the field name.
::rust::Vec<std::int32_t> require_i32(const molrs::FrameRef& src,
                                      const char* block, const char* col) {
    auto v = molrs::frame_column_i32(src, block, col);
    if (v.empty()) {
        throw std::runtime_error(std::string("cpu::materialize: required "
                                             "field '")
                                 + block + "." + col
                                 + "' is missing or empty in the molrs.Frame");
    }
    return v;
}

// ── Optional-field copy helpers (silent zero on absent column) ──────────

/// @brief Copy an optional f64 column into `dst[0..n)`; leave zeros if absent.
void copy_optional_f64(const molrs::FrameRef& src, const char* block,
                       const char* col, xt::xtensor<double, 1>& dst, int n) {
    auto v = molrs::frame_column_f64(src, block, col);
    if (static_cast<int>(v.size()) != n) {
        return; // absent or shape mismatch -> keep the zero-initialized array
    }
    for (int i = 0; i < n; ++i) {
        dst[i] = v[i];
    }
}

/// @brief Materialize an optional i32 topology column into `dst`.
xt::xtensor<int, 1> read_topology_column(const molrs::FrameRef& src,
                                         const char* block, const char* col,
                                         int n) {
    auto v = molrs::frame_column_i32(src, block, col);
    xt::xtensor<int, 1> out = xt::zeros<int>({n});
    if (static_cast<int>(v.size()) == n) {
        for (int i = 0; i < n; ++i) {
            out[i] = v[i];
        }
    }
    return out;
}

} // namespace

Frame materialize(const molrs::FrameRef& src) {
    // ── Required: atom positions + type ────────────────────────────────
    auto x = require_f64(src, "atoms", "x");
    auto y = require_f64(src, "atoms", "y");
    auto z = require_f64(src, "atoms", "z");
    auto type = require_i32(src, "atoms", "type");

    const int n = static_cast<int>(x.size());
    if (static_cast<int>(y.size()) != n || static_cast<int>(z.size()) != n) {
        throw std::runtime_error("cpu::materialize: atoms.x / atoms.y / "
                                 "atoms.z have inconsistent row counts");
    }
    if (static_cast<int>(type.size()) != n) {
        throw std::runtime_error("cpu::materialize: atoms.type row count "
                                 "does not match atoms.x");
    }

    Frame frame(n);
    for (int i = 0; i < n; ++i) {
        frame.atoms.x[i] = x[i];
        frame.atoms.y[i] = y[i];
        frame.atoms.z[i] = z[i];
        frame.atoms.type[i] = type[i];
    }

    // ── Optional per-atom fields (silent zero when absent) ─────────────
    copy_optional_f64(src, "atoms", "vx", frame.atoms.vx, n);
    copy_optional_f64(src, "atoms", "vy", frame.atoms.vy, n);
    copy_optional_f64(src, "atoms", "vz", frame.atoms.vz, n);
    copy_optional_f64(src, "atoms", "fx", frame.atoms.fx, n);
    copy_optional_f64(src, "atoms", "fy", frame.atoms.fy, n);
    copy_optional_f64(src, "atoms", "fz", frame.atoms.fz, n);
    copy_optional_f64(src, "atoms", "mass", frame.atoms.mass, n);
    copy_optional_f64(src, "atoms", "charge", frame.atoms.charge, n);

    // Optional global tag (i32 in molrs -> i64 in cpu::Frame). When absent
    // the Frame ctor's identity default (tag(i) == i) is kept.
    {
        auto tag = molrs::frame_column_i32(src, "atoms", "tag");
        if (static_cast<int>(tag.size()) == n) {
            for (int i = 0; i < n; ++i) {
                frame.atoms.tag[i] = static_cast<std::int64_t>(tag[i]);
            }
        }
    }

    // Optional per-axis PBC image counts (i32 in molrs -> i8 [n,3] in
    // cpu::Frame). Absent columns leave the zero-initialized image array.
    {
        const char* axis_cols[3] = {"ix", "iy", "iz"};
        for (int d = 0; d < 3; ++d) {
            auto img = molrs::frame_column_i32(src, "atoms", axis_cols[d]);
            if (static_cast<int>(img.size()) == n) {
                for (int i = 0; i < n; ++i) {
                    frame.atoms.image(i, d) = static_cast<std::int8_t>(img[i]);
                }
            }
        }
    }

    // ── Optional topology blocks ───────────────────────────────────────
    {
        const int nb = static_cast<int>(molrs::frame_block_nrows(src, "bonds"));
        if (nb > 0) {
            frame.bonds.atomi = read_topology_column(src, "bonds", "atomi", nb);
            frame.bonds.atomj = read_topology_column(src, "bonds", "atomj", nb);
            frame.bonds.type = read_topology_column(src, "bonds", "type", nb);
        }
    }
    {
        const int na =
            static_cast<int>(molrs::frame_block_nrows(src, "angles"));
        if (na > 0) {
            frame.angles.atomi =
                read_topology_column(src, "angles", "atomi", na);
            frame.angles.atomj =
                read_topology_column(src, "angles", "atomj", na);
            frame.angles.atomk =
                read_topology_column(src, "angles", "atomk", na);
            frame.angles.type = read_topology_column(src, "angles", "type", na);
        }
    }
    {
        const int nd =
            static_cast<int>(molrs::frame_block_nrows(src, "dihedrals"));
        if (nd > 0) {
            frame.dihedrals.atomi =
                read_topology_column(src, "dihedrals", "atomi", nd);
            frame.dihedrals.atomj =
                read_topology_column(src, "dihedrals", "atomj", nd);
            frame.dihedrals.atomk =
                read_topology_column(src, "dihedrals", "atomk", nd);
            frame.dihedrals.atoml =
                read_topology_column(src, "dihedrals", "atoml", nd);
            frame.dihedrals.type =
                read_topology_column(src, "dihedrals", "type", nd);
        }
    }

    // ── Optional simulation box (9 row-major H entries) ────────────────
    {
        auto h = molrs::frame_simbox(src);
        if (h.size() == 9u) {
            Mat3 cell;
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    cell(i, j) = static_cast<Real>(h[i * 3 + j]);
                }
            }
            frame.box.set_cell(cell);
        }
    }

    return frame;
}

void write_back(const Frame& src, molrs::FrameRef& dst) {
    const int n = src.atoms.size();
    auto slice = [n](const xt::xtensor<double, 1>& a) {
        return ::rust::Slice<const double>(a.data(),
                                           static_cast<std::size_t>(n));
    };

    // ── Per-atom results ───────────────────────────────────────────────
    molrs::frame_set_column_f64(dst, "atoms", "x", slice(src.atoms.x));
    molrs::frame_set_column_f64(dst, "atoms", "y", slice(src.atoms.y));
    molrs::frame_set_column_f64(dst, "atoms", "z", slice(src.atoms.z));
    molrs::frame_set_column_f64(dst, "atoms", "vx", slice(src.atoms.vx));
    molrs::frame_set_column_f64(dst, "atoms", "vy", slice(src.atoms.vy));
    molrs::frame_set_column_f64(dst, "atoms", "vz", slice(src.atoms.vz));
    molrs::frame_set_column_f64(dst, "atoms", "fx", slice(src.atoms.fx));
    molrs::frame_set_column_f64(dst, "atoms", "fy", slice(src.atoms.fy));
    molrs::frame_set_column_f64(dst, "atoms", "fz", slice(src.atoms.fz));
    molrs::frame_set_column_f64(dst, "atoms", "charge",
                                slice(src.atoms.charge));

    // ── System scalars -> dedicated molrs blocks ───────────────────────
    //
    // A molrs block requires every column to share the same row count, so
    // the three scalar quantities go into a 1-row "scalars" block while the
    // two 6-component Voigt tensors each get their own 6-row block. Every
    // tensor lives under a single "value" column.
    const double energy[1] = {src.scalars.energy};
    const double kinetic[1] = {src.scalars.kinetic};
    const double conserved[1] = {src.scalars.conserved};
    molrs::frame_set_column_f64(dst, "scalars", "energy",
                                ::rust::Slice<const double>(energy, 1));
    molrs::frame_set_column_f64(dst, "scalars", "kinetic",
                                ::rust::Slice<const double>(kinetic, 1));
    molrs::frame_set_column_f64(dst, "scalars", "conserved",
                                ::rust::Slice<const double>(conserved, 1));
    molrs::frame_set_column_f64(
        dst, "stress", "value",
        ::rust::Slice<const double>(src.scalars.stress, 6));
    molrs::frame_set_column_f64(
        dst, "stress_kinetic", "value",
        ::rust::Slice<const double>(src.scalars.stress_kinetic, 6));
}

} // namespace cpu
} // namespace ATV
