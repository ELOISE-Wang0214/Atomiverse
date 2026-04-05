/// @file python/cpu_frame_wrapper.h
/// @brief CpuFrame -- Python-visible wrapper around cpu::Frame.
#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include "atomiverse/cpu/frame.h"
#include "atomiverse/core/real.h"
#include "atomiverse/core/int.h"

#include <cstring>
#include <stdexcept>

namespace nb = nanobind;
using namespace ATV;

using NpReal = nb::ndarray<Real, nb::ndim<1>, nb::device::cpu, nb::c_contig>;
using NpInt  = nb::ndarray<Int,  nb::ndim<1>, nb::device::cpu, nb::c_contig>;

/// @brief Python-visible wrapper around cpu::Frame.
///
/// Owns a cpu::Frame. Provides populate/extract methods for numpy arrays.
struct CpuFrame {
    cpu::Frame frame;

    explicit CpuFrame(Int n_atoms) : frame(n_atoms) {}

    Int n_atoms() const { return frame.atoms.size(); }

    // ── Populate from numpy arrays ──

    void set_positions(NpReal x, NpReal y, NpReal z) {
        Int n = frame.atoms.size();
        if (static_cast<Int>(x.shape(0)) != n ||
            static_cast<Int>(y.shape(0)) != n ||
            static_cast<Int>(z.shape(0)) != n)
            throw std::invalid_argument("position array size mismatch");
        std::memcpy(frame.atoms.x.data(), x.data(), n * sizeof(Real));
        std::memcpy(frame.atoms.y.data(), y.data(), n * sizeof(Real));
        std::memcpy(frame.atoms.z.data(), z.data(), n * sizeof(Real));
    }

    void set_velocities(NpReal vx, NpReal vy, NpReal vz) {
        Int n = frame.atoms.size();
        std::memcpy(frame.atoms.vx.data(), vx.data(), n * sizeof(Real));
        std::memcpy(frame.atoms.vy.data(), vy.data(), n * sizeof(Real));
        std::memcpy(frame.atoms.vz.data(), vz.data(), n * sizeof(Real));
    }

    void set_mass(NpReal mass) {
        Int n = frame.atoms.size();
        if (static_cast<Int>(mass.shape(0)) != n)
            throw std::invalid_argument("mass array size mismatch");
        std::memcpy(frame.atoms.mass.data(), mass.data(), n * sizeof(Real));
    }

    void set_type(NpInt type) {
        Int n = frame.atoms.size();
        if (static_cast<Int>(type.shape(0)) != n)
            throw std::invalid_argument("type array size mismatch");
        std::memcpy(frame.atoms.type.data(), type.data(), n * sizeof(Int));
    }

    void set_box(nb::ndarray<Real, nb::shape<3, 3>, nb::device::cpu, nb::c_contig> box) {
        const Real* bp = box.data();
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                frame.box(i, j) = bp[i * 3 + j];
    }

    // ── Topology ──

    void set_bonds(NpInt atomi, NpInt atomj, NpInt type) {
        Int nb_ = static_cast<Int>(atomi.shape(0));
        frame.bonds.atomi = xt::zeros<Int>({nb_});
        frame.bonds.atomj = xt::zeros<Int>({nb_});
        frame.bonds.type  = xt::zeros<Int>({nb_});
        std::memcpy(frame.bonds.atomi.data(), atomi.data(), nb_ * sizeof(Int));
        std::memcpy(frame.bonds.atomj.data(), atomj.data(), nb_ * sizeof(Int));
        std::memcpy(frame.bonds.type.data(),  type.data(),  nb_ * sizeof(Int));
    }

    void set_angles(NpInt atomi, NpInt atomj, NpInt atomk, NpInt type) {
        Int na = static_cast<Int>(atomi.shape(0));
        frame.angles.atomi = xt::zeros<Int>({na});
        frame.angles.atomj = xt::zeros<Int>({na});
        frame.angles.atomk = xt::zeros<Int>({na});
        frame.angles.type  = xt::zeros<Int>({na});
        std::memcpy(frame.angles.atomi.data(), atomi.data(), na * sizeof(Int));
        std::memcpy(frame.angles.atomj.data(), atomj.data(), na * sizeof(Int));
        std::memcpy(frame.angles.atomk.data(), atomk.data(), na * sizeof(Int));
        std::memcpy(frame.angles.type.data(),  type.data(),  na * sizeof(Int));
    }

    void set_dihedrals(NpInt atomi, NpInt atomj, NpInt atomk, NpInt atoml, NpInt type) {
        Int nd = static_cast<Int>(atomi.shape(0));
        frame.dihedrals.atomi = xt::zeros<Int>({nd});
        frame.dihedrals.atomj = xt::zeros<Int>({nd});
        frame.dihedrals.atomk = xt::zeros<Int>({nd});
        frame.dihedrals.atoml = xt::zeros<Int>({nd});
        frame.dihedrals.type  = xt::zeros<Int>({nd});
        std::memcpy(frame.dihedrals.atomi.data(), atomi.data(), nd * sizeof(Int));
        std::memcpy(frame.dihedrals.atomj.data(), atomj.data(), nd * sizeof(Int));
        std::memcpy(frame.dihedrals.atomk.data(), atomk.data(), nd * sizeof(Int));
        std::memcpy(frame.dihedrals.atoml.data(), atoml.data(), nd * sizeof(Int));
        std::memcpy(frame.dihedrals.type.data(),  type.data(),  nd * sizeof(Int));
    }

    // ── Extract helpers ──

    template <typename T>
    static nb::ndarray<nb::numpy, T, nb::ndim<1>> copy_to_numpy(const T* src, Int n) {
        T* buf = new T[n];
        std::memcpy(buf, src, n * sizeof(T));
        nb::capsule owner(buf, [](void* p) noexcept { delete[] static_cast<T*>(p); });
        return nb::ndarray<nb::numpy, T, nb::ndim<1>>(buf, {static_cast<size_t>(n)}, owner);
    }

    nb::ndarray<nb::numpy, Real, nb::ndim<1>> get_x()  const { return copy_to_numpy(frame.atoms.x.data(),  frame.atoms.size()); }
    nb::ndarray<nb::numpy, Real, nb::ndim<1>> get_y()  const { return copy_to_numpy(frame.atoms.y.data(),  frame.atoms.size()); }
    nb::ndarray<nb::numpy, Real, nb::ndim<1>> get_z()  const { return copy_to_numpy(frame.atoms.z.data(),  frame.atoms.size()); }
    nb::ndarray<nb::numpy, Real, nb::ndim<1>> get_vx() const { return copy_to_numpy(frame.atoms.vx.data(), frame.atoms.size()); }
    nb::ndarray<nb::numpy, Real, nb::ndim<1>> get_vy() const { return copy_to_numpy(frame.atoms.vy.data(), frame.atoms.size()); }
    nb::ndarray<nb::numpy, Real, nb::ndim<1>> get_vz() const { return copy_to_numpy(frame.atoms.vz.data(), frame.atoms.size()); }
    nb::ndarray<nb::numpy, Real, nb::ndim<1>> get_fx() const { return copy_to_numpy(frame.atoms.fx.data(), frame.atoms.size()); }
    nb::ndarray<nb::numpy, Real, nb::ndim<1>> get_fy() const { return copy_to_numpy(frame.atoms.fy.data(), frame.atoms.size()); }
    nb::ndarray<nb::numpy, Real, nb::ndim<1>> get_fz() const { return copy_to_numpy(frame.atoms.fz.data(), frame.atoms.size()); }
    nb::ndarray<nb::numpy, Real, nb::ndim<1>> get_mass() const { return copy_to_numpy(frame.atoms.mass.data(), frame.atoms.size()); }
    nb::ndarray<nb::numpy, Int,  nb::ndim<1>> get_type() const { return copy_to_numpy(frame.atoms.type.data(), frame.atoms.size()); }

    Real get_energy()  const { return frame.energy; }
    Real get_kinetic() const { return frame.kinetic; }

    nb::ndarray<nb::numpy, Real, nb::shape<3, 3>> get_box() const {
        Real* buf = new Real[9];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                buf[i * 3 + j] = frame.box(i, j);
        nb::capsule owner(buf, [](void* p) noexcept { delete[] static_cast<Real*>(p); });
        return nb::ndarray<nb::numpy, Real, nb::shape<3, 3>>(buf, {3, 3}, owner);
    }
};
