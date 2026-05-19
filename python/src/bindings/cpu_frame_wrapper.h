/// @file python/cpu_frame_wrapper.h
/// @brief CpuFrame -- Python-visible wrapper around cpu::Frame.
#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include "atomiverse/cpu/frame.h"

#include <cstring>
#include <stdexcept>

namespace nb = nanobind;
using namespace ATV;

using NpReal = nb::ndarray<double, nb::ndim<1>, nb::device::cpu, nb::c_contig>;
using NpInt  = nb::ndarray<int,  nb::ndim<1>, nb::device::cpu, nb::c_contig>;

/// @brief Python-visible wrapper around cpu::Frame.
///
/// Owns a cpu::Frame. Provides populate/extract methods for numpy arrays.
struct CpuFrame {
    cpu::Frame frame;

    explicit CpuFrame(int n_atoms) : frame(n_atoms) {}

    int n_atoms() const { return frame.atoms.size(); }

    // ── Populate from numpy arrays ──

    void set_positions(NpReal x, NpReal y, NpReal z) {
        int n = frame.atoms.size();
        if (static_cast<int>(x.shape(0)) != n ||
            static_cast<int>(y.shape(0)) != n ||
            static_cast<int>(z.shape(0)) != n)
            throw std::invalid_argument("position array size mismatch");
        std::memcpy(frame.atoms.x.data(), x.data(), n * sizeof(double));
        std::memcpy(frame.atoms.y.data(), y.data(), n * sizeof(double));
        std::memcpy(frame.atoms.z.data(), z.data(), n * sizeof(double));
    }

    void set_velocities(NpReal vx, NpReal vy, NpReal vz) {
        int n = frame.atoms.size();
        std::memcpy(frame.atoms.vx.data(), vx.data(), n * sizeof(double));
        std::memcpy(frame.atoms.vy.data(), vy.data(), n * sizeof(double));
        std::memcpy(frame.atoms.vz.data(), vz.data(), n * sizeof(double));
    }

    void set_mass(NpReal mass) {
        int n = frame.atoms.size();
        if (static_cast<int>(mass.shape(0)) != n)
            throw std::invalid_argument("mass array size mismatch");
        std::memcpy(frame.atoms.mass.data(), mass.data(), n * sizeof(double));
    }

    void set_type(NpInt type) {
        int n = frame.atoms.size();
        if (static_cast<int>(type.shape(0)) != n)
            throw std::invalid_argument("type array size mismatch");
        std::memcpy(frame.atoms.type.data(), type.data(), n * sizeof(int));
    }

    void set_box(nb::ndarray<double, nb::shape<3, 3>, nb::device::cpu, nb::c_contig> box) {
        const double* bp = box.data();
        Mat3 cell;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                cell(i, j) = bp[i * 3 + j];
        frame.box = cell;
    }

    // ── Topology ──

    void set_bonds(NpInt atomi, NpInt atomj, NpInt type) {
        int nb_ = static_cast<int>(atomi.shape(0));
        frame.bonds.atomi = xt::zeros<int>({nb_});
        frame.bonds.atomj = xt::zeros<int>({nb_});
        frame.bonds.type  = xt::zeros<int>({nb_});
        std::memcpy(frame.bonds.atomi.data(), atomi.data(), nb_ * sizeof(int));
        std::memcpy(frame.bonds.atomj.data(), atomj.data(), nb_ * sizeof(int));
        std::memcpy(frame.bonds.type.data(),  type.data(),  nb_ * sizeof(int));
    }

    void set_angles(NpInt atomi, NpInt atomj, NpInt atomk, NpInt type) {
        int na = static_cast<int>(atomi.shape(0));
        frame.angles.atomi = xt::zeros<int>({na});
        frame.angles.atomj = xt::zeros<int>({na});
        frame.angles.atomk = xt::zeros<int>({na});
        frame.angles.type  = xt::zeros<int>({na});
        std::memcpy(frame.angles.atomi.data(), atomi.data(), na * sizeof(int));
        std::memcpy(frame.angles.atomj.data(), atomj.data(), na * sizeof(int));
        std::memcpy(frame.angles.atomk.data(), atomk.data(), na * sizeof(int));
        std::memcpy(frame.angles.type.data(),  type.data(),  na * sizeof(int));
    }

    void set_dihedrals(NpInt atomi, NpInt atomj, NpInt atomk, NpInt atoml, NpInt type) {
        int nd = static_cast<int>(atomi.shape(0));
        frame.dihedrals.atomi = xt::zeros<int>({nd});
        frame.dihedrals.atomj = xt::zeros<int>({nd});
        frame.dihedrals.atomk = xt::zeros<int>({nd});
        frame.dihedrals.atoml = xt::zeros<int>({nd});
        frame.dihedrals.type  = xt::zeros<int>({nd});
        std::memcpy(frame.dihedrals.atomi.data(), atomi.data(), nd * sizeof(int));
        std::memcpy(frame.dihedrals.atomj.data(), atomj.data(), nd * sizeof(int));
        std::memcpy(frame.dihedrals.atomk.data(), atomk.data(), nd * sizeof(int));
        std::memcpy(frame.dihedrals.atoml.data(), atoml.data(), nd * sizeof(int));
        std::memcpy(frame.dihedrals.type.data(),  type.data(),  nd * sizeof(int));
    }

    // ── Extract helpers ──

    template <typename T>
    static nb::ndarray<nb::numpy, T, nb::ndim<1>> copy_to_numpy(const T* src, int n) {
        T* buf = new T[n];
        std::memcpy(buf, src, n * sizeof(T));
        nb::capsule owner(buf, [](void* p) noexcept { delete[] static_cast<T*>(p); });
        return nb::ndarray<nb::numpy, T, nb::ndim<1>>(buf, {static_cast<size_t>(n)}, owner);
    }

    nb::ndarray<nb::numpy, double, nb::ndim<1>> get_x()  const { return copy_to_numpy(frame.atoms.x.data(),  frame.atoms.size()); }
    nb::ndarray<nb::numpy, double, nb::ndim<1>> get_y()  const { return copy_to_numpy(frame.atoms.y.data(),  frame.atoms.size()); }
    nb::ndarray<nb::numpy, double, nb::ndim<1>> get_z()  const { return copy_to_numpy(frame.atoms.z.data(),  frame.atoms.size()); }
    nb::ndarray<nb::numpy, double, nb::ndim<1>> get_vx() const { return copy_to_numpy(frame.atoms.vx.data(), frame.atoms.size()); }
    nb::ndarray<nb::numpy, double, nb::ndim<1>> get_vy() const { return copy_to_numpy(frame.atoms.vy.data(), frame.atoms.size()); }
    nb::ndarray<nb::numpy, double, nb::ndim<1>> get_vz() const { return copy_to_numpy(frame.atoms.vz.data(), frame.atoms.size()); }
    nb::ndarray<nb::numpy, double, nb::ndim<1>> get_fx() const { return copy_to_numpy(frame.atoms.fx.data(), frame.atoms.size()); }
    nb::ndarray<nb::numpy, double, nb::ndim<1>> get_fy() const { return copy_to_numpy(frame.atoms.fy.data(), frame.atoms.size()); }
    nb::ndarray<nb::numpy, double, nb::ndim<1>> get_fz() const { return copy_to_numpy(frame.atoms.fz.data(), frame.atoms.size()); }
    nb::ndarray<nb::numpy, double, nb::ndim<1>> get_mass() const { return copy_to_numpy(frame.atoms.mass.data(), frame.atoms.size()); }
    nb::ndarray<nb::numpy, int,  nb::ndim<1>> get_type() const { return copy_to_numpy(frame.atoms.type.data(), frame.atoms.size()); }

    double get_energy()  const { return frame.energy; }
    double get_kinetic() const { return frame.kinetic; }

    nb::ndarray<nb::numpy, double, nb::shape<3, 3>> get_box() const {
        double* buf = new double[9];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                buf[i * 3 + j] = frame.box(i, j);
        nb::capsule owner(buf, [](void* p) noexcept { delete[] static_cast<double*>(p); });
        return nb::ndarray<nb::numpy, double, nb::shape<3, 3>>(buf, {3, 3}, owner);
    }
};
