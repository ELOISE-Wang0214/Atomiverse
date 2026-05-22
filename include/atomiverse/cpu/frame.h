/// @file cpu/frame.h
/// @brief cpu::Frame -- CPU backend execution frame backed by xtensor.
///
/// `cpu::Frame` is a pure execution product, not a user-facing data type.
/// The canonical exchange layer is `molrs.Frame`: a driver's user-facing
/// `run(molrs::FrameRef&, ...)` entry materializes a `molrs.Frame` into a
/// transient `cpu::Frame` (via `cpu::materialize`), executes against it, and
/// writes results back (via `cpu::write_back`). The `cpu::Frame` lives only
/// for the duration of that one `run()` call. CPU is treated as a device,
/// symmetric with CUDA; transfer between device frames uses cpu::to_cuda()
/// and cuda::to_cpu().
///
/// Data is organized into sub-structs (atoms, bonds, angles, dihedrals)
/// matching the molcrafts naming convention (molpy Block model).
///
/// @see molrs.hpp for the canonical molrs.Frame exchange layer
/// @see cpu::materialize / cpu::write_back for the molrs.Frame <-> cpu::Frame bridge
/// @see cuda::Frame for the GPU device counterpart
/// @see cpu::to_cuda() for CPU -> CUDA transfer
/// @see cuda::to_cpu() for CUDA -> CPU transfer
#pragma once

#include <cstdint>

#include <xtensor/containers/xtensor.hpp>
#include <xtensor/generators/xbuilder.hpp>
#include <xtensor/views/xview.hpp>

#include "atomiverse/core/mat3.h"
#include "atomiverse/core/system_scalars.h"
#include "atomiverse/cpu/box.h"

namespace ATV {
namespace cpu {

/// @brief CPU backend execution frame backed by xtensor arrays.
///
/// A transient execution product owned by a driver's `run()` call, not a
/// user-facing data type. Materialized from a `molrs.Frame` on entry and
/// written back on exit. All arrays use SoA layout (one xtensor per
/// dimension). Fields are organized into sub-structs matching molpy's Block
/// model: atoms, bonds, angles, dihedrals.
///
/// Internal usage (driver / test code):
/// @code{.cpp}
/// cpu::Frame frame(1000);
/// frame.atoms.x(0) = 0.0;
/// frame.atoms.mass(0) = 12.011;
/// frame.box = Mat3::cubic(10.0);
/// @endcode
class Frame {
  public:
    /// @brief Construct a Frame for the given number of atoms.
    ///
    /// All per-atom arrays are allocated and zero-initialized.
    /// Topology blocks (bonds, angles, dihedrals) default to empty.
    ///
    /// @param n_atoms Number of atoms in the system.
    explicit Frame(int n_atoms) : atoms(n_atoms) {}

    // ── Per-atom block (molcrafts: "atoms" namespace) ──

    /// @brief Per-atom data block matching molpy atoms namespace.
    struct Atoms {
        /// @brief Construct Atoms block for the given number of atoms.
        /// @param n Number of atoms.
        explicit Atoms(int n)
            : x(xt::zeros<double>({n})), y(xt::zeros<double>({n})),
              z(xt::zeros<double>({n})), vx(xt::zeros<double>({n})),
              vy(xt::zeros<double>({n})), vz(xt::zeros<double>({n})),
              fx(xt::zeros<double>({n})), fy(xt::zeros<double>({n})),
              fz(xt::zeros<double>({n})), mass(xt::zeros<double>({n})),
              charge(xt::zeros<double>({n})), type(xt::zeros<int>({n})),
              tag(xt::arange<std::int64_t>(0, n)),
              image(xt::zeros<std::int8_t>({n, 3})), n_(n) {}

        xt::xtensor<double, 1> x;    ///< x positions [Angstrom].
        xt::xtensor<double, 1> y;    ///< y positions [Angstrom].
        xt::xtensor<double, 1> z;    ///< z positions [Angstrom].
        xt::xtensor<double, 1> vx;   ///< x velocities [Angstrom/fs].
        xt::xtensor<double, 1> vy;   ///< y velocities [Angstrom/fs].
        xt::xtensor<double, 1> vz;   ///< z velocities [Angstrom/fs].
        xt::xtensor<double, 1> fx;   ///< x forces [eV/Angstrom].
        xt::xtensor<double, 1> fy;   ///< y forces [eV/Angstrom].
        xt::xtensor<double, 1> fz;   ///< z forces [eV/Angstrom].
        xt::xtensor<double, 1> mass; ///< Masses [amu].

        /// @brief Per-atom partial charge [e] (elementary charge).
        ///
        /// Zero-initialized. Mirrors `Capability::Charge`: components that
        /// declare that capability write here. Treated as a fundamental
        /// physical property (analogue of LAMMPS atom->q).
        xt::xtensor<double, 1> charge;

        xt::xtensor<int, 1> type; ///< Atom type indices.

        /// @brief Global atom identity, persistent across MPI migration and restart.
        ///
        /// Default-initialized to identity: `tag(i) == i` after `Frame(n)`. Required
        /// for bonded topology to survive `AtomMigration`. Treated as a fundamental
        /// physical property (analogue of LAMMPS atom->tag).
        xt::xtensor<std::int64_t, 1> tag;

        /// @brief Per-atom PBC crossing counts, shape [n, 3].
        ///
        /// `image(i, d)` increments by ±1 when atom i crosses the box face along
        /// dimension d. Zero-initialized. Used to reconstruct unwrapped trajectories
        /// (diffusion, MSD). Treated as a fundamental physical property
        /// (analogue of LAMMPS atom->image).
        xt::xtensor<std::int8_t, 2> image;

        /// @brief Number of atoms.
        /// @return Atom count.
        int size() const { return n_; }

      private:
        int n_;
    } atoms;

    // ── Bond topology block (molcrafts: "bonds" namespace) ──

    /// @brief Bond topology block matching molpy bonds namespace.
    struct Bonds {
        xt::xtensor<int, 1> atomi; ///< First atom index for each bond.
        xt::xtensor<int, 1> atomj; ///< Second atom index for each bond.
        xt::xtensor<int, 1> type;  ///< Bond type index.

        /// @brief Number of bonds.
        /// @return Bond count.
        int size() const { return static_cast<int>(atomi.size()); }
    } bonds;

    // ── Angle topology block (molcrafts: "angles" namespace) ──

    /// @brief Angle topology block matching molpy angles namespace.
    struct Angles {
        xt::xtensor<int, 1> atomi; ///< First atom index for each angle.
        xt::xtensor<int, 1> atomj; ///< Center atom index for each angle.
        xt::xtensor<int, 1> atomk; ///< Third atom index for each angle.
        xt::xtensor<int, 1> type;  ///< Angle type index.

        /// @brief Number of angles.
        /// @return Angle count.
        int size() const { return static_cast<int>(atomi.size()); }
    } angles;

    // ── Dihedral topology block (molcrafts: "dihedrals" namespace) ──

    /// @brief Dihedral topology block matching molpy dihedrals namespace.
    struct Dihedrals {
        xt::xtensor<int, 1> atomi; ///< First atom index for each dihedral.
        xt::xtensor<int, 1> atomj; ///< Second atom index for each dihedral.
        xt::xtensor<int, 1> atomk; ///< Third atom index for each dihedral.
        xt::xtensor<int, 1> atoml; ///< Fourth atom index for each dihedral.
        xt::xtensor<int, 1> type;  ///< Dihedral type index.

        /// @brief Number of dihedrals.
        /// @return Dihedral count.
        int size() const { return static_cast<int>(atomi.size()); }
    } dihedrals;

    // ── System scalars ──

    /// @brief Per-frame system scalar results (energy, kinetic, conserved,
    /// stress, kinetic stress). See `core/system_scalars.h` for the
    /// field-by-field documentation; the shape is shared with the CUDA and
    /// Torch frames.
    SystemScalars scalars{};

    cpu::Box
        box{}; ///< Simulation cell [Angstrom]; canonical Box at backend root (see cpu/box.h).
};

} // namespace cpu
} // namespace ATV
