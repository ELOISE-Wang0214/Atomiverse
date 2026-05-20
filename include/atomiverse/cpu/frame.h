/// @file cpu/frame.h
/// @brief cpu::Frame -- CPU device system state backed by xtensor.
///
/// This is the primary data container users interact with. CPU is treated
/// as a device, symmetric with CUDA. Transfer between devices uses
/// cpu::to_cuda() and cuda::to_cpu().
///
/// Data is organized into sub-structs (atoms, bonds, angles, dihedrals)
/// matching the molcrafts naming convention (molpy Block model).
///
/// @see cuda::Frame for the GPU device counterpart
/// @see cpu::to_cuda() for CPU -> CUDA transfer
/// @see cuda::to_cpu() for CUDA -> CPU transfer
#pragma once

#include <cstdint>

#include <xtensor/containers/xtensor.hpp>
#include <xtensor/generators/xbuilder.hpp>
#include <xtensor/views/xview.hpp>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/box.h"

namespace ATV {
namespace cpu {

/// @brief CPU device system state backed by xtensor arrays.
///
/// All arrays use SoA layout (one xtensor per dimension).
/// Fields are organized into sub-structs matching molpy's Block model:
/// atoms, bonds, angles, dihedrals.
///
/// Usage:
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
    explicit Frame(int n_atoms)
        : atoms(n_atoms) {}

    // ── Per-atom block (molcrafts: "atoms" namespace) ──

    /// @brief Per-atom data block matching molpy atoms namespace.
    struct Atoms {
        /// @brief Construct Atoms block for the given number of atoms.
        /// @param n Number of atoms.
        explicit Atoms(int n)
            : x(xt::zeros<double>({n})),
              y(xt::zeros<double>({n})),
              z(xt::zeros<double>({n})),
              vx(xt::zeros<double>({n})),
              vy(xt::zeros<double>({n})),
              vz(xt::zeros<double>({n})),
              fx(xt::zeros<double>({n})),
              fy(xt::zeros<double>({n})),
              fz(xt::zeros<double>({n})),
              mass(xt::zeros<double>({n})),
              type(xt::zeros<int>({n})),
              tag(xt::arange<std::int64_t>(0, n)),
              image(xt::zeros<std::int8_t>({n, 3})),
              n_(n) {}

        xt::xtensor<double, 1> x;      ///< x positions [Angstrom].
        xt::xtensor<double, 1> y;      ///< y positions [Angstrom].
        xt::xtensor<double, 1> z;      ///< z positions [Angstrom].
        xt::xtensor<double, 1> vx;     ///< x velocities [Angstrom/fs].
        xt::xtensor<double, 1> vy;     ///< y velocities [Angstrom/fs].
        xt::xtensor<double, 1> vz;     ///< z velocities [Angstrom/fs].
        xt::xtensor<double, 1> fx;     ///< x forces [eV/Angstrom].
        xt::xtensor<double, 1> fy;     ///< y forces [eV/Angstrom].
        xt::xtensor<double, 1> fz;     ///< z forces [eV/Angstrom].
        xt::xtensor<double, 1> mass;   ///< Masses [amu].
        xt::xtensor<int, 1>  type;     ///< Atom type indices.

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
        xt::xtensor<int, 1> atomi;   ///< First atom index for each bond.
        xt::xtensor<int, 1> atomj;   ///< Second atom index for each bond.
        xt::xtensor<int, 1> type;    ///< Bond type index.

        /// @brief Number of bonds.
        /// @return Bond count.
        int size() const { return static_cast<int>(atomi.size()); }
    } bonds;

    // ── Angle topology block (molcrafts: "angles" namespace) ──

    /// @brief Angle topology block matching molpy angles namespace.
    struct Angles {
        xt::xtensor<int, 1> atomi;   ///< First atom index for each angle.
        xt::xtensor<int, 1> atomj;   ///< Center atom index for each angle.
        xt::xtensor<int, 1> atomk;   ///< Third atom index for each angle.
        xt::xtensor<int, 1> type;    ///< Angle type index.

        /// @brief Number of angles.
        /// @return Angle count.
        int size() const { return static_cast<int>(atomi.size()); }
    } angles;

    // ── Dihedral topology block (molcrafts: "dihedrals" namespace) ──

    /// @brief Dihedral topology block matching molpy dihedrals namespace.
    struct Dihedrals {
        xt::xtensor<int, 1> atomi;   ///< First atom index for each dihedral.
        xt::xtensor<int, 1> atomj;   ///< Second atom index for each dihedral.
        xt::xtensor<int, 1> atomk;   ///< Third atom index for each dihedral.
        xt::xtensor<int, 1> atoml;   ///< Fourth atom index for each dihedral.
        xt::xtensor<int, 1> type;    ///< Dihedral type index.

        /// @brief Number of dihedrals.
        /// @return Dihedral count.
        int size() const { return static_cast<int>(atomi.size()); }
    } dihedrals;

    // ── System scalars ──

    double energy = 0;     ///< Potential energy [eV].
    double kinetic = 0;    ///< Kinetic energy [eV].
    /// @brief Integrator's conserved quantity [eV].
    ///
    /// For an extended-Lagrangian integrator this is the conserved
    /// *extended* Hamiltonian — e.g. `NoseHooverNPT` writes the full
    /// MTK Hamiltonian (atom KE+PE + barostat KE + Nosé-Hoover chain
    /// KE/PE + p·V). For a plain symplectic integrator it is just
    /// `KE+PE`. The integrator writes it once per step; it stays `0`
    /// when no integrator populates it. A flat `conserved` trace is the
    /// canonical correctness check that an integrator is symplectic;
    /// `Thermo` logs it as a column.
    double conserved = 0;
    /// @brief Pair / bonded virial tensor Σ r·F [eV], Voigt order xx yy zz xy xz yz.
    ///
    /// Written by pair potentials (LJ, Coulomb, ...) and bonded
    /// evaluators at `post_force`. **Does not** include the kinetic
    /// contribution Σ m·v·v -- that lives in `stress_kinetic` so each
    /// piece is consistent with its own time slice (virial is fixed
    /// once forces are computed, kinetic varies as velocities are
    /// updated through the Trotter step). The full pressure-tensor
    /// sum is `stress[i] + stress_kinetic[i]`, and dividing the trace
    /// by `3·V` gives the internal pressure under the LAMMPS sign
    /// convention (positive ⇔ compressive).
    double stress[6] = {};
    /// @brief Kinetic stress tensor Σ m·v·v [eV], Voigt order xx yy zz xy xz yz.
    ///
    /// Written by `VirialStress` at `post_force` using the velocities
    /// at that point. Time-dependent: integrators like
    /// `NoseHooverNPT` that update velocities between `post_force`
    /// and pressure evaluation recompute the kinetic stress fresh
    /// from the post-update velocities rather than reusing this
    /// stored value, so this field is *informational* for any
    /// downstream observer/recorder that wants the post-force
    /// pressure snapshot.
    double stress_kinetic[6] = {};

    cpu::Box box{};  ///< Simulation cell [Angstrom]; canonical Box at backend root (see cpu/box.h).
};

} // namespace cpu
} // namespace ATV
