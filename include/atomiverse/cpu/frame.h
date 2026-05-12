/// @file cpu/frame.h
/// @brief cpu::Frame -- CPU device system state backed by xtensor.
///
/// This is the primary data container users interact with. CPU is treated
/// as a device, symmetric with CUDA. Transfer between devices uses
/// cpu::to_cuda() and cuda::to_cpu().
///
/// @see cuda::Frame for the GPU device counterpart
/// @see cpu::to_cuda() for CPU -> CUDA transfer
/// @see cuda::to_cpu() for CUDA -> CPU transfer
#pragma once

#include <xtensor/xtensor.hpp>

#include "atomiverse/core/int.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/core/real.h"

namespace ATV {
namespace cpu {

/// @brief CPU device system state backed by xtensor arrays.
///
/// All arrays use SoA layout (one xtensor per dimension).
///
/// Usage:
/// @code{.cpp}
/// cpu::Frame frame(1000);
/// frame.x(0) = 0.0;
/// frame.mass(0) = 12.011;
/// frame.box = Mat3::diagonal(10.0, 10.0, 10.0);
/// @endcode
class Frame {
public:
    /// @brief Construct a Frame for the given number of atoms.
    ///
    /// All arrays are allocated and zero-initialized.
    ///
    /// @param n_atoms Number of atoms in the system.
    explicit Frame(Int n_atoms)
        : x(xt::zeros<Real>({n_atoms})),
          y(xt::zeros<Real>({n_atoms})),
          z(xt::zeros<Real>({n_atoms})),
          vx(xt::zeros<Real>({n_atoms})),
          vy(xt::zeros<Real>({n_atoms})),
          vz(xt::zeros<Real>({n_atoms})),
          fx(xt::zeros<Real>({n_atoms})),
          fy(xt::zeros<Real>({n_atoms})),
          fz(xt::zeros<Real>({n_atoms})),
          mass(xt::zeros<Real>({n_atoms})),
          type_id(xt::zeros<Int>({n_atoms})),
          n_atoms_(n_atoms) {}

    /// @brief Number of atoms in the system.
    /// @return Atom count fixed at construction.
    Int n_atoms() const { return n_atoms_; }

    // ── Per-atom arrays ──

    xt::xtensor<Real, 1> x;             ///< x positions [Angstrom].
    xt::xtensor<Real, 1> y;             ///< y positions [Angstrom].
    xt::xtensor<Real, 1> z;             ///< z positions [Angstrom].
    xt::xtensor<Real, 1> vx;            ///< x velocities [Angstrom/fs].
    xt::xtensor<Real, 1> vy;            ///< y velocities [Angstrom/fs].
    xt::xtensor<Real, 1> vz;            ///< z velocities [Angstrom/fs].
    xt::xtensor<Real, 1> fx;            ///< x forces [eV/Angstrom].
    xt::xtensor<Real, 1> fy;            ///< y forces [eV/Angstrom].
    xt::xtensor<Real, 1> fz;            ///< z forces [eV/Angstrom].
    xt::xtensor<Real, 1> mass;          ///< Masses [amu].
    xt::xtensor<Int, 1>  type_id;       ///< Atom type indices.

    // ── Topology ──

    xt::xtensor<Int, 1> bond_i;            ///< First atom index for each bond.
    xt::xtensor<Int, 1> bond_j;            ///< Second atom index for each bond.
    xt::xtensor<Int, 1> bond_type;         ///< Bond type index.

    xt::xtensor<Int, 1> angle_i;           ///< First atom index for each angle.
    xt::xtensor<Int, 1> angle_j;           ///< Center atom index for each angle.
    xt::xtensor<Int, 1> angle_k;           ///< Third atom index for each angle.
    xt::xtensor<Int, 1> angle_type;        ///< Angle type index.

    xt::xtensor<Int, 1> dihedral_i;        ///< First atom index for each dihedral.
    xt::xtensor<Int, 1> dihedral_j;        ///< Second atom index for each dihedral.
    xt::xtensor<Int, 1> dihedral_k;        ///< Third atom index for each dihedral.
    xt::xtensor<Int, 1> dihedral_l;        ///< Fourth atom index for each dihedral.
    xt::xtensor<Int, 1> dihedral_type;     ///< Dihedral type index.

    // ── System scalars ──

    Real energy = 0;     ///< Potential energy [eV].
    Real kinetic = 0;    ///< Kinetic energy [eV].
    Real stress[6] = {}; ///< Stress tensor [eV/Angstrom^3], Voigt: xx yy zz xy xz yz.

    Mat3 box{};  ///< Simulation cell [Angstrom], 3x3 row-major.

    // ── Component communication ──

    const void* pair_list = nullptr; ///< Non-owning pointer to active PairList (set by NeighborList).

private:
    Int n_atoms_; ///< Cached atom count set at construction.
};

} // namespace cpu
} // namespace ATV
