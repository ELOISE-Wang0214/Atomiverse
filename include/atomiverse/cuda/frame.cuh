/// @file cuda/frame.cuh
/// @brief cuda::Frame -- GPU-side system state (SoA device pointers).
///
/// ALL array data lives on the GPU. Host-cached metadata (n_atoms, n_bonds,
/// etc.) avoids D2H transfers during the simulation loop.
/// CUDA-free header (contains only pointer declarations).
///
/// @see cpu::Frame for the host-side counterpart
/// @see cpu::to_cuda / cuda::to_cpu for transfer
#pragma once

#include "atomiverse/core/int.h"
#include "atomiverse/core/real.h"

namespace ATV {
namespace cuda {

/// @brief Neighbor pair list -- all arrays on device.
///
/// Allocated by NeighborList::build(), populated by NeighborList::launch(),
/// freed by NeighborList::destroy(). Starts as nullptr in a fresh Frame.
struct PairList {
    Int* i = nullptr;         ///< First atom index [device].
    Int* j = nullptr;         ///< Second atom index [device].
    Real* dx = nullptr;       ///< x-displacement (j - i) [Angstrom, device].
    Real* dy = nullptr;       ///< y-displacement (j - i) [Angstrom, device].
    Real* dz = nullptr;       ///< z-displacement (j - i) [Angstrom, device].
    Real* dist = nullptr;     ///< Pair distance [Angstrom, device].
    Int* n_pairs = nullptr;   ///< Current pair count [device, 1 element].
    Int* capacity = nullptr;  ///< Pre-allocated slot count [device, 1 element].
};

/// @brief Bond topology list -- all arrays on device.
struct BondList {
    Int* i = nullptr;        ///< First atom index [device].
    Int* j = nullptr;        ///< Second atom index [device].
    Int* type = nullptr;     ///< Bond type index [device].
    Int* n_bonds = nullptr;  ///< Bond count [device, 1 element].
};

/// @brief Angle topology list -- all arrays on device.
struct AngleList {
    Int* i = nullptr;         ///< First atom index [device].
    Int* j = nullptr;         ///< Center atom index [device].
    Int* k = nullptr;         ///< Third atom index [device].
    Int* type = nullptr;      ///< Angle type index [device].
    Int* n_angles = nullptr;  ///< Angle count [device, 1 element].
};

/// @brief Dihedral topology list -- all arrays on device.
struct DihedralList {
    Int* i = nullptr;            ///< First atom index [device].
    Int* j = nullptr;            ///< Second atom index [device].
    Int* k = nullptr;            ///< Third atom index [device].
    Int* l = nullptr;            ///< Fourth atom index [device].
    Int* type = nullptr;         ///< Dihedral type index [device].
    Int* n_dihedrals = nullptr;  ///< Dihedral count [device, 1 element].
};

/// @brief GPU-side system state -- device pointers + host-cached metadata.
///
/// Device pointer fields are GPU-resident. Metadata fields (n_atoms_, n_bonds_,
/// etc.) are host-cached counts set at materialize time to avoid D2H transfers.
///
/// Created by cpu::to_cuda(), freed by cuda::destroy().
struct Frame {
    // ── Per-atom arrays (device pointers) ──

    Real* x = nullptr;      ///< x positions [Angstrom, device].
    Real* y = nullptr;      ///< y positions [Angstrom, device].
    Real* z = nullptr;      ///< z positions [Angstrom, device].
    Real* vx = nullptr;     ///< x velocities [Angstrom/fs, device].
    Real* vy = nullptr;     ///< y velocities [Angstrom/fs, device].
    Real* vz = nullptr;     ///< z velocities [Angstrom/fs, device].
    Real* fx = nullptr;     ///< x forces [eV/Angstrom, device].
    Real* fy = nullptr;     ///< y forces [eV/Angstrom, device].
    Real* fz = nullptr;     ///< z forces [eV/Angstrom, device].
    Real* mass = nullptr;   ///< Masses [amu, device].
    Int* type_id = nullptr; ///< Atom types [device].

    // ── System scalars (device pointers) ──

    Real* energy = nullptr;                       ///< Potential energy [eV, device, 1 element].
    Real* kinetic = nullptr;                      ///< Kinetic energy [eV, device, 1 element].
    Real* stress = nullptr;                       ///< Stress tensor [eV/Angstrom^3, device, 6 elements, Voigt].
    Real* box = nullptr;                          ///< Simulation cell 3x3 row-major [Angstrom, device, 9 elements].

    // ── Sub-structures ──

    PairList pairs;                               ///< Neighbor pair list (allocated by NeighborList::build).
    BondList bonds;                               ///< Bond topology (transferred by cpu::to_cuda).
    AngleList angles;                             ///< Angle topology (transferred by cpu::to_cuda).
    DihedralList dihedrals;                       ///< Dihedral topology (transferred by cpu::to_cuda).

    /// @brief SCF electronic-structure workspace -- all arrays on device.
    ///
    /// Allocated by SCF components in build(), freed in destroy().
    /// Column-major dense matrices (cuBLAS convention), dimension n x n.
    struct SCFWorkspace {
        // ── Shared (toy-scf + RHF) ──
        Real* H0    = nullptr;  ///< Base Hamiltonian (toy) [device, n*n].
        Real* H     = nullptr;  ///< Effective Hamiltonian / Fock matrix [device, n*n].
        Real* P     = nullptr;  ///< Density matrix [device, n*n].
        Real* P_old = nullptr;  ///< Previous density matrix [device, n*n].
        Real* C     = nullptr;  ///< Eigenvectors / MO coefficients [device, n*n, col-major].
        Real* eig   = nullptr;  ///< Eigenvalues / orbital energies [device, n].
        Real* work  = nullptr;  ///< Scratch buffer [device, n*n].
        Int   n     = 0;        ///< Matrix dimension (host-cached).
        Int   n_occ = 0;        ///< Occupation count (host-cached).
        bool  converged = false; ///< Convergence flag (host-cached).
        int   iteration = 0;     ///< Iteration that converged (host-cached).
        Real  residual  = 0;     ///< Last convergence residual (host-cached).

        // ── RHF-specific ──
        Real* Hcore   = nullptr; ///< Core Hamiltonian T+V [device, n*n].
        Real* X       = nullptr; ///< Orthogonalization matrix S^{-1/2} [device, n*n].
        Real* F_prime = nullptr; ///< Transformed Fock X^T*F*X [device, n*n].
        Real* ERI     = nullptr; ///< Two-electron integrals [device, n^4].
        Real  V_nuc   = 0;      ///< Nuclear repulsion energy (host-cached).

        // ── Basis data for GPU-native integral kernels ──
        Real* d_shell_alpha  = nullptr; ///< Exponents [device, n_shells * n_prim].
        Real* d_shell_coeff  = nullptr; ///< Normalized coefficients [device, n_shells * n_prim].
        Real* d_shell_center = nullptr; ///< Shell centers [device, n_shells * 3] (Bohr).
        Int*  d_atom_Z       = nullptr; ///< Nuclear charges [device, natm].
        Real* d_atom_pos     = nullptr; ///< Atom positions [device, natm * 3] (Bohr).
        Int   n_prim = 0;              ///< Primitives per shell (host-cached).
        Int   n_shells = 0;            ///< Number of shells (host-cached).
        Int   natm = 0;                ///< Number of atoms (host-cached).
    };
    SCFWorkspace scf;                             ///< SCF workspace (allocated by SCF components).

    Int* d_n_atoms = nullptr;                     ///< Atom count [device, 1 element].

    // ── Host-cached metadata (set at materialize, avoids D2H in loop) ──

    /// @brief Number of atoms cached on the host side.
    /// @return Atom count.
    Int n_atoms() const { return n_atoms_; }
    /// @brief Number of bonds cached on the host side.
    /// @return Bond count.
    Int n_bonds() const { return n_bonds_; }
    /// @brief Number of angles cached on the host side.
    /// @return Angle count.
    Int n_angles() const { return n_angles_; }
    /// @brief Number of dihedrals cached on the host side.
    /// @return Dihedral count.
    Int n_dihedrals() const { return n_dihedrals_; }

    Int n_atoms_ = 0;      ///< Host-cached atom count.
    Int n_bonds_ = 0;      ///< Host-cached bond count.
    Int n_angles_ = 0;     ///< Host-cached angle count.
    Int n_dihedrals_ = 0;  ///< Host-cached dihedral count.
};

} // namespace cuda
} // namespace ATV
