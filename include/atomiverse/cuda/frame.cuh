/// @file cuda/frame.cuh
/// @brief cuda::Frame -- GPU-side system state (SoA device pointers).
///
/// ALL array data lives on the GPU. Host-cached metadata (n_atoms, n_bonds,
/// etc.) avoids D2H transfers during the simulation loop.
/// CUDA-free header (contains only pointer declarations).
///
/// Data is organized into sub-structs (atoms, bonds, angles, dihedrals)
/// matching the molcrafts naming convention (molpy Block model).
///
/// @see cpu::Frame for the host-side counterpart
/// @see cpu::to_cuda / cuda::to_cpu for transfer
#pragma once

#include "atomiverse/cuda/types.cuh"

namespace ATV {
namespace cuda {

/// @brief Neighbor pair list -- all arrays on device.
///
/// Allocated by NeighborList::build(), populated by NeighborList::launch(),
/// freed by NeighborList::destroy(). Starts as nullptr in a fresh Frame.
struct PairList {
    int* i = nullptr;         ///< First atom index [device].
    int* j = nullptr;         ///< Second atom index [device].
    real_t* dx = nullptr;     ///< x-displacement (j - i) [Angstrom, device].
    real_t* dy = nullptr;     ///< y-displacement (j - i) [Angstrom, device].
    real_t* dz = nullptr;     ///< z-displacement (j - i) [Angstrom, device].
    real_t* dist = nullptr;   ///< Pair distance [Angstrom, device].
    int* n_pairs = nullptr;   ///< Current pair count [device, 1 element].
    int* capacity = nullptr;  ///< Pre-allocated slot count [device, 1 element].
};

/// @brief Bond topology list -- all arrays on device, molcrafts naming.
struct BondList {
    int* atomi = nullptr;    ///< First atom index [device].
    int* atomj = nullptr;    ///< Second atom index [device].
    int* type = nullptr;     ///< Bond type index [device].
    int* n_bonds = nullptr;  ///< Bond count [device, 1 element].
};

/// @brief Angle topology list -- all arrays on device, molcrafts naming.
struct AngleList {
    int* atomi = nullptr;     ///< First atom index [device].
    int* atomj = nullptr;     ///< Center atom index [device].
    int* atomk = nullptr;     ///< Third atom index [device].
    int* type = nullptr;      ///< Angle type index [device].
    int* n_angles = nullptr;  ///< Angle count [device, 1 element].
};

/// @brief Dihedral topology list -- all arrays on device, molcrafts naming.
struct DihedralList {
    int* atomi = nullptr;        ///< First atom index [device].
    int* atomj = nullptr;        ///< Second atom index [device].
    int* atomk = nullptr;        ///< Third atom index [device].
    int* atoml = nullptr;        ///< Fourth atom index [device].
    int* type = nullptr;         ///< Dihedral type index [device].
    int* n_dihedrals = nullptr;  ///< Dihedral count [device, 1 element].
};

/// @brief GPU-side system state -- device pointers + host-cached metadata.
///
/// Device pointer fields are GPU-resident. Metadata fields (n_atoms_, n_bonds_,
/// etc.) are host-cached counts set at materialize time to avoid D2H transfers.
///
/// Created by cpu::to_cuda(), freed by cuda::destroy().
struct Frame {
    // ── Per-atom arrays (device pointers, molcrafts "atoms" block) ──

    /// @brief Per-atom data block matching molpy atoms namespace.
    struct Atoms {
        real_t* x = nullptr;      ///< x positions [Angstrom, device].
        real_t* y = nullptr;      ///< y positions [Angstrom, device].
        real_t* z = nullptr;      ///< z positions [Angstrom, device].
        real_t* vx = nullptr;     ///< x velocities [Angstrom/fs, device].
        real_t* vy = nullptr;     ///< y velocities [Angstrom/fs, device].
        real_t* vz = nullptr;     ///< z velocities [Angstrom/fs, device].
        real_t* fx = nullptr;     ///< x forces [eV/Angstrom, device].
        real_t* fy = nullptr;     ///< y forces [eV/Angstrom, device].
        real_t* fz = nullptr;     ///< z forces [eV/Angstrom, device].
        real_t* mass = nullptr;   ///< Masses [amu, device].
        int* type = nullptr;      ///< Atom types [device].
        int* d_size = nullptr;    ///< Atom count [device, 1 element].
    } atoms;

    // ── System scalars (device pointers) ──

    prec_t* energy = nullptr;                     ///< Potential energy [eV, device, 1 element].
    prec_t* kinetic = nullptr;                    ///< Kinetic energy [eV, device, 1 element].
    prec_t* stress = nullptr;                     ///< Stress tensor [eV/Angstrom^3, device, 6 elements, Voigt].
    prec_t* box = nullptr;                        ///< Simulation cell 3x3 row-major [Angstrom, device, 9 elements].

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
        prec_t* H0    = nullptr;  ///< Base Hamiltonian (toy) [device, n*n].
        prec_t* H     = nullptr;  ///< Effective Hamiltonian / Fock matrix [device, n*n].
        prec_t* P     = nullptr;  ///< Density matrix [device, n*n].
        prec_t* P_old = nullptr;  ///< Previous density matrix [device, n*n].
        prec_t* C     = nullptr;  ///< Eigenvectors / MO coefficients [device, n*n, col-major].
        prec_t* eig   = nullptr;  ///< Eigenvalues / orbital energies [device, n].
        prec_t* work  = nullptr;  ///< Scratch buffer [device, n*n].
        int     n     = 0;        ///< Matrix dimension (host-cached).
        int     n_occ = 0;        ///< Occupation count (host-cached).
        bool    converged = false; ///< Convergence flag (host-cached).
        int     iteration = 0;     ///< Iteration that converged (host-cached).
        double  residual  = 0;     ///< Last convergence residual (host-cached).

        // ── RHF-specific ──
        prec_t* Hcore   = nullptr; ///< Core Hamiltonian T+V [device, n*n].
        prec_t* X       = nullptr; ///< Orthogonalization matrix S^{-1/2} [device, n*n].
        prec_t* F_prime = nullptr; ///< Transformed Fock X^T*F*X [device, n*n].
        prec_t* ERI     = nullptr; ///< Two-electron integrals [device, n^4].
        double  V_nuc   = 0;      ///< Nuclear repulsion energy (host-cached).

        // ── Basis data for GPU-native integral kernels ──
        prec_t* d_shell_alpha  = nullptr; ///< Exponents [device, n_shells * n_prim].
        prec_t* d_shell_coeff  = nullptr; ///< Normalized coefficients [device, n_shells * n_prim].
        prec_t* d_shell_center = nullptr; ///< Shell centers [device, n_shells * 3] (Bohr).
        int*    d_atom_Z       = nullptr; ///< Nuclear charges [device, natm].
        prec_t* d_atom_pos     = nullptr; ///< Atom positions [device, natm * 3] (Bohr).
        int     n_prim = 0;              ///< Primitives per shell (host-cached).
        int     n_shells = 0;            ///< Number of shells (host-cached).
        int     natm = 0;                ///< Number of atoms (host-cached).
    };
    SCFWorkspace scf;                             ///< SCF workspace (allocated by SCF components).

    // ── Host-cached metadata (set at materialize, avoids D2H in loop) ──

    /// @brief Number of atoms cached on the host side.
    /// @return Atom count.
    int n_atoms() const { return n_atoms_; }
    /// @brief Number of bonds cached on the host side.
    /// @return Bond count.
    int n_bonds() const { return n_bonds_; }
    /// @brief Number of angles cached on the host side.
    /// @return Angle count.
    int n_angles() const { return n_angles_; }
    /// @brief Number of dihedrals cached on the host side.
    /// @return Dihedral count.
    int n_dihedrals() const { return n_dihedrals_; }

    int n_atoms_ = 0;      ///< Host-cached atom count.
    int n_bonds_ = 0;      ///< Host-cached bond count.
    int n_angles_ = 0;     ///< Host-cached angle count.
    int n_dihedrals_ = 0;  ///< Host-cached dihedral count.
};

} // namespace cuda
} // namespace ATV
