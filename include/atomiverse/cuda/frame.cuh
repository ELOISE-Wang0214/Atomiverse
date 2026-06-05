/// @file cuda/frame.cuh
/// @brief cuda::Frame -- GPU-side system state (SoA device pointers).
///
/// ALL array data lives on the GPU. Host-cached metadata (n_atoms, n_bonds,
/// etc.) avoids D2H transfers during the simulation loop.
/// CUDA-free header (contains only pointer declarations).
///
/// Data is organized into sub-structs (atoms, bonds, angles, dihedrals)
/// matching the molcrafts naming convention (molpy Block model). The field
/// set mirrors `cpu::Frame` and `torch::Frame` exactly -- the three
/// execution frames differ only in device representation (xtensor /
/// raw device pointers / at::Tensor), never in which physical quantities
/// they hold. Per the CLAUDE.md Frame purity rule the frame carries only
/// fundamental physical quantities: method-specific scratch (the SCF
/// workspace, the neighbor pair list) is owned by the components that
/// compute it, not by the frame.
///
/// @see cpu::Frame for the host-side counterpart
/// @see core/system_scalars.h for the shared system-scalar bundle
/// @see cuda/scf/scf_workspace.cuh for the relocated SCF scratch
/// @see cuda/neighbor_list.cuh for the relocated neighbor pair list
/// @see cpu::to_cuda / cuda::to_cpu for transfer
#pragma once

#include <cstdint>

#include "atomiverse/core/system_scalars.h"
#include "atomiverse/cuda/types.cuh"

namespace ATV {
namespace cuda {

/// @brief Bond topology list -- all arrays on device, molcrafts naming.
struct BondList {
    int* atomi = nullptr;   ///< First atom index [device].
    int* atomj = nullptr;   ///< Second atom index [device].
    int* type = nullptr;    ///< Bond type index [device].
    int* n_bonds = nullptr; ///< Bond count [device, 1 element].
};

/// @brief Angle topology list -- all arrays on device, molcrafts naming.
struct AngleList {
    int* atomi = nullptr;    ///< First atom index [device].
    int* atomj = nullptr;    ///< Center atom index [device].
    int* atomk = nullptr;    ///< Third atom index [device].
    int* type = nullptr;     ///< Angle type index [device].
    int* n_angles = nullptr; ///< Angle count [device, 1 element].
};

/// @brief Dihedral topology list -- all arrays on device, molcrafts naming.
struct DihedralList {
    int* atomi = nullptr;       ///< First atom index [device].
    int* atomj = nullptr;       ///< Second atom index [device].
    int* atomk = nullptr;       ///< Third atom index [device].
    int* atoml = nullptr;       ///< Fourth atom index [device].
    int* type = nullptr;        ///< Dihedral type index [device].
    int* n_dihedrals = nullptr; ///< Dihedral count [device, 1 element].
};

/// @brief GPU-side system state -- device pointers + host-cached metadata.
///
/// Device pointer fields are GPU-resident. Metadata fields (n_atoms_,
/// n_bonds_, etc.) are host-cached counts set at materialize time to avoid
/// D2H transfers.
///
/// Created by cpu::to_cuda() / cuda::materialize(), freed by cuda::destroy().
struct Frame {
    // -- Per-atom arrays (device pointers, molcrafts "atoms" block) --

    /// @brief Per-atom data block matching molpy atoms namespace.
    struct Atoms {
        real_t* x = nullptr;     ///< x positions [Angstrom, device].
        real_t* y = nullptr;     ///< y positions [Angstrom, device].
        real_t* z = nullptr;     ///< z positions [Angstrom, device].
        real_t* vx = nullptr;    ///< x velocities [Angstrom/fs, device].
        real_t* vy = nullptr;    ///< y velocities [Angstrom/fs, device].
        real_t* vz = nullptr;    ///< z velocities [Angstrom/fs, device].
        real_t* fx = nullptr;    ///< x forces [eV/Angstrom, device].
        real_t* fy = nullptr;    ///< y forces [eV/Angstrom, device].
        real_t* fz = nullptr;    ///< z forces [eV/Angstrom, device].
        real_t* mass = nullptr;  ///< Masses [amu, device].
        real_t* charge = nullptr; ///< Partial charges [e, device].
        int* type = nullptr;     ///< Atom types [device].

        /// @brief Global atom identity [device], persistent across MPI
        /// migration and restart. Mirrors `cpu::Frame::Atoms::tag`.
        std::int64_t* tag = nullptr;

        /// @brief Per-atom PBC crossing counts, flat row-major [device,
        /// 3 * n_atoms], `image[i * 3 + d]` for atom i along dimension d.
        /// Mirrors `cpu::Frame::Atoms::image`.
        std::int8_t* image = nullptr;

        int* d_size = nullptr; ///< Atom count [device, 1 element].
    } atoms;

    // -- System scalars (device-resident SystemScalars) --

    /// @brief Per-frame system scalars (energy, kinetic, conserved, stress,
    /// kinetic stress) -- one device-resident `SystemScalars` struct.
    ///
    /// Embeds `core/system_scalars.h`'s `SystemScalars` exactly as
    /// `cpu::Frame` does; the only difference is the device representation
    /// (a raw device pointer here, an inline value on the CPU frame). See
    /// that header for the field-by-field documentation.
    SystemScalars* scalars = nullptr;

    prec_t* box =
        nullptr; ///< Simulation cell 3x3 row-major [Angstrom, device, 9].

    // -- Topology sub-structures --

    BondList bonds;     ///< Bond topology (transferred by cpu::to_cuda).
    AngleList angles;   ///< Angle topology (transferred by cpu::to_cuda).
    DihedralList
        dihedrals; ///< Dihedral topology (transferred by cpu::to_cuda).

    // -- Host-cached metadata (set at materialize, avoids D2H in loop) --

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

    int n_atoms_ = 0;     ///< Host-cached atom count.
    int n_bonds_ = 0;     ///< Host-cached bond count.
    int n_angles_ = 0;    ///< Host-cached angle count.
    int n_dihedrals_ = 0; ///< Host-cached dihedral count.

    // -- Domain-decomposition counts (single-process multi-GPU) --
    // Owned atoms occupy [0, n_local); halo ghosts received from peer GPUs
    // occupy [n_local, n_local + n_ghost); the per-atom arrays are allocated to
    // `capacity` (>= n_local + n_ghost) so ghosts can be appended each step
    // without reallocation. SINGLE-GPU invariant: n_local == n_atoms,
    // n_ghost == 0, capacity == n_atoms (set at materialize) — so existing
    // kernels that read n_atoms() are unchanged.

    /// Owned atoms this subdomain integrates (force kernels write [0,n_local)).
    int n_local() const { return n_local_; }
    /// Ghost atoms received from peer GPUs (read-only halo copies).
    int n_ghost() const { return n_ghost_; }
    /// Owned + ghost — the range force kernels must READ positions over.
    int n_all() const { return n_local_ + n_ghost_; }
    /// Allocated per-atom array length (headroom for ghosts).
    int capacity() const { return capacity_; }

    int n_local_ = 0;
    int n_ghost_ = 0;
    int capacity_ = 0;
};

} // namespace cuda
} // namespace ATV
