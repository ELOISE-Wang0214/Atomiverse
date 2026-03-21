/// @file cuda/transfer.cu
/// @brief cpu::to_cuda(), cuda::to_cpu(), cuda::destroy() implementations.

#include "atomiverse/cpu/transfer.h"
#include "atomiverse/cuda/transfer.cuh"

#include <cuda_runtime.h>

namespace ATV {

// ── cpu::to_cuda ──────────────────────────────────────────────────────

namespace cpu {

cuda::Frame to_cuda(const Frame& frame) {
    Int n = frame.n_atoms();
    cuda::Frame d{};

    // Host-cached metadata
    d.n_atoms_ = n;
    d.n_bonds_ = static_cast<Int>(frame.bond_i.size());
    d.n_angles_ = static_cast<Int>(frame.angle_i.size());
    d.n_dihedrals_ = static_cast<Int>(frame.dihedral_i.size());

    // Scalar metadata on device
    cudaMalloc(&d.d_n_atoms, sizeof(Int));
    cudaMemcpy(d.d_n_atoms, &n, sizeof(Int), cudaMemcpyHostToDevice);

    // Per-atom arrays
    cudaMalloc(&d.x, n * sizeof(Real));
    cudaMalloc(&d.y, n * sizeof(Real));
    cudaMalloc(&d.z, n * sizeof(Real));
    cudaMalloc(&d.vx, n * sizeof(Real));
    cudaMalloc(&d.vy, n * sizeof(Real));
    cudaMalloc(&d.vz, n * sizeof(Real));
    cudaMalloc(&d.fx, n * sizeof(Real));
    cudaMalloc(&d.fy, n * sizeof(Real));
    cudaMalloc(&d.fz, n * sizeof(Real));
    cudaMalloc(&d.mass, n * sizeof(Real));
    cudaMalloc(&d.type_id, n * sizeof(Int));
    cudaMalloc(&d.energy, sizeof(Real));
    cudaMalloc(&d.kinetic, sizeof(Real));

    // Stress tensor (6 Voigt components)
    cudaMalloc(&d.stress, 6 * sizeof(Real));
    cudaMemset(d.stress, 0, 6 * sizeof(Real));

    cudaMemcpy(d.x, frame.x.data(), n * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(d.y, frame.y.data(), n * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(d.z, frame.z.data(), n * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(d.vx, frame.vx.data(), n * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(d.vy, frame.vy.data(), n * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(d.vz, frame.vz.data(), n * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(d.mass, frame.mass.data(), n * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(d.type_id, frame.type_id.data(), n * sizeof(Int), cudaMemcpyHostToDevice);

    cudaMemset(d.fx, 0, n * sizeof(Real));
    cudaMemset(d.fy, 0, n * sizeof(Real));
    cudaMemset(d.fz, 0, n * sizeof(Real));
    Real zero = 0;
    cudaMemcpy(d.energy, &zero, sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(d.kinetic, &zero, sizeof(Real), cudaMemcpyHostToDevice);

    // Box: 3x3 row-major on device
    cudaMalloc(&d.box, 9 * sizeof(Real));
    Real box_data[9];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            box_data[i * 3 + j] = frame.box(i, j);
    cudaMemcpy(d.box, box_data, 9 * sizeof(Real), cudaMemcpyHostToDevice);

    // Pair list -- empty (allocated by NeighborList on first build)
    d.pairs = {};

    // Bond topology
    Int nb = d.n_bonds_;
    if (nb > 0) {
        cudaMalloc(&d.bonds.n_bonds, sizeof(Int));
        cudaMemcpy(d.bonds.n_bonds, &nb, sizeof(Int), cudaMemcpyHostToDevice);
        cudaMalloc(&d.bonds.i, nb * sizeof(Int));
        cudaMalloc(&d.bonds.j, nb * sizeof(Int));
        cudaMalloc(&d.bonds.type, nb * sizeof(Int));
        cudaMemcpy(d.bonds.i, frame.bond_i.data(), nb * sizeof(Int), cudaMemcpyHostToDevice);
        cudaMemcpy(d.bonds.j, frame.bond_j.data(), nb * sizeof(Int), cudaMemcpyHostToDevice);
        cudaMemcpy(d.bonds.type, frame.bond_type.data(), nb * sizeof(Int), cudaMemcpyHostToDevice);
    }

    // Angle topology
    Int na = d.n_angles_;
    if (na > 0) {
        cudaMalloc(&d.angles.n_angles, sizeof(Int));
        cudaMemcpy(d.angles.n_angles, &na, sizeof(Int), cudaMemcpyHostToDevice);
        cudaMalloc(&d.angles.i, na * sizeof(Int));
        cudaMalloc(&d.angles.j, na * sizeof(Int));
        cudaMalloc(&d.angles.k, na * sizeof(Int));
        cudaMalloc(&d.angles.type, na * sizeof(Int));
        cudaMemcpy(d.angles.i, frame.angle_i.data(), na * sizeof(Int), cudaMemcpyHostToDevice);
        cudaMemcpy(d.angles.j, frame.angle_j.data(), na * sizeof(Int), cudaMemcpyHostToDevice);
        cudaMemcpy(d.angles.k, frame.angle_k.data(), na * sizeof(Int), cudaMemcpyHostToDevice);
        cudaMemcpy(d.angles.type, frame.angle_type.data(), na * sizeof(Int), cudaMemcpyHostToDevice);
    }

    // Dihedral topology
    Int nd = d.n_dihedrals_;
    if (nd > 0) {
        cudaMalloc(&d.dihedrals.n_dihedrals, sizeof(Int));
        cudaMemcpy(d.dihedrals.n_dihedrals, &nd, sizeof(Int), cudaMemcpyHostToDevice);
        cudaMalloc(&d.dihedrals.i, nd * sizeof(Int));
        cudaMalloc(&d.dihedrals.j, nd * sizeof(Int));
        cudaMalloc(&d.dihedrals.k, nd * sizeof(Int));
        cudaMalloc(&d.dihedrals.l, nd * sizeof(Int));
        cudaMalloc(&d.dihedrals.type, nd * sizeof(Int));
        cudaMemcpy(d.dihedrals.i, frame.dihedral_i.data(), nd * sizeof(Int), cudaMemcpyHostToDevice);
        cudaMemcpy(d.dihedrals.j, frame.dihedral_j.data(), nd * sizeof(Int), cudaMemcpyHostToDevice);
        cudaMemcpy(d.dihedrals.k, frame.dihedral_k.data(), nd * sizeof(Int), cudaMemcpyHostToDevice);
        cudaMemcpy(d.dihedrals.l, frame.dihedral_l.data(), nd * sizeof(Int), cudaMemcpyHostToDevice);
        cudaMemcpy(d.dihedrals.type, frame.dihedral_type.data(), nd * sizeof(Int), cudaMemcpyHostToDevice);
    }

    return d;
}

} // namespace cpu

// ── cuda::to_cpu / cuda::destroy ──────────────────────────────────────

namespace cuda {

void to_cpu(const Frame& d, cpu::Frame& frame) {
    Int n = frame.n_atoms();
    cudaMemcpy(frame.x.data(), d.x, n * sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(frame.y.data(), d.y, n * sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(frame.z.data(), d.z, n * sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(frame.vx.data(), d.vx, n * sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(frame.vy.data(), d.vy, n * sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(frame.vz.data(), d.vz, n * sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(frame.fx.data(), d.fx, n * sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(frame.fy.data(), d.fy, n * sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(frame.fz.data(), d.fz, n * sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(&frame.energy, d.energy, sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(&frame.kinetic, d.kinetic, sizeof(Real), cudaMemcpyDeviceToHost);
    if (d.stress)
        cudaMemcpy(frame.stress, d.stress, 6 * sizeof(Real), cudaMemcpyDeviceToHost);

    // Box (may have changed via Deform or Barostat)
    if (d.box) {
        Real box_data[9];
        cudaMemcpy(box_data, d.box, 9 * sizeof(Real), cudaMemcpyDeviceToHost);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                frame.box(i, j) = box_data[i * 3 + j];
    }
}

void destroy(Frame& d) {
    cudaFree(d.x); cudaFree(d.y); cudaFree(d.z);
    cudaFree(d.vx); cudaFree(d.vy); cudaFree(d.vz);
    cudaFree(d.fx); cudaFree(d.fy); cudaFree(d.fz);
    cudaFree(d.mass); cudaFree(d.type_id);
    cudaFree(d.energy); cudaFree(d.kinetic);
    cudaFree(d.stress);
    cudaFree(d.box); cudaFree(d.d_n_atoms);

    cudaFree(d.pairs.i); cudaFree(d.pairs.j);
    cudaFree(d.pairs.dx); cudaFree(d.pairs.dy);
    cudaFree(d.pairs.dz); cudaFree(d.pairs.dist);
    cudaFree(d.pairs.n_pairs); cudaFree(d.pairs.capacity);

    cudaFree(d.bonds.i); cudaFree(d.bonds.j);
    cudaFree(d.bonds.type); cudaFree(d.bonds.n_bonds);

    cudaFree(d.angles.i); cudaFree(d.angles.j); cudaFree(d.angles.k);
    cudaFree(d.angles.type); cudaFree(d.angles.n_angles);

    cudaFree(d.dihedrals.i); cudaFree(d.dihedrals.j);
    cudaFree(d.dihedrals.k); cudaFree(d.dihedrals.l);
    cudaFree(d.dihedrals.type); cudaFree(d.dihedrals.n_dihedrals);

    // SCF workspace (shared)
    cudaFree(d.scf.H0); cudaFree(d.scf.H);
    cudaFree(d.scf.P); cudaFree(d.scf.P_old);
    cudaFree(d.scf.C); cudaFree(d.scf.eig);
    cudaFree(d.scf.work);
    // SCF workspace (RHF)
    cudaFree(d.scf.Hcore); cudaFree(d.scf.X);
    cudaFree(d.scf.F_prime); cudaFree(d.scf.ERI);
    // SCF workspace (GPU integral basis data)
    cudaFree(d.scf.d_shell_alpha); cudaFree(d.scf.d_shell_coeff);
    cudaFree(d.scf.d_shell_center);
    cudaFree(d.scf.d_atom_Z); cudaFree(d.scf.d_atom_pos);

    d = {};
}

} // namespace cuda
} // namespace ATV
