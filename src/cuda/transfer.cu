/// @file cuda/transfer.cu
/// @brief cpu::to_cuda(), cuda::to_cpu(), cuda::destroy() implementations.

#include "atomiverse/cpu/transfer.h"
#include "atomiverse/cuda/transfer.cuh"

#include <cuda_runtime.h>

namespace ATV {

// ── cpu::to_cuda ──────────────────────────────────────────────────────

namespace cpu {

cuda::Frame to_cuda(const Frame& frame) {
    Int n = frame.atoms.size();
    cuda::Frame d{};

    // Host-cached metadata
    d.n_atoms_ = n;
    d.n_bonds_ = frame.bonds.size();
    d.n_angles_ = frame.angles.size();
    d.n_dihedrals_ = frame.dihedrals.size();

    // Scalar metadata on device
    cudaMalloc(&d.atoms.d_size, sizeof(Int));
    cudaMemcpy(d.atoms.d_size, &n, sizeof(Int), cudaMemcpyHostToDevice);

    // Per-atom arrays
    cudaMalloc(&d.atoms.x, n * sizeof(Real));
    cudaMalloc(&d.atoms.y, n * sizeof(Real));
    cudaMalloc(&d.atoms.z, n * sizeof(Real));
    cudaMalloc(&d.atoms.vx, n * sizeof(Real));
    cudaMalloc(&d.atoms.vy, n * sizeof(Real));
    cudaMalloc(&d.atoms.vz, n * sizeof(Real));
    cudaMalloc(&d.atoms.fx, n * sizeof(Real));
    cudaMalloc(&d.atoms.fy, n * sizeof(Real));
    cudaMalloc(&d.atoms.fz, n * sizeof(Real));
    cudaMalloc(&d.atoms.mass, n * sizeof(Real));
    cudaMalloc(&d.atoms.type, n * sizeof(Int));
    cudaMalloc(&d.energy, sizeof(Real));
    cudaMalloc(&d.kinetic, sizeof(Real));

    // Stress tensor (6 Voigt components)
    cudaMalloc(&d.stress, 6 * sizeof(Real));
    cudaMemset(d.stress, 0, 6 * sizeof(Real));

    cudaMemcpy(d.atoms.x, frame.atoms.x.data(), n * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(d.atoms.y, frame.atoms.y.data(), n * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(d.atoms.z, frame.atoms.z.data(), n * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(d.atoms.vx, frame.atoms.vx.data(), n * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(d.atoms.vy, frame.atoms.vy.data(), n * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(d.atoms.vz, frame.atoms.vz.data(), n * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(d.atoms.mass, frame.atoms.mass.data(), n * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(d.atoms.type, frame.atoms.type.data(), n * sizeof(Int), cudaMemcpyHostToDevice);

    cudaMemset(d.atoms.fx, 0, n * sizeof(Real));
    cudaMemset(d.atoms.fy, 0, n * sizeof(Real));
    cudaMemset(d.atoms.fz, 0, n * sizeof(Real));
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
        cudaMalloc(&d.bonds.atomi, nb * sizeof(Int));
        cudaMalloc(&d.bonds.atomj, nb * sizeof(Int));
        cudaMalloc(&d.bonds.type, nb * sizeof(Int));
        cudaMemcpy(d.bonds.atomi, frame.bonds.atomi.data(), nb * sizeof(Int), cudaMemcpyHostToDevice);
        cudaMemcpy(d.bonds.atomj, frame.bonds.atomj.data(), nb * sizeof(Int), cudaMemcpyHostToDevice);
        cudaMemcpy(d.bonds.type, frame.bonds.type.data(), nb * sizeof(Int), cudaMemcpyHostToDevice);
    }

    // Angle topology
    Int na = d.n_angles_;
    if (na > 0) {
        cudaMalloc(&d.angles.n_angles, sizeof(Int));
        cudaMemcpy(d.angles.n_angles, &na, sizeof(Int), cudaMemcpyHostToDevice);
        cudaMalloc(&d.angles.atomi, na * sizeof(Int));
        cudaMalloc(&d.angles.atomj, na * sizeof(Int));
        cudaMalloc(&d.angles.atomk, na * sizeof(Int));
        cudaMalloc(&d.angles.type, na * sizeof(Int));
        cudaMemcpy(d.angles.atomi, frame.angles.atomi.data(), na * sizeof(Int), cudaMemcpyHostToDevice);
        cudaMemcpy(d.angles.atomj, frame.angles.atomj.data(), na * sizeof(Int), cudaMemcpyHostToDevice);
        cudaMemcpy(d.angles.atomk, frame.angles.atomk.data(), na * sizeof(Int), cudaMemcpyHostToDevice);
        cudaMemcpy(d.angles.type, frame.angles.type.data(), na * sizeof(Int), cudaMemcpyHostToDevice);
    }

    // Dihedral topology
    Int nd = d.n_dihedrals_;
    if (nd > 0) {
        cudaMalloc(&d.dihedrals.n_dihedrals, sizeof(Int));
        cudaMemcpy(d.dihedrals.n_dihedrals, &nd, sizeof(Int), cudaMemcpyHostToDevice);
        cudaMalloc(&d.dihedrals.atomi, nd * sizeof(Int));
        cudaMalloc(&d.dihedrals.atomj, nd * sizeof(Int));
        cudaMalloc(&d.dihedrals.atomk, nd * sizeof(Int));
        cudaMalloc(&d.dihedrals.atoml, nd * sizeof(Int));
        cudaMalloc(&d.dihedrals.type, nd * sizeof(Int));
        cudaMemcpy(d.dihedrals.atomi, frame.dihedrals.atomi.data(), nd * sizeof(Int), cudaMemcpyHostToDevice);
        cudaMemcpy(d.dihedrals.atomj, frame.dihedrals.atomj.data(), nd * sizeof(Int), cudaMemcpyHostToDevice);
        cudaMemcpy(d.dihedrals.atomk, frame.dihedrals.atomk.data(), nd * sizeof(Int), cudaMemcpyHostToDevice);
        cudaMemcpy(d.dihedrals.atoml, frame.dihedrals.atoml.data(), nd * sizeof(Int), cudaMemcpyHostToDevice);
        cudaMemcpy(d.dihedrals.type, frame.dihedrals.type.data(), nd * sizeof(Int), cudaMemcpyHostToDevice);
    }

    return d;
}

} // namespace cpu

// ── cuda::to_cpu / cuda::destroy ──────────────────────────────────────

namespace cuda {

void to_cpu(const Frame& d, cpu::Frame& frame) {
    Int n = frame.atoms.size();
    cudaMemcpy(frame.atoms.x.data(), d.atoms.x, n * sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(frame.atoms.y.data(), d.atoms.y, n * sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(frame.atoms.z.data(), d.atoms.z, n * sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(frame.atoms.vx.data(), d.atoms.vx, n * sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(frame.atoms.vy.data(), d.atoms.vy, n * sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(frame.atoms.vz.data(), d.atoms.vz, n * sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(frame.atoms.fx.data(), d.atoms.fx, n * sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(frame.atoms.fy.data(), d.atoms.fy, n * sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(frame.atoms.fz.data(), d.atoms.fz, n * sizeof(Real), cudaMemcpyDeviceToHost);
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
    cudaFree(d.atoms.x); cudaFree(d.atoms.y); cudaFree(d.atoms.z);
    cudaFree(d.atoms.vx); cudaFree(d.atoms.vy); cudaFree(d.atoms.vz);
    cudaFree(d.atoms.fx); cudaFree(d.atoms.fy); cudaFree(d.atoms.fz);
    cudaFree(d.atoms.mass); cudaFree(d.atoms.type);
    cudaFree(d.energy); cudaFree(d.kinetic);
    cudaFree(d.stress);
    cudaFree(d.box); cudaFree(d.atoms.d_size);

    cudaFree(d.pairs.i); cudaFree(d.pairs.j);
    cudaFree(d.pairs.dx); cudaFree(d.pairs.dy);
    cudaFree(d.pairs.dz); cudaFree(d.pairs.dist);
    cudaFree(d.pairs.n_pairs); cudaFree(d.pairs.capacity);

    cudaFree(d.bonds.atomi); cudaFree(d.bonds.atomj);
    cudaFree(d.bonds.type); cudaFree(d.bonds.n_bonds);

    cudaFree(d.angles.atomi); cudaFree(d.angles.atomj); cudaFree(d.angles.atomk);
    cudaFree(d.angles.type); cudaFree(d.angles.n_angles);

    cudaFree(d.dihedrals.atomi); cudaFree(d.dihedrals.atomj);
    cudaFree(d.dihedrals.atomk); cudaFree(d.dihedrals.atoml);
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
