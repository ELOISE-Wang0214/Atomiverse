/// @file cuda/transfer.cu
/// @brief cpu::to_cuda(), cuda::to_cpu(), cuda::destroy() implementations.
///
/// Implements double<->real_t/prec_t conversion with temporary buffers
/// when the device storage type differs from the host Real type.

#include "atomiverse/cpu/transfer.h"
#include "atomiverse/cuda/transfer.cuh"
#include "atomiverse/cuda/types.cuh"

#include <cuda_runtime.h>
#include <vector>

namespace ATV {

using cuda::prec_t;
using cuda::real_t;
using cuda::work_t;

// ── Narrowing helpers ────────────────────────────────────────────────
//
// cpu::Frame stores arrays as Real (double by default).
// cuda::Frame stores per-atom arrays as real_t (float by default).
// When real_t != Real, we must convert through a temporary buffer.

namespace {

/// @brief Copy host double array to device real_t array with narrowing.
void copy_h2d_narrow(real_t* d_dst, const Real* h_src, int n) {
    if constexpr (std::is_same_v<real_t, Real>) {
        cudaMemcpy(d_dst, h_src, n * sizeof(real_t), cudaMemcpyHostToDevice);
    } else {
        std::vector<real_t> tmp(n);
        for (int i = 0; i < n; ++i)
            tmp[i] = static_cast<real_t>(h_src[i]);
        cudaMemcpy(d_dst, tmp.data(), n * sizeof(real_t),
                   cudaMemcpyHostToDevice);
    }
}

/// @brief Copy device real_t array to host double array with widening.
void copy_d2h_widen(Real* h_dst, const real_t* d_src, int n) {
    if constexpr (std::is_same_v<real_t, Real>) {
        cudaMemcpy(h_dst, d_src, n * sizeof(real_t), cudaMemcpyDeviceToHost);
    } else {
        std::vector<real_t> tmp(n);
        cudaMemcpy(tmp.data(), d_src, n * sizeof(real_t),
                   cudaMemcpyDeviceToHost);
        for (int i = 0; i < n; ++i)
            h_dst[i] = static_cast<Real>(tmp[i]);
    }
}

/// @brief Copy host double array to device prec_t array with narrowing.
void copy_h2d_prec(prec_t* d_dst, const Real* h_src, int n) {
    if constexpr (std::is_same_v<prec_t, Real>) {
        cudaMemcpy(d_dst, h_src, n * sizeof(prec_t), cudaMemcpyHostToDevice);
    } else {
        std::vector<prec_t> tmp(n);
        for (int i = 0; i < n; ++i)
            tmp[i] = static_cast<prec_t>(h_src[i]);
        cudaMemcpy(d_dst, tmp.data(), n * sizeof(prec_t),
                   cudaMemcpyHostToDevice);
    }
}

/// @brief Copy device prec_t array to host double array with widening.
void copy_d2h_prec(Real* h_dst, const prec_t* d_src, int n) {
    if constexpr (std::is_same_v<prec_t, Real>) {
        cudaMemcpy(h_dst, d_src, n * sizeof(prec_t), cudaMemcpyDeviceToHost);
    } else {
        std::vector<prec_t> tmp(n);
        cudaMemcpy(tmp.data(), d_src, n * sizeof(prec_t),
                   cudaMemcpyDeviceToHost);
        for (int i = 0; i < n; ++i)
            h_dst[i] = static_cast<Real>(tmp[i]);
    }
}

} // anonymous namespace

// ── cpu::to_cuda ──────────────────────────────────────────────────────

namespace cpu {

cuda::Frame to_cuda(const Frame& frame) {
    int n = frame.atoms.size();
    cuda::Frame d{};

    // Host-cached metadata
    d.n_atoms_ = n;
    d.n_bonds_ = frame.bonds.size();
    d.n_angles_ = frame.angles.size();
    d.n_dihedrals_ = frame.dihedrals.size();

    // Scalar metadata on device
    cudaMalloc(&d.atoms.d_size, sizeof(int));
    cudaMemcpy(d.atoms.d_size, &n, sizeof(int), cudaMemcpyHostToDevice);

    // Per-atom arrays (real_t storage)
    cudaMalloc(&d.atoms.x, n * sizeof(real_t));
    cudaMalloc(&d.atoms.y, n * sizeof(real_t));
    cudaMalloc(&d.atoms.z, n * sizeof(real_t));
    cudaMalloc(&d.atoms.vx, n * sizeof(real_t));
    cudaMalloc(&d.atoms.vy, n * sizeof(real_t));
    cudaMalloc(&d.atoms.vz, n * sizeof(real_t));
    cudaMalloc(&d.atoms.fx, n * sizeof(real_t));
    cudaMalloc(&d.atoms.fy, n * sizeof(real_t));
    cudaMalloc(&d.atoms.fz, n * sizeof(real_t));
    cudaMalloc(&d.atoms.mass, n * sizeof(real_t));
    cudaMalloc(&d.atoms.type, n * sizeof(int));

    // Precision-critical scalars (prec_t storage)
    cudaMalloc(&d.energy, sizeof(prec_t));
    cudaMalloc(&d.kinetic, sizeof(prec_t));

    // Stress tensor (6 Voigt components, prec_t)
    cudaMalloc(&d.stress, 6 * sizeof(prec_t));
    cudaMemset(d.stress, 0, 6 * sizeof(prec_t));

    // H2D with narrowing conversion for per-atom arrays
    copy_h2d_narrow(d.atoms.x, frame.atoms.x.data(), n);
    copy_h2d_narrow(d.atoms.y, frame.atoms.y.data(), n);
    copy_h2d_narrow(d.atoms.z, frame.atoms.z.data(), n);
    copy_h2d_narrow(d.atoms.vx, frame.atoms.vx.data(), n);
    copy_h2d_narrow(d.atoms.vy, frame.atoms.vy.data(), n);
    copy_h2d_narrow(d.atoms.vz, frame.atoms.vz.data(), n);
    copy_h2d_narrow(d.atoms.mass, frame.atoms.mass.data(), n);
    cudaMemcpy(d.atoms.type, frame.atoms.type.data(), n * sizeof(int),
               cudaMemcpyHostToDevice);

    cudaMemset(d.atoms.fx, 0, n * sizeof(real_t));
    cudaMemset(d.atoms.fy, 0, n * sizeof(real_t));
    cudaMemset(d.atoms.fz, 0, n * sizeof(real_t));
    prec_t zero = 0;
    cudaMemcpy(d.energy, &zero, sizeof(prec_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d.kinetic, &zero, sizeof(prec_t), cudaMemcpyHostToDevice);

    // Box: 3x3 row-major on device (prec_t)
    cudaMalloc(&d.box, 9 * sizeof(prec_t));
    prec_t box_data[9];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            box_data[i * 3 + j] = static_cast<prec_t>(frame.box(i, j));
    cudaMemcpy(d.box, box_data, 9 * sizeof(prec_t), cudaMemcpyHostToDevice);

    // Pair list -- empty (allocated by NeighborList on first build)
    d.pairs = {};

    // Bond topology
    int nb = d.n_bonds_;
    if (nb > 0) {
        cudaMalloc(&d.bonds.n_bonds, sizeof(int));
        cudaMemcpy(d.bonds.n_bonds, &nb, sizeof(int), cudaMemcpyHostToDevice);
        cudaMalloc(&d.bonds.atomi, nb * sizeof(int));
        cudaMalloc(&d.bonds.atomj, nb * sizeof(int));
        cudaMalloc(&d.bonds.type, nb * sizeof(int));
        cudaMemcpy(d.bonds.atomi, frame.bonds.atomi.data(), nb * sizeof(int),
                   cudaMemcpyHostToDevice);
        cudaMemcpy(d.bonds.atomj, frame.bonds.atomj.data(), nb * sizeof(int),
                   cudaMemcpyHostToDevice);
        cudaMemcpy(d.bonds.type, frame.bonds.type.data(), nb * sizeof(int),
                   cudaMemcpyHostToDevice);
    }

    // Angle topology
    int na = d.n_angles_;
    if (na > 0) {
        cudaMalloc(&d.angles.n_angles, sizeof(int));
        cudaMemcpy(d.angles.n_angles, &na, sizeof(int), cudaMemcpyHostToDevice);
        cudaMalloc(&d.angles.atomi, na * sizeof(int));
        cudaMalloc(&d.angles.atomj, na * sizeof(int));
        cudaMalloc(&d.angles.atomk, na * sizeof(int));
        cudaMalloc(&d.angles.type, na * sizeof(int));
        cudaMemcpy(d.angles.atomi, frame.angles.atomi.data(), na * sizeof(int),
                   cudaMemcpyHostToDevice);
        cudaMemcpy(d.angles.atomj, frame.angles.atomj.data(), na * sizeof(int),
                   cudaMemcpyHostToDevice);
        cudaMemcpy(d.angles.atomk, frame.angles.atomk.data(), na * sizeof(int),
                   cudaMemcpyHostToDevice);
        cudaMemcpy(d.angles.type, frame.angles.type.data(), na * sizeof(int),
                   cudaMemcpyHostToDevice);
    }

    // Dihedral topology
    int nd = d.n_dihedrals_;
    if (nd > 0) {
        cudaMalloc(&d.dihedrals.n_dihedrals, sizeof(int));
        cudaMemcpy(d.dihedrals.n_dihedrals, &nd, sizeof(int),
                   cudaMemcpyHostToDevice);
        cudaMalloc(&d.dihedrals.atomi, nd * sizeof(int));
        cudaMalloc(&d.dihedrals.atomj, nd * sizeof(int));
        cudaMalloc(&d.dihedrals.atomk, nd * sizeof(int));
        cudaMalloc(&d.dihedrals.atoml, nd * sizeof(int));
        cudaMalloc(&d.dihedrals.type, nd * sizeof(int));
        cudaMemcpy(d.dihedrals.atomi, frame.dihedrals.atomi.data(),
                   nd * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d.dihedrals.atomj, frame.dihedrals.atomj.data(),
                   nd * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d.dihedrals.atomk, frame.dihedrals.atomk.data(),
                   nd * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d.dihedrals.atoml, frame.dihedrals.atoml.data(),
                   nd * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d.dihedrals.type, frame.dihedrals.type.data(),
                   nd * sizeof(int), cudaMemcpyHostToDevice);
    }

    return d;
}

} // namespace cpu

// ── cuda::to_cpu / cuda::destroy ──────────────────────────────────────

namespace cuda {

void to_cpu(const Frame& d, cpu::Frame& frame) {
    int n = frame.atoms.size();

    // D2H with widening conversion for per-atom arrays
    copy_d2h_widen(frame.atoms.x.data(), d.atoms.x, n);
    copy_d2h_widen(frame.atoms.y.data(), d.atoms.y, n);
    copy_d2h_widen(frame.atoms.z.data(), d.atoms.z, n);
    copy_d2h_widen(frame.atoms.vx.data(), d.atoms.vx, n);
    copy_d2h_widen(frame.atoms.vy.data(), d.atoms.vy, n);
    copy_d2h_widen(frame.atoms.vz.data(), d.atoms.vz, n);
    copy_d2h_widen(frame.atoms.fx.data(), d.atoms.fx, n);
    copy_d2h_widen(frame.atoms.fy.data(), d.atoms.fy, n);
    copy_d2h_widen(frame.atoms.fz.data(), d.atoms.fz, n);

    // Precision-critical scalars
    copy_d2h_prec(&frame.scalars.energy, d.energy, 1);
    copy_d2h_prec(&frame.scalars.kinetic, d.kinetic, 1);
    if (d.stress)
        copy_d2h_prec(frame.scalars.stress, d.stress, 6);

    // Box (may have changed via Deform or Barostat)
    if (d.box) {
        prec_t box_data[9];
        cudaMemcpy(box_data, d.box, 9 * sizeof(prec_t), cudaMemcpyDeviceToHost);
        Mat3 m;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                m(i, j) = static_cast<Real>(box_data[i * 3 + j]);
        frame.box.set_cell(m);
    }
}

void destroy(Frame& d) {
    cudaFree(d.atoms.x);
    cudaFree(d.atoms.y);
    cudaFree(d.atoms.z);
    cudaFree(d.atoms.vx);
    cudaFree(d.atoms.vy);
    cudaFree(d.atoms.vz);
    cudaFree(d.atoms.fx);
    cudaFree(d.atoms.fy);
    cudaFree(d.atoms.fz);
    cudaFree(d.atoms.mass);
    cudaFree(d.atoms.type);
    cudaFree(d.energy);
    cudaFree(d.kinetic);
    cudaFree(d.stress);
    cudaFree(d.box);
    cudaFree(d.atoms.d_size);

    cudaFree(d.pairs.i);
    cudaFree(d.pairs.j);
    cudaFree(d.pairs.dx);
    cudaFree(d.pairs.dy);
    cudaFree(d.pairs.dz);
    cudaFree(d.pairs.dist);
    cudaFree(d.pairs.n_pairs);
    cudaFree(d.pairs.capacity);

    cudaFree(d.bonds.atomi);
    cudaFree(d.bonds.atomj);
    cudaFree(d.bonds.type);
    cudaFree(d.bonds.n_bonds);

    cudaFree(d.angles.atomi);
    cudaFree(d.angles.atomj);
    cudaFree(d.angles.atomk);
    cudaFree(d.angles.type);
    cudaFree(d.angles.n_angles);

    cudaFree(d.dihedrals.atomi);
    cudaFree(d.dihedrals.atomj);
    cudaFree(d.dihedrals.atomk);
    cudaFree(d.dihedrals.atoml);
    cudaFree(d.dihedrals.type);
    cudaFree(d.dihedrals.n_dihedrals);

    // SCF workspace (shared)
    cudaFree(d.scf.H0);
    cudaFree(d.scf.H);
    cudaFree(d.scf.P);
    cudaFree(d.scf.P_old);
    cudaFree(d.scf.C);
    cudaFree(d.scf.eig);
    cudaFree(d.scf.work);
    // SCF workspace (RHF)
    cudaFree(d.scf.Hcore);
    cudaFree(d.scf.X);
    cudaFree(d.scf.F_prime);
    cudaFree(d.scf.ERI);
    // SCF workspace (GPU integral basis data)
    cudaFree(d.scf.d_shell_alpha);
    cudaFree(d.scf.d_shell_coeff);
    cudaFree(d.scf.d_shell_center);
    cudaFree(d.scf.d_atom_Z);
    cudaFree(d.scf.d_atom_pos);

    d = {};
}

} // namespace cuda
} // namespace ATV
