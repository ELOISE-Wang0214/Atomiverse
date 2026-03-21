/// @file bench_rhf.cu
/// @brief RHF integral benchmarks: native, CUDA, and optionally libcint.

#ifdef ATV_HAS_RHF

#include <benchmark/benchmark.h>
#include <cmath>
#include <vector>

#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/native_integrals.h"
#include "atomiverse/core/constants.h"

#ifdef ATV_HAS_LIBCINT
extern "C" {
#include <cint.h>
#include <cint_funcs.h>
}
#undef atm
#undef bas
#endif // ATV_HAS_LIBCINT

using namespace ATV;
using namespace ATV::cpu::scf;

// ── Shared constants ──────────────────────────────────────────────────

static constexpr double kA2B = ATV::kAngstromToBohr;
static constexpr double sto3g_exp[3] = {3.42525091, 0.62391373, 0.16885540};
static constexpr double sto3g_coeff[3] = {0.15432897, 0.53532814, 0.44463454};

#ifdef ATV_HAS_LIBCINT
// ── libcint environment (reusable) ────────────────────────────────────

struct LibcintEnv {
    std::vector<FINT> atm, bas;
    std::vector<double> env;
    int natm, nbas;

    LibcintEnv() : natm(2), nbas(2) {
        int off = PTR_ENV_START;
        atm.resize(natm * ATM_SLOTS, 0);
        bas.resize(nbas * BAS_SLOTS, 0);
        env.resize(off + 12, 0.0);

        atm[CHARGE_OF] = 1; atm[PTR_COORD] = off;
        env[off] = 0; env[off+1] = 0; env[off+2] = 0; off += 3;
        atm[CHARGE_OF + ATM_SLOTS] = 1; atm[PTR_COORD + ATM_SLOTS] = off;
        env[off] = 0.74 * kA2B; env[off+1] = 0; env[off+2] = 0; off += 3;

        int pe = off;
        for (int p = 0; p < 3; ++p) env[off++] = sto3g_exp[p];
        int pc = off;
        for (int p = 0; p < 3; ++p) env[off++] = sto3g_coeff[p] * CINTgto_norm(0, sto3g_exp[p]);

        for (int s = 0; s < nbas; ++s) {
            bas[ATOM_OF + BAS_SLOTS*s] = s;
            bas[ANG_OF + BAS_SLOTS*s] = 0;
            bas[NPRIM_OF + BAS_SLOTS*s] = 3;
            bas[NCTR_OF + BAS_SLOTS*s] = 1;
            bas[PTR_EXP + BAS_SLOTS*s] = pe;
            bas[PTR_COEFF + BAS_SLOTS*s] = pc;
        }
    }
};
#endif // ATV_HAS_LIBCINT

// ── Native shell setup (reusable) ─────────────────────────────────────

struct NativeEnv {
    double norm_coeff[3];
    double centers[2][3];
    native::Shell shells[2];
    native::Atom atoms[2];

    NativeEnv() {
        for (int p = 0; p < 3; ++p)
            norm_coeff[p] = sto3g_coeff[p] * native::phys_norm(sto3g_exp[p]);
        centers[0][0] = 0; centers[0][1] = 0; centers[0][2] = 0;
        centers[1][0] = 0.74 * kA2B; centers[1][1] = 0; centers[1][2] = 0;
        for (int s = 0; s < 2; ++s) {
            shells[s] = {3, sto3g_exp, norm_coeff, {centers[s][0], centers[s][1], centers[s][2]}};
            atoms[s] = {1, {centers[s][0], centers[s][1], centers[s][2]}};
        }
    }
};

#ifdef ATV_HAS_LIBCINT
// ── GtoNorm ───────────────────────────────────────────────────────────

static void BM_GtoNorm_Libcint(benchmark::State& state) {
    for (auto _ : state)
        for (int p = 0; p < 3; ++p)
            benchmark::DoNotOptimize(CINTgto_norm(0, sto3g_exp[p]));
}
BENCHMARK(BM_GtoNorm_Libcint);
#endif // ATV_HAS_LIBCINT

static void BM_GtoNorm_Native(benchmark::State& state) {
    for (auto _ : state)
        for (int p = 0; p < 3; ++p)
            benchmark::DoNotOptimize(native::gto_norm(0, sto3g_exp[p]));
}
BENCHMARK(BM_GtoNorm_Native);

#ifdef ATV_HAS_LIBCINT
// ── Overlap ───────────────────────────────────────────────────────────

static void BM_Overlap_Libcint(benchmark::State& state) {
    LibcintEnv lc;
    for (auto _ : state) {
        FINT shls[2];
        double buf[1];
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j) {
                shls[0] = i; shls[1] = j;
                int1e_ovlp_sph(buf, nullptr, shls, lc.atm.data(), lc.natm,
                               lc.bas.data(), lc.nbas, lc.env.data(), nullptr, nullptr);
                benchmark::DoNotOptimize(buf[0]);
            }
    }
}
BENCHMARK(BM_Overlap_Libcint);
#endif // ATV_HAS_LIBCINT

static void BM_Overlap_Native(benchmark::State& state) {
    NativeEnv ne;
    for (auto _ : state) {
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                benchmark::DoNotOptimize(native::overlap_ss(ne.shells[i], ne.shells[j]));
    }
}
BENCHMARK(BM_Overlap_Native);

#ifdef ATV_HAS_LIBCINT
// ── Kinetic ───────────────────────────────────────────────────────────

static void BM_Kinetic_Libcint(benchmark::State& state) {
    LibcintEnv lc;
    for (auto _ : state) {
        FINT shls[2]; double buf[1];
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j) {
                shls[0] = i; shls[1] = j;
                int1e_kin_sph(buf, nullptr, shls, lc.atm.data(), lc.natm,
                             lc.bas.data(), lc.nbas, lc.env.data(), nullptr, nullptr);
                benchmark::DoNotOptimize(buf[0]);
            }
    }
}
BENCHMARK(BM_Kinetic_Libcint);
#endif // ATV_HAS_LIBCINT

static void BM_Kinetic_Native(benchmark::State& state) {
    NativeEnv ne;
    for (auto _ : state) {
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                benchmark::DoNotOptimize(native::kinetic_ss(ne.shells[i], ne.shells[j]));
    }
}
BENCHMARK(BM_Kinetic_Native);

#ifdef ATV_HAS_LIBCINT
// ── Nuclear ───────────────────────────────────────────────────────────

static void BM_Nuclear_Libcint(benchmark::State& state) {
    LibcintEnv lc;
    for (auto _ : state) {
        FINT shls[2]; double buf[1];
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j) {
                shls[0] = i; shls[1] = j;
                int1e_nuc_sph(buf, nullptr, shls, lc.atm.data(), lc.natm,
                             lc.bas.data(), lc.nbas, lc.env.data(), nullptr, nullptr);
                benchmark::DoNotOptimize(buf[0]);
            }
    }
}
BENCHMARK(BM_Nuclear_Libcint);
#endif // ATV_HAS_LIBCINT

static void BM_Nuclear_Native(benchmark::State& state) {
    NativeEnv ne;
    for (auto _ : state) {
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                benchmark::DoNotOptimize(native::nuclear_ss(ne.shells[i], ne.shells[j], ne.atoms, 2));
    }
}
BENCHMARK(BM_Nuclear_Native);

#ifdef ATV_HAS_LIBCINT
// ── ERI ───────────────────────────────────────────────────────────────

static void BM_ERI_Libcint(benchmark::State& state) {
    LibcintEnv lc;
    for (auto _ : state) {
        FINT shls[4]; double buf[1];
        for (int i = 0; i < 2; ++i)
          for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 2; ++k)
              for (int l = 0; l < 2; ++l) {
                shls[0]=i; shls[1]=j; shls[2]=k; shls[3]=l;
                cint2e_sph(buf, shls, lc.atm.data(), lc.natm,
                          lc.bas.data(), lc.nbas, lc.env.data(), nullptr);
                benchmark::DoNotOptimize(buf[0]);
              }
    }
}
BENCHMARK(BM_ERI_Libcint);
#endif // ATV_HAS_LIBCINT

static void BM_ERI_Native(benchmark::State& state) {
    NativeEnv ne;
    for (auto _ : state) {
        for (int i = 0; i < 2; ++i)
          for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 2; ++k)
              for (int l = 0; l < 2; ++l)
                benchmark::DoNotOptimize(
                    native::eri_ssss(ne.shells[i], ne.shells[j], ne.shells[k], ne.shells[l]));
    }
}
BENCHMARK(BM_ERI_Native);

// ── CUDA kernels ──────────────────────────────────────────────────────

#include "atomiverse/cuda/scf/device_integrals.cuh"

struct CudaBasisEnv {
    Real* d_alpha = nullptr;
    Real* d_coeff = nullptr;
    Real* d_centers = nullptr;
    Int* d_Z = nullptr;
    Real* d_pos = nullptr;
    Real* d_S = nullptr;
    Real* d_T = nullptr;
    Real* d_V = nullptr;
    Real* d_ERI = nullptr;

    CudaBasisEnv() {
        Real h_alpha[6], h_coeff[6], h_centers[6], h_pos[6];
        Int h_Z[2] = {1, 1};
        double cx[2] = {0, 0.74 * kA2B};
        for (int s = 0; s < 2; ++s) {
            for (int p = 0; p < 3; ++p) {
                h_alpha[s*3+p] = Real(sto3g_exp[p]);
                h_coeff[s*3+p] = Real(sto3g_coeff[p] * native::phys_norm(sto3g_exp[p]));
            }
            h_centers[s*3] = Real(cx[s]); h_centers[s*3+1] = 0; h_centers[s*3+2] = 0;
            h_pos[s*3] = Real(cx[s]); h_pos[s*3+1] = 0; h_pos[s*3+2] = 0;
        }
        cudaMalloc(&d_alpha, 6*sizeof(Real)); cudaMemcpy(d_alpha, h_alpha, 6*sizeof(Real), cudaMemcpyHostToDevice);
        cudaMalloc(&d_coeff, 6*sizeof(Real)); cudaMemcpy(d_coeff, h_coeff, 6*sizeof(Real), cudaMemcpyHostToDevice);
        cudaMalloc(&d_centers, 6*sizeof(Real)); cudaMemcpy(d_centers, h_centers, 6*sizeof(Real), cudaMemcpyHostToDevice);
        cudaMalloc(&d_Z, 2*sizeof(Int)); cudaMemcpy(d_Z, h_Z, 2*sizeof(Int), cudaMemcpyHostToDevice);
        cudaMalloc(&d_pos, 6*sizeof(Real)); cudaMemcpy(d_pos, h_pos, 6*sizeof(Real), cudaMemcpyHostToDevice);
        cudaMalloc(&d_S, 4*sizeof(Real)); cudaMalloc(&d_T, 4*sizeof(Real)); cudaMalloc(&d_V, 4*sizeof(Real));
        cudaMalloc(&d_ERI, 16*sizeof(Real));
    }
    ~CudaBasisEnv() {
        cudaFree(d_alpha); cudaFree(d_coeff); cudaFree(d_centers);
        cudaFree(d_Z); cudaFree(d_pos);
        cudaFree(d_S); cudaFree(d_T); cudaFree(d_V); cudaFree(d_ERI);
    }
};

static void BM_1e_CUDA_Kernel(benchmark::State& state) {
    CudaBasisEnv ce;
    // warmup
    cuda::scf::launch_1e_integrals(ce.d_S, ce.d_T, ce.d_V, ce.d_alpha, ce.d_coeff, ce.d_centers,
                                   ce.d_Z, ce.d_pos, 2, 3, 2, nullptr);
    cudaDeviceSynchronize();
    for (auto _ : state) {
        cuda::scf::launch_1e_integrals(ce.d_S, ce.d_T, ce.d_V, ce.d_alpha, ce.d_coeff, ce.d_centers,
                                       ce.d_Z, ce.d_pos, 2, 3, 2, nullptr);
        cudaDeviceSynchronize();
    }
}
BENCHMARK(BM_1e_CUDA_Kernel);

static void BM_ERI_CUDA_Kernel(benchmark::State& state) {
    CudaBasisEnv ce;
    cuda::scf::launch_eri(ce.d_ERI, ce.d_alpha, ce.d_coeff, ce.d_centers, 2, 3, nullptr);
    cudaDeviceSynchronize();
    for (auto _ : state) {
        cuda::scf::launch_eri(ce.d_ERI, ce.d_alpha, ce.d_coeff, ce.d_centers, 2, 3, nullptr);
        cudaDeviceSynchronize();
    }
}
BENCHMARK(BM_ERI_CUDA_Kernel);

// ── Full build ────────────────────────────────────────────────────────

static cpu::Frame make_h2() {
    cpu::Frame frame(2);
    frame.x[0] = Real(0);     frame.y[0] = Real(0); frame.z[0] = Real(0);
    frame.x[1] = Real(0.74);  frame.y[1] = Real(0); frame.z[1] = Real(0);
    frame.mass[0] = Real(1.008); frame.mass[1] = Real(1.008);
    frame.type_id[0] = 1; frame.type_id[1] = 1;
    frame.box = Mat3::cubic(100.0);
    return frame;
}

static void BM_FullBuild_Native(benchmark::State& state) {
    auto frame = make_h2();
    for (auto _ : state) {
        auto pack = build_integrals_sto3g(frame);
        benchmark::DoNotOptimize(pack.V_nuc);
    }
}
BENCHMARK(BM_FullBuild_Native)->Unit(benchmark::kMicrosecond);

// ── Stage-separated end-to-end benchmarks ────────────────────────────

#include "atomiverse/cpu/scf/rhf_solve.h"
#include "atomiverse/cuda/scf/scf.cuh"
#include "atomiverse/cuda/scf/rhf_fock.cuh"
#include "atomiverse/cuda/scf/rhf_density.cuh"

static void BM_CPURHFSolve(benchmark::State& state) {
    auto frame = make_h2();
    auto integrals = build_integrals_sto3g(frame);
    Real tol = (sizeof(Real) == sizeof(float)) ? Real(1e-5) : Real(1e-10);
    for (auto _ : state) {
        auto result = cpu::scf::rhf_solve(integrals, tol, 100);
        benchmark::DoNotOptimize(result.energy);
    }
}
BENCHMARK(BM_CPURHFSolve)->Unit(benchmark::kMicrosecond);

static void BM_CUDARHFSolve(benchmark::State& state) {
    auto frame = make_h2();
    auto integrals = build_integrals_sto3g(frame);
    Real tol = (sizeof(Real) == sizeof(float)) ? Real(1e-5) : Real(1e-10);
    for (auto _ : state) {
        // Includes materialize + solve + retrieve (full Driver lifecycle)
        scf::SCF solver;
        solver.add(scf::RHFFock{integrals});
        solver.add(scf::EigenSolver{});
        solver.add(scf::RHFDensity{integrals.n_occ()});
        solver.add(scf::ConvergenceCheck{tol, false});
        solver.build();
        solver.run(frame, 100);
        benchmark::DoNotOptimize(frame.energy);
    }
}
BENCHMARK(BM_CUDARHFSolve)->Unit(benchmark::kMicrosecond);

// ── LiH benchmarks ───────────────────────────────────────────────────

static cpu::Frame make_lih() {
    cpu::Frame frame(2);
    frame.x[0] = Real(0);     frame.y[0] = Real(0); frame.z[0] = Real(0);
    frame.x[1] = Real(1.595); frame.y[1] = Real(0); frame.z[1] = Real(0);
    frame.mass[0] = Real(6.941); frame.mass[1] = Real(1.008);
    frame.type_id[0] = 3; frame.type_id[1] = 1;
    frame.box = Mat3::cubic(100.0);
    return frame;
}

static void BM_CPURHFSolve_LiH(benchmark::State& state) {
    auto frame = make_lih();
    auto integrals = build_integrals_sto3g(frame);
    Real tol = (sizeof(Real) == sizeof(float)) ? Real(1e-5) : Real(1e-10);
    for (auto _ : state) {
        auto result = cpu::scf::rhf_solve(integrals, tol, 100);
        benchmark::DoNotOptimize(result.energy);
    }
}
BENCHMARK(BM_CPURHFSolve_LiH)->Unit(benchmark::kMicrosecond);

static void BM_CUDARHFSolve_LiH(benchmark::State& state) {
    auto frame = make_lih();
    auto integrals = build_integrals_sto3g(frame);
    Real tol = (sizeof(Real) == sizeof(float)) ? Real(1e-5) : Real(1e-10);
    for (auto _ : state) {
        scf::SCF solver;
        solver.add(scf::RHFFock{integrals});
        solver.add(scf::EigenSolver{});
        solver.add(scf::RHFDensity{integrals.n_occ()});
        solver.add(scf::ConvergenceCheck{tol, false});
        solver.build();
        solver.run(frame, 100);
        benchmark::DoNotOptimize(frame.energy);
    }
}
BENCHMARK(BM_CUDARHFSolve_LiH)->Unit(benchmark::kMicrosecond);

BENCHMARK_MAIN();

#else
int main() { return 0; }
#endif
