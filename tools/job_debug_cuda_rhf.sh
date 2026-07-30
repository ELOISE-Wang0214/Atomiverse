#!/bin/bash -l
#SBATCH -p gpu
#SBATCH -A naiss2026-4-715-gpu
#SBATCH --gpus 1
#SBATCH -c 8
#SBATCH -t 00:30:00
#SBATCH -J atv-dbg-rhf
#SBATCH -o tools/debug_cuda_rhf_%j.out
#SBATCH -e tools/debug_cuda_rhf_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
module load GPU/buildenv-gcccuda/2026.03-cu13.0
export OPENBLAS_ROOT="${OPENBLASDIR:-}"
unset RUSTUP_TOOLCHAIN || true
export RUSTUP_TOOLCHAIN=stable-aarch64-unknown-linux-gnu
export LD_LIBRARY_PATH="${LIBRARY_PATH:-}:${LD_LIBRARY_PATH:-}"

BD=build-cuda-bench
if [[ ! -f "$BD/CMakeCache.txt" ]]; then
  cmake -B "$BD" -DCMAKE_BUILD_TYPE=Release \
    -DATV_BUILD_CUDA=ON -DATV_BUILD_HF=ON -DATV_BUILD_TESTS=ON \
    -DATV_BUILD_PW=OFF -DATV_USE_MPI=OFF -DBLA_VENDOR=OpenBLAS \
    -DATV_BUILD_REGRESSIONS=OFF -DCMAKE_CUDA_ARCHITECTURES=90
fi
cmake --build "$BD" -j"$(nproc)" --target atomiverse

cat > /tmp/dbg_rhf.cu <<'CU'
#include <cstdio>
#include <cmath>
#include <vector>
#include <cuda_runtime.h>
#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/units.h"
#include "atomiverse/core/gto.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/basis_data/accessors.h"
#include "atomiverse/cpu/scf/internal/integrals.h"
#include "atomiverse/cuda/frame/transfer.cuh"
#include "atomiverse/cpu/frame/transfer.h"
#include "atomiverse/cuda/scf/device_integrals.cuh"
#include "atomiverse/cuda/scf/rhf.cuh"

using namespace ATV;

int main() {
  cpu::Frame f(2);
  f.span(kAtom_x)[0]=0; f.span(kAtom_y)[0]=0; f.span(kAtom_z)[0]=0;
  f.span(kAtom_x)[1]=0.74; f.span(kAtom_y)[1]=0; f.span(kAtom_z)[1]=0;
  f.span(kAtom_type)[0]=1; f.span(kAtom_type)[1]=1;
  f.span(kAtom_mass)[0]=1.008; f.span(kAtom_mass)[1]=1.008;
  f.box = cpu::Box::cubic(100.0);

  // Host reference integrals
  cpu::scf::internal::IntegralWorkspace ref{f, cpu::scf::BasisSet::STO3G, 0, true};
  std::printf("HOST V_nuc=%.10f\n", ref.V_nuc);
  std::printf("HOST S00=%.10f S01=%.10f S11=%.10f\n",
    ref.S(0,0), ref.S(0,1), ref.S(1,1));
  std::printf("HOST H00=%.10f H01=%.10f H11=%.10f\n",
    ref.Hcore(0,0), ref.Hcore(0,1), ref.Hcore(1,1));

  // Device pure-s integrals (closed form)
  const auto& h1s = cpu::scf::basis_data::pure_s_sto3g_shell(1);
  const int np = (int)h1s.alpha.size();
  const int ns = 2;
  std::vector<double> ha(ns*np), hc(ns*np), hctr(ns*3), hpos(6);
  int hZ[2]={1,1};
  const double a2b = units::kAngstromInBohr;
  for (int s=0;s<ns;s++) {
    for (int p=0;p<np;p++) {
      ha[s*np+p]=h1s.alpha[p];
      hc[s*np+p]=h1s.coeff[p]*ATV::phys_norm_s(h1s.alpha[p]);
    }
    hctr[s*3+0]=(s==0?0.0:0.74*a2b); hctr[s*3+1]=0; hctr[s*3+2]=0;
    hpos[s*3+0]=hctr[s*3+0]; hpos[s*3+1]=0; hpos[s*3+2]=0;
  }
  double *dA,*dC,*dCtr,*dS,*dT,*dV,*dPos; int *dZ;
  cudaMalloc(&dA,ha.size()*8); cudaMalloc(&dC,hc.size()*8);
  cudaMalloc(&dCtr,ns*3*8); cudaMalloc(&dS,4*8); cudaMalloc(&dT,4*8);
  cudaMalloc(&dV,4*8); cudaMalloc(&dPos,6*8); cudaMalloc(&dZ,2*4);
  cudaMemcpy(dA,ha.data(),ha.size()*8,cudaMemcpyHostToDevice);
  cudaMemcpy(dC,hc.data(),hc.size()*8,cudaMemcpyHostToDevice);
  cudaMemcpy(dCtr,hctr.data(),ns*3*8,cudaMemcpyHostToDevice);
  cudaMemcpy(dPos,hpos.data(),6*8,cudaMemcpyHostToDevice);
  cudaMemcpy(dZ,hZ,8,cudaMemcpyHostToDevice);
  cuda::scf::launch_1e_integrals(dS,dT,dV,dA,dC,dCtr,dZ,dPos,ns,np,2,nullptr);
  cudaDeviceSynchronize();
  double hS[4],hT[4],hVmat[4];
  cudaMemcpy(hS,dS,32,cudaMemcpyDeviceToHost);
  cudaMemcpy(hT,dT,32,cudaMemcpyDeviceToHost);
  cudaMemcpy(hVmat,dV,32,cudaMemcpyDeviceToHost);
  std::printf("DEV_SS S00=%.10f S01=%.10f S11=%.10f\n", hS[0],hS[1],hS[3]);
  std::printf("DEV_SS T00=%.10f T01=%.10f\n", hT[0],hT[1]);
  std::printf("DEV_SS V00=%.10f V01=%.10f\n", hVmat[0],hVmat[1]);
  std::printf("DEV_SS H00=%.10f H01=%.10f\n", hT[0]+hVmat[0], hT[1]+hVmat[1]);

  // Full product path
  cuda::scf::RHF rhf{cpu::scf::BasisSet::STO3G, 0, 1e-10, 100};
  rhf.with_gradient(false);
  rhf.run(f);
  std::printf("PRODUCT E=%.10f Ha  conv=%d iters=%d\n",
    rhf.energy(), (int)rhf.converged(), rhf.iterations());
  std::printf("PRODUCT meta_eV=%.10f\n", (double)f.meta.get(kEnergy));

  // Host-fill CUDA path (upload matrices)
  {
    cpu::Frame f2=f;
    // force rebuild
    xt::xtensor<double,2> S,H; xt::xtensor<double,4> ERI; double Vn=0; int nb=0,ne=0;
    cpu::scf::internal::IntegralWorkspace::fill_rhf_matrices(
      f2, cpu::scf::BasisSet::STO3G, 0, true, S, H, ERI, Vn, nb, ne);
    std::printf("HOST_FILL V_nuc=%.10f n_basis=%d n_elec=%d\n", Vn, nb, ne);
    // Use RHF host path energy as oracle
    cpu::scf::RHF host{cpu::scf::BasisSet::STO3G, 0, 1e-10, 100};
    host.run(f2);
    std::printf("CPU_RHF E=%.10f Ha\n", host.energy());
  }
  return 0;
}
CU

# Link like the test binary: pull flags from existing target if possible
nvcc -O3 -std=c++20 -arch=sm_90 -DATV_HAS_HF -DATV_HAS_CUDA \
  -Iinclude -I"$BD" -I"$BD/_deps" \
  $(find "$BD/_deps" -type d -name include 2>/dev/null | sed 's/^/-I/' | tr '\n' ' ') \
  /tmp/dbg_rhf.cu -o /tmp/dbg_rhf \
  "$BD/libatomiverse.a" \
  -lcublas -lcusolver -lcurand -Xlinker -lgomp 2>&1 | tail -50

if [[ -x /tmp/dbg_rhf ]]; then
  /tmp/dbg_rhf
else
  echo "link failed"
  # fallback: try with more libs
  ls "$BD"/lib*.a "$BD"/**/lib*.a 2>/dev/null | head
fi
echo DONE
