# Acceptance — const-pot-bo-parallel-plate

## Build (binding)

MD-only configure (HF/PW/GPW/CUDA/Torch **OFF**):

```bash
bash -lc '
  module load buildenv-gcc/2026.03-mpich
  export OPENBLAS_ROOT="$OPENBLASDIR"
  export RUSTUP_TOOLCHAIN=stable-x86_64-unknown-linux-gnu
  export LD_LIBRARY_PATH="$LIBRARY_PATH:$LD_LIBRARY_PATH"
  cmake -B build -DCMAKE_BUILD_TYPE=Release \
    -DATV_BUILD_MD=ON -DATV_BUILD_KSPACE=ON \
    -DATV_BUILD_TESTS=ON -DATV_BUILD_REGRESSIONS=ON \
    -DATV_BUILD_HF=OFF -DATV_BUILD_PW=OFF -DATV_BUILD_GPW=OFF \
    -DATV_BUILD_CUDA=OFF -DATV_BUILD_TORCH=OFF -DBLA_VENDOR=OpenBLAS
  cmake --build build -j$(nproc) --target atomiverse reg_const-pot-bo-parallel-plate
'
```

## Oracle

| Field | Value |
| --- | --- |
| oracle | `lammps` |
| oracle_version | `3 Nov 2022` ELECTRODE madelung example log + plate_cap Madelung literature |
| pin path | `regressions/md/const-pot-bo-parallel-plate/REFERENCE.md` |
| live third-party in ctest | **forbidden** |

## Observable and tolerances

| ID | Check | Tol |
| --- | --- | --- |
| ac-001 | Public API is `cpu::md::ConstPotBO` Component; **no** `*Config` / `make_*` | grep / compile |
| ac-002 | Solves on `pre_force`; writes `kAtom_charge` for electrode types only | unit + e2e |
| ac-003 | Electroneutrality \(\lvert\sum_i Q_i\rvert \le 1\times10^{-10}\) e | e2e |
| ac-004 | Madelung twin: \(\lvert Q_{\mathrm{top}}\rvert\) vs LAMMPS offline pin \(0.279214485147238\) e | **rel ≤ 1%** (v1 pure-real/2D-image \(A\) vs Ewald electrode; documented) |
| ac-005 | \(\mathrm{sign}(Q_{\mathrm{top}})=\mathrm{sign}(\Psi_{\mathrm{top}}-\Psi_{\mathrm{bot}})\) for \(\Delta\Psi>0\) | e2e |
| ac-006 | Residual \(\lVert A Q - (\Psi-b)\rVert_\infty\) (unconstrained residual after symm proj) finite; projected residual small | \(\le 1\times10^{-8}\) on projected system or documented |
| ac-007 | MD-only build flags as above | configure log |
| ac-008 | Dual-gate: correctness pass + (`within_10x` or `perf_exempt`) | REFERENCE.md |

## Perf

No ELECTRODE-enabled `lmp` on worker PATH → **perf_exempt=true**, reason recorded.
When a twin binary exists, fair wall on the same 4-atom madelung geometry; ratio ≤ 10.

## Out of scope (must not block)

Freeze atoms, charging-work meta, CG solve, full Ewald electrode matrix, multi-GPU, SCF.
