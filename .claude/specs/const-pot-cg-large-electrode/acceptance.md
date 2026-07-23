# Acceptance — const-pot-cg-large-electrode

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
  cmake --build build -j$(nproc) --target atomiverse \
        test_cpu_md_const_pot_cg reg_const-pot-cg-large-electrode
'
```

## Oracle

| Field | Value |
| --- | --- |
| oracle | `metalwalls` |
| oracle_version | MetalWalls-style CG electrode solve (map E2); dual twin = in-tree `ConstPotBO` dense path on same \(A\) |
| pin path | `regressions/md/const-pot-cg-large-electrode/REFERENCE.md` |
| live third-party in ctest | **forbidden** |

## Observable and tolerances

| ID | Check | Tol |
| --- | --- | --- |
| ac-001 | Public API is `cpu::md::ConstPotCG` Component; **no** `*Config` / `make_*` | grep / compile |
| ac-002 | Solves on `pre_force`; writes `kAtom_charge` for electrode types only | unit + e2e |
| ac-003 | Electroneutrality \(\lvert\sum_i Q_i\rvert \le 1\times10^{-10}\) e | e2e |
| ac-004 | **Large plate dual twin:** \(\max_i \lvert Q_i^{\mathrm{CG}}-Q_i^{\mathrm{BO}}\rvert \le 1\times10^{-8}\) e (same geometry, \(\eta\), \(\Psi\), images) | e2e |
| ac-005 | Group totals \(\lvert Q_{\mathrm{bot}}^{\mathrm{CG}}-Q_{\mathrm{bot}}^{\mathrm{BO}}\rvert\) and top | \(\le 1\times10^{-8}\) e |
| ac-006 | Small madelung primary-cell: CG \(Q_{\mathrm{top}}\) matches dense BO analytic pin \(0.27175356\) e | abs \(\le 1\times10^{-6}\) e |
| ac-007 | MD-only build flags as above | configure log |
| ac-008 | Dual-gate: correctness pass + (`within_10x` or `perf_exempt`) | REFERENCE.md |
| ac-009 | \(N_{\mathrm{elec}} \ge 18\) on the large-plate regression (many sites) | e2e |

## Perf

No MetalWalls runtime binary on worker PATH → **perf_exempt=true**, reason
recorded. Optional internal wall: \(T_{\mathrm{CG}}\) vs \(T_{\mathrm{BO}}\) on
same multi-site plate (informational, not the 10× oracle gate).

## Out of scope (must not block)

ASPC predictor, freeze, charging-work, full Ewald electrode matrix, PIM, TF,
CUDA.
