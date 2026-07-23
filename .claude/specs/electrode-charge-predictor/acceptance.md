# Acceptance — electrode-charge-predictor

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
    test_cpu_electrode_charge_aspc reg_electrode-charge-predictor
'
```

## Oracle

| Field | Value |
| --- | --- |
| oracle | `metalwalls` |
| oracle_version | algorithm: MetalWalls const-pot BO/CG + Kolafa ASPC (J. Comput. Chem. 25, 335 (2004); Coretti et al. JCP 157, 184801 (2022) / MetalWalls JOSS) |
| pin path | `regressions/md/electrode-charge-predictor/REFERENCE.md` |
| live third-party in ctest | **forbidden** |

## Observable and tolerances

| ID | Check | Tol |
| --- | --- | --- |
| ac-001 | Public `cpu::electrodes::ChargeASPC` + `ConstPotBO` CG/predictor members; **no** `*Config` / `make_*` | grep / compile |
| ac-002 | ASPC \(K=2\) Lagrange coeffs \(\omega=\{1,-3,3\}\) (oldest→newest) | unit exact |
| ac-003 | Linear history \(Q_n=n\): predict next = \(n+1\) (exact for order ≥1 polynomials of deg ≤K) | unit abs ≤ 1e-12 |
| ac-004 | CG (`use_cg=true`) \(Q\) matches dense Cholesky on fixed madelung-like RHS | abs ≤ 1e-8 e |
| ac-005 | Short MD / multi-RHS: mean `last_iters` with `predictor_order=2` **strictly less** than cold (`predictor_order=0`) after history fill | e2e |
| ac-006 | Electroneutrality \(\lvert\sum Q_i\rvert \le 1\times10^{-8}\) e on CG+ASPC path | e2e |
| ac-007 | MD-only build flags as above | configure log |
| ac-008 | Dual-gate: correctness pass + (`within_10x` or `perf_exempt`) | REFERENCE.md |

## Perf

No MetalWalls twin binary on PATH → **perf_exempt=true**, reason recorded.
When a twin exists, fair wall on the same short-MD geometry; ratio ≤ 10.

## Out of scope (must not block)

Preconditioned large-electrode CG product, piston, TF metallicity, SCF density ASPC.
