# Acceptance — `piston-electrode-pz`

Binding checklist. `acceptance_met=true` only when **every** box is verified.

## Build (binding — MD-only)

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
        test_cpu_md_piston_electrode reg_piston-electrode-pz
'
```

- [ ] HF/PW/GPW/CUDA/Torch **OFF**; MD+KSPACE+TESTS+REGRESSIONS **ON**

## Oracle and pin

- [ ] `oracle: metalwalls` in `regressions/md/piston-electrode-pz/REFERENCE.md`
- [ ] `oracle_version:` recorded (MetalWalls algorithm / JOSS citation when no binary)
- [ ] Offline pin path documented under `tools/metalwalls_oracle/` and/or REFERENCE
- [ ] No live MetalWalls / LAMMPS launch at ctest time

## Observable and tolerances

| ID | Check | Tol |
| --- | --- | --- |
| ac-001 | Public API is `cpu::md::PistonElectrode` Component; **no** `*Config` / `make_*` | compile / grep |
| ac-002 | Rigid \(z\) motion: all piston atoms share \(\Delta z\) and \(v_z\); \(v_x=v_y=0\) | unit |
| ac-003 | \(F_{\mathrm{ext}} = \pm P_{\mathrm{target}} A\) with \(A=L_x L_y\); sign by `Side` | unit |
| ac-004 | Instantaneous \(P_z = \pm F_z^{\mathrm{atoms}}/A\) matches sum of wall \(F_z\) / area | unit ≤ 1e-12 rel |
| ac-005 | Force balance: with \(P_{\mathrm{target}} = P_z^{\mathrm{atoms}}\) at fixed geometry, \(\lvert F_{\mathrm{net}}\rvert / (\lvert F_z^{\mathrm{atoms}}\rvert+\varepsilon) \le 1\times10^{-10}\) | unit + e2e |
| ac-006 | E2E parallel plate: production mean \(\langle P_z\rangle\) within **10%** of \(P_{\mathrm{target}}\) (or absolute pin in REFERENCE) after damped equilibration | regression |
| ac-007 | Lower wall (if frozen) stays fixed; piston type is **not** also `FixFreeze` | e2e |
| ac-008 | MD-only build flags as above | configure |
| ac-009 | Dual-gate: correctness pass + (`within_10x` **or** `perf_exempt`) | REFERENCE |

## Architecture

- [ ] Component on `cpu::md::MD` — not a new Driver category
- [ ] No public `*Config`; no `make_*` / `build_*`
- [ ] Frame I/O only; charge solve not inside `KSpace*`
- [ ] Parameters = ctor args / members (`type_id`, `P_target`, `Side`, optional `tau_damp`)

## Tests

- [ ] Unit test `tests/cpu/dynamics/test_piston_electrode.cpp` passes
- [ ] Regression `reg_piston-electrode-pz` passes vs offline pin

## Perf

- [ ] Fair wall vs MetalWalls twin if binary available; ratio ≤ 10.0
- [ ] Else `perf_exempt=true` with reason in REFERENCE (expected on worker: no MW binary)

## Delivery

- [ ] Branch `wf/piston-electrode-pz`; never commit on `dev`; never push
- [ ] Commit includes spec, component, unit test, regression, REFERENCE, CMake wire
