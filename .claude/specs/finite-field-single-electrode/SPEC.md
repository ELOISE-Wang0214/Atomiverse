# finite-field-single-electrode

**Branch:** `wf/finite-field-single-electrode`  
**Chapter:** E3 (external field + electrode; Dufils single-electrode 3D)  
**Oracle:** `metalwalls` (MW-only physics) — single electrode + finite field under
3D PBC is **not** the two-electrode LAMMPS `ffield` path (which requires
exactly two electrodes). Offline pin = analytic Dufils/Reed Gaussian BO with
external-field potential term (same kernels as MetalWalls / ConstPotBO).  
**Status:** acceptance met (2026-07-24)

## Goal

Ship a minimal **finite-field + single electrode (3D PBC)** foothold: one
constant-potential electrode wall whose BO charges respond to a uniform
external electric field, composed with `FixEfield` for electrolyte (and
electrode) forces, plus a classical polarization observer. Frame I/O only;
Components on `cpu::md::MD`; no public `*Config` bags; no `make_*` / `build_*`.

## Non-goals

- Two-electrode LAMMPS `ffield` twin (ΔΨ/Lz internal efield) — separate
- CG / ASPC electrode predictors (other wave-2 slugs)
- Thomas–Fermi metallicity, PIM joint solve
- CUDA path; SCF / PW / GPW
- Full Ewald electrode matrix (E2 residual)

## Physics (binding)

Dufils *et al.*, Phys. Rev. Lett. **123**, 195501 (2019): finite electric field
on a **single** constant-potential electrode under 3D PBC.

Electrode BO (Gaussian elastance, same as `ConstPotBO`):

\[
(AQ)_i = \Psi_i - b_i + \mathbf{E}\cdot\mathbf{r}_i
\]

with \(\Phi_{\mathrm{ext}} = -\mathbf{E}\cdot\mathbf{r}\) so
\(\Psi = (AQ) + b + \Phi_{\mathrm{ext}}\). Electrolyte contribution \(b_i\)
unchanged from E1. Self term
\(A_{ii} = k_e\,\eta\sqrt{2/\pi}\), off-diagonal
\(A_{ij}=k_e\,\mathrm{erf}(\eta_{\mathrm{eff}} r)/r\).

Single-electrode free charge: **electroneutrality off** by default for the
single-wall foothold (total \(Q\) is the free charge that screens \(\mathbf{E}\);
two-electrode `symm on` remains available when two groups are registered).

Forces on charged atoms: compose existing `cpu::md::FixEfield{Ex,Ey,Ez}` at
`post_force` (metal units \(F=qE\), \(U=-\sum q\mathbf{E}\cdot\mathbf{r}^{\mathrm{u}}\)).

Polarization (observable):

\[
\mathbf{P} = \frac{1}{V}\sum_i q_i\,\mathbf{r}_i^{\mathrm{u}}
\quad [e/\mathrm{Å}^2]
\]

## Public API (binding)

```cpp
// Finite-field BO on a single electrode wall (members, no Config bag)
cpu::md::ConstPotBO conp(/*eta [1/Å]=*/2.0);
conp.add_electrode(/*type=*/1, /*Psi [V]=*/0.0);  // one wall
conp.Ex = 0.0; conp.Ey = 0.0; conp.Ez = 0.1;     // [V/Å]
conp.electroneutrality = false;                   // free total Q
// optional: conp.xy_image_range = 0;               // primary-cell foothold

// Forces: existing FixEfield (dep ✓)
cpu::md::FixEfield field(/*Ex=*/0.0, /*Ey=*/0.0, /*Ez=*/0.1);

// Polarization observer (classical)
cpu::md::ClassicalPolarization pol;  // last_P() [e/Å²]

cpu::md::MD sim;
sim.add(std::move(conp));
sim.add(std::move(field));
sim.add(std::move(pol));
```

No new Driver category; charge solve stays outside `KSpace*`.

## Tasks

1. [x] Extend `ConstPotBO` with public `Ex,Ey,Ez` [V/Å] in the BO RHS
2. [x] `ClassicalPolarization` Observer: \(P = V^{-1}\sum q r^u\) [e/Å²]
3. [x] Unit tests: vacuum single-site \(Q\) vs analytic; field sign; pol formula
4. [x] Offline pin + `REFERENCE.md` (`oracle: metalwalls` / literature analytic)
5. [x] E2E regression `regressions/md/finite-field-single-electrode/`
6. [x] Wire CMake regression + unit test targets
7. [x] MD-only build; regression PASS
8. [x] Perf gate: `perf_exempt` if no MetalWalls twin binary

## Acceptance

| ID | Criterion | Tol |
| -- | --------- | --- |
| A1 | Oracle name documented: metalwalls / literature analytic Dufils+Reed | — |
| A2 | Offline pin path under `tools/metalwalls_oracle/` or REFERENCE formulas | — |
| A3 | `REFERENCE.md` has `oracle:` and `oracle_version:` lines | — |
| A4 | Vacuum single-site \(Q\) matches analytic pin | \(\lvert Q-Q_{\mathrm{pin}}\rvert\le 10^{-10}\) e |
| A5 | Dual-gate polarization: \(P_z = Q z / V\) matches pin | \(\lvert P_z-P_{\mathrm{pin}}\rvert\le 10^{-12}\) e/Å² |
| A6 | Field forces: electrolyte \(F_z = q E_z\) via `FixEfield` | \(10^{-12}\) eV/Å |
| A7 | MD-only cmake: HF/PW/GPW/CUDA/Torch **OFF** | — |
| A8 | Unit tests under `tests/` only; one e2e regression | — |
| A9 | Perf: ≤10× vs oracle binary, or `perf_exempt` with reason | ≤10× |

## Observable

- Electrode total charge \(Q\) [e]
- Polarization \(\mathbf{P}\) [e/Å²]
- Field forces on electrolyte [eV/Å]

## Build (binding)

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
        -DATV_BUILD_CUDA=OFF -DATV_BUILD_TORCH=OFF \
        -DBLA_VENDOR=OpenBLAS
  cmake --build build -j$(nproc)
'
```
