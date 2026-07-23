# thomas-fermi-electrode

**Branch:** `wf/thomas-fermi-electrode`  
**Chapter:** E7 (electrode metallicity — Thomas–Fermi finite screening)  
**Oracle:** `metalwalls` (MW-only physics; LAMMPS ELECTRODE is perfect conductor only)  
**Status:** acceptance_met (verified on wf/thomas-fermi-electrode)  
**Resume:** wave-2 assign; no prior branch/spec

## Goal

Ship **Thomas–Fermi (TF) finite metallicity** on the constant-potential BO electrode
solve: a single screening length \(l_{\mathrm{TF}}\) [Å] plus electrode site density
\(d\) [1/Å³] that adds the Scalfi / MetalWalls quadratic hardness to the elastance
matrix, reducing capacitance and total electrode charge \(Q\) relative to the perfect
metal (\(l_{\mathrm{TF}}=0\)).

Public API: members on `cpu::md::ConstPotBO` (Component on `cpu::md::MD` at
`pre_force`). Frame I/O only; **no** public `*Config` bags; **no** `make_*` /
`build_*` factories; **not** a new Driver category; charge solve **not** inside
`KSpace*`.

## Science (binding)

Scalfi, Dufils, Reeves, Rotenberg, Salanne, *J. Chem. Phys.* **153**, 174704
(2020), DOI [10.1063/5.0028232](https://doi.org/10.1063/5.0028232)
(arXiv:1910.13341). Implemented in MetalWalls.

Second-order TF kinetic energy (metal units \(e\), Å, eV):

\[
U_{\mathrm{TF}}^{(2)} = \frac{l_{\mathrm{TF}}^{2}\, d}{2\epsilon_0}\sum_i Q_i^{2}
\quad\Rightarrow\quad
A_{ii}\;\leftarrow\; A_{ii} + \frac{l_{\mathrm{TF}}^{2}\, d}{\epsilon_0}.
\]

With \(k_e = 1/(4\pi\epsilon_0) =\) `units::coulomb`:

\[
\frac{1}{\epsilon_0} = 4\pi\,k_e
\quad\Rightarrow\quad
A_{ii}\;\mathrel{+}= 4\pi\cdot\texttt{units::coulomb}\cdot l_{\mathrm{TF}}^{2}\cdot d.
\]

\(l_{\mathrm{TF}}=0\) recovers perfect-metal `ConstPotBO` bit-for-bit (same path as
`const-pot-bo-parallel-plate`). Electroneutrality projector unchanged (acts on the
hardened \(A\)).

## Non-goals

- CG / iterative electrode solve (E2 `const-pot-cg-large-electrode`)
- Full 2D Ewald electrode matrix beyond existing 2D Gaussian images
- Electrolyte structure response / long MD (only algebraic \(Q\), \(C\))
- CUDA path; SCF / PW / GPW
- AIM / DAIM / ML electrodes (E11–E12)

## Public API (binding)

```cpp
cpu::md::ConstPotBO conp(/*eta [1/Å]=*/1.0 / 0.56);  // MetalWalls η = 1/w
conp.add_electrode(/*type=*/1, /*Psi [V]=*/0.0);
conp.add_electrode(/*type=*/2, /*Psi [V]=*/1.0);
conp.electroneutrality = true;
conp.xy_image_range = 50;
conp.l_tf = 2.0;              // Thomas–Fermi length [Å]; 0 = perfect metal
conp.site_density = 0.0593;   // electrode site density d [1/Å³] (Au FCC ~4/a0³)
// … build / launch or MD add_ref …
double Q = conp.total_charge(1);           // [e]
double C = std::abs(Q) / conp.delta_psi(); // [e/V] integral capacitance
```

Parameters are **public members** (and ctor `eta`) — no `ThomasFermiConfig`.

## Tasks

1. [ ] SPEC (this file) + index row
2. [ ] Add `l_tf` [Å] and `site_density` [1/Å³] to `ConstPotBO`; apply TF diagonal in `build_matrix_`
3. [ ] Unit tests: (a) \(l_{\mathrm{TF}}=0\) matches existing primary-cell pin; (b) finite \(l_{\mathrm{TF}}\) reduces \(|Q|\); (c) analytic diagonal hardness on 2-site vacuum capacitor
4. [ ] Offline pin (literature / continuum + discrete offline script) under `tools/metalwalls_oracle/` + `REFERENCE.md` (`oracle: metalwalls`)
5. [ ] E2E regression `regressions/md/thomas-fermi-electrode/` — gold capacitor toy: \(C\), \(Q\) response
6. [ ] Wire CMake regression target
7. [ ] MD-only build; unit + regression PASS
8. [ ] Perf gate: MW binary twin unavailable → `perf_exempt` with reason

## Acceptance

| ID | Criterion | Tol |
| -- | --------- | --- |
| A1 | Oracle is MetalWalls TF model (not LAMMPS perfect conductor); pin path offline | — |
| A2 | `REFERENCE.md` has `oracle: metalwalls` and `oracle_version:` | — |
| A3 | Offline pin documented under `tools/metalwalls_oracle/` (inputs and/or pin script notes) | — |
| A4 | Gold-capacitor toy: \(Q_{\mathrm{top}}(l_{\mathrm{TF}}=0)\) and \(Q_{\mathrm{top}}(l_{\mathrm{TF}}>0)\) match offline pin within **1%** relative (dual-gate \(Q\)) | 1% rel |
| A5 | Dual-gate capacitance: \(C(l_{\mathrm{TF}}>0) < C(l_{\mathrm{TF}}=0)\) and \(C = \|Q\|/\Delta\Psi\) matches pin within **1%** | 1% rel |
| A6 | \(l_{\mathrm{TF}}=0\) path does not change existing `const-pot-bo` primary-cell unit pin | 1e-12 |
| A7 | MD-only cmake: HF/PW/GPW/CUDA/Torch **OFF** | — |
| A8 | Unit tests under `tests/` only; one e2e under `regressions/md/thomas-fermi-electrode/` | — |
| A9 | Perf: \(\le 10\times\) vs MW twin, or `perf_exempt` if no MW binary | ≤10× or exempt |

## Observable

- Total electrode charge \(Q\) [e] on each wall
- Integral capacitance \(C = |Q| / \Delta\Psi\) [e/V] (equivalently F when converted)
- Sign: \(\mathrm{sign}(Q_{\mathrm{top}})=+1\) for \(\Delta\Psi = +1\) V

## Geometry (gold capacitor toy)

Empty parallel-plate capacitor, gold-like site density \(d = 4 / a_0^3\) with
\(a_0 = 4.07\) Å (\(d \approx 0.0593\) Å⁻³). Minimal electrode stacks (few atomic
planes), vacuum gap, \(\Delta\Psi = 1\) V, Gaussian width matching MetalWalls-class
\(\eta^{-1} \approx 0.56\) Å. No electrolyte (isolates metallicity → \(C,Q\)).

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

## Perf

MetalWalls TF twin binary is not required on the login node. If unavailable:
`perf_exempt: true` with reason. Status green requires `within_10x` or
`perf_exempt`.
