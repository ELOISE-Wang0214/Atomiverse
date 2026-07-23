# four-site-tip4p-foothold

**Branch:** `wf/four-site-tip4p-foothold`  
**Chapter:** E9 (molecular water / virtual sites)  
**Oracle:** `lammps` (mandatory — `pair_style lj/cut/tip4p/cut`)  
**Status:** open

## Goal

Ship a minimal **four-site TIP4P-class** water nonbonded foothold: Lennard-Jones on
oxygen + Coulomb on the off-oxygen M virtual site (and H), with force redistribution
to real atoms (Feenstra partitioning), as a public `cpu::` pair Potential Component.
Frame I/O only; no public `*Config` bags; no `make_*` / `build_*` factories.

## Non-goals

- Long-range TIP4P (`tip4p/long` + PPPM/tip4p)
- Flexible bonded FF calibration / GAFF water
- MetalWalls electrode + water coupling (later E9 follow-on)
- CUDA path
- SCF / PW / GPW

## Public API (binding)

```cpp
cpu::PairLJCutTip4pCut pair(
    /*type_O=*/0, /*type_H=*/1,
    /*qdist [Å]=*/0.15,
    /*r_oh [Å]=*/0.9572,
    /*theta_hoh [rad]=*/104.52 * π/180,
    /*eps_OO [eV]=*/..., /*sigma_OO [Å]=*/3.15365,
    /*cut_lj [Å]=*/10.0, /*cut_coul [Å]=*/10.0);
pair.set_special_bonds(0.0, 0.0, 0.0);  // water excludes 1-2/1-3/1-4
```

Constructor arguments establish all parameters (no `Tip4pConfig`).  
M-site geometry follows LAMMPS:

\[
\alpha = \frac{q_{\mathrm{dist}}}{\cos(\theta/2)\, r_{\mathrm{OH}}},\quad
\mathbf{r}_M = \mathbf{r}_O + \alpha\cdot\tfrac12\bigl((\mathbf{r}_{H1}-\mathbf{r}_O)+(\mathbf{r}_{H2}-\mathbf{r}_O)\bigr)
\]

H atoms of oxygen with tag `t` are tags `t+1` and `t+2` (LAMMPS tip4p convention).

## Tasks

1. [x] `PairLJCutTip4pCut` header + impl under `cpu/potentials/`
2. [x] Unit tests: M-site geometry, dimer energy vs pin, force Newton III, force pin
3. [x] Offline LAMMPS pin + `REFERENCE.md` (`oracle: lammps`)
4. [x] E2E regression `regressions/md/four-site-tip4p-foothold/`
5. [x] Wire CMake regression + unit test targets
6. [x] MD-only build; regression PASS
7. [x] Perf gate vs same LAMMPS input (ratio ≈ 0.34 ≤ 10)

## Acceptance

| ID | Criterion | Tol |
| -- | --------- | --- |
| A1 | Oracle is LAMMPS `lj/cut/tip4p/cut`; pin path documented offline | — |
| A2 | Offline pin files under `tools/metalwalls_oracle/inputs/lammps/four-site-tip4p-foothold/` | — |
| A3 | `REFERENCE.md` has `oracle: lammps` and `oracle_version:` | — |
| A4 | Dimer total energy matches pin within **1e-8 eV** (dual-gate energy) | 1e-8 eV |
| A5 | Dual-gate geometry: reconstructed OM distance = `qdist` within **1e-10 Å** | 1e-10 Å |
| A6 | Forces match LAMMPS pin within **1e-6 eV/Å** (per-component on all atoms) | 1e-6 eV/Å |
| A7 | MD-only cmake: HF/PW/GPW/CUDA/Torch **OFF** | — |
| A8 | Unit tests under `tests/` only; one e2e regression | — |
| A9 | Perf: \(T_{\mathrm{atv}}/T_{\mathrm{oracle}} \le 10\) on fair wall, or `perf_exempt` | ≤10× |

## Observable

- Total pair energy \(E\) [eV]
- M-site OM distance [Å]
- Forces on O/H [eV/Å]

## Build (binding)

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release \
  -DATV_BUILD_MD=ON -DATV_BUILD_KSPACE=ON \
  -DATV_BUILD_TESTS=ON -DATV_BUILD_REGRESSIONS=ON \
  -DATV_BUILD_HF=OFF -DATV_BUILD_PW=OFF -DATV_BUILD_GPW=OFF \
  -DATV_BUILD_CUDA=OFF -DATV_BUILD_TORCH=OFF \
  -DBLA_VENDOR=OpenBLAS
```
