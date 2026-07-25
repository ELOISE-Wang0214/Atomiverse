---
slug: rsh-eri-cam-b3lyp
status: approved
created: 2026-07-24
revised: 2026-07-24
domain: scf / dft
queue: pyscf-engine-gaps
deferred: true
depends_on: []
unblocks: [cam-b3lyp]
grilled: false
---

# Spec: `rsh-eri-cam-b3lyp` — erf-attenuated ERI + CAM-B3LYP foothold

**Slug:** `rsh-eri-cam-b3lyp`  
**Domain:** Gaussian-basis SCF / DFT (CPU)  
**Queue:** PySCF engine gaps (deferred — not auto-impl)  
**Status:** approved, **implement later**  
**Unblocks queue row:** `cam-b3lyp`  
**Sibling:** `periodic-rhf-gamma` (independent engine; no dep)

## Summary

Ship the missing **range-separated hybrid (RSH) Coulomb engine** and a
product foothold for **CAM-B3LYP** on CPU. Today `rys_eri_block` only
implements the bare \(1/r\) operator (documented as libcint-equivalent
*non*-range-separated Coulomb). Global hybrids (B3LYP / PBE0) reuse that
operator with a single scalar \(\alpha_\mathrm{HF}\); CAM-B3LYP needs a
partitioned exact-exchange operator with parameters \((\alpha,\beta,\mu)\).
When this spec lands, offline PySCF can pin molecular CAM-B3LYP energy and
the `cam-b3lyp` queue row can leave `blocked`.

## Domain basis

**Coulomb range separation** (Yukawa / error-function partition used by
libcint / PySCF / most RSH DFT):

\[
\frac{1}{r}
  = \frac{\mathrm{erfc}(\mu r)}{r}
  + \frac{\mathrm{erf}(\mu r)}{r}
  \equiv v^\mathrm{SR}(r;\mu) + v^\mathrm{LR}(r;\mu).
\]

**CAM-B3LYP** (Yanai, Tew, Handy, *Chem. Phys. Lett.* **393**, 51 (2004);
DOI [10.1016/j.cplett.2004.06.011](https://doi.org/10.1016/j.cplett.2004.06.011))
mixes HF exchange as

\[
E_\mathrm{x}^\mathrm{HF}
  = \alpha\, E_\mathrm{x}^\mathrm{HF}[1/r]
  + \beta\, E_\mathrm{x}^\mathrm{HF}[\mathrm{erf}(\mu r)/r]
\]

with the standard parameter set used by PySCF / Gaussian:

| Param | Value | Role |
| ----- | ----- | ---- |
| \(\alpha\) | 0.19 | short-range HF weight on bare Coulomb |
| \(\beta\) | 0.46 | long-range HF weight on erf-attenuated Coulomb |
| \(\mu\) | 0.33 \(a_0^{-1}\) | range-separation parameter |

Semi-local B88 exchange and VWN/LYP correlation pieces match the existing
B3LYP semi-local path; only the **HF-exchange operator** changes from a
single \(\alpha_\mathrm{HF}\cdot K[1/r]\) to the two-operator CAM form
above. Energy unit in SCF internals: Hartree; Frame public energy: eV
(`core/units.h`).

**Integral target (engine):**

\[
(\mu\nu|\lambda\sigma)_{\omega}
  = \int\!\!\int
    \chi_\mu(\mathbf{r}_1)\chi_\nu(\mathbf{r}_1)\,
    \frac{\mathrm{erf}(\omega r_{12})}{r_{12}}\,
    \chi_\lambda(\mathbf{r}_2)\chi_\sigma(\mathbf{r}_2)\,
    d\mathbf{r}_1\,d\mathbf{r}_2
\]

with \(\omega=\mu\) for CAM long-range. Short-range can be recovered as
\((1/r) - \mathrm{erf}/r\) or implemented with `erfc` if cheaper; v1 may
compute bare + erf and form SR by difference (matches many engines).

## Design

### Entities

| Piece | Location | Change |
| ----- | -------- | ------ |
| Rys ERI engine | `cpu/scf/rys_eri.{h,cpp}` | add attenuated block API |
| Boys / Rys roots | `cpu/scf/rys_roots.*`, `native_integrals` | extend Boys \(F_n(T)\) path for erf-Coulomb (or dedicated attenuated recurrence) |
| Hamiltonian builder | `cpu/scf/internal/hamiltonian_builder.*` | RSH exact-exchange admixture |
| XC catalog | `cpu/scf/xc_functional.{h,cpp}` | `XCType::CAMB3LYP` + param accessors |
| RKS product | `cpu/scf/rks.h` + stages | accept CAM; no public Config bag |
| Oracle / regression | `tools/pyscf_oracle/`, `regressions/scf/` | offline pin |

### Public API (binding shape)

```cpp
// --- Engine (native namespace; may stay internal to scf if only Fock uses it) ---
// Attenuated Coulomb (ab|cd)_ω with operator erf(ω r)/r, ω in 1/Bohr.
// ω = 0 must reproduce rys_eri_block within 1e-12 (identity gate).
void rys_eri_erf_block(const ShellBlock& A, const ShellBlock& B,
                       const ShellBlock& C, const ShellBlock& D,
                       double omega, double* out);

// --- Functional ---
enum class XCType { /* existing */, CAMB3LYP };

// CAM / RSH parameters (0 for non-RSH types)
struct RshExchangeParams {
    double alpha = 0;  // bare-Coulomb HF weight
    double beta  = 0;  // erf-Coulomb HF weight
    double mu    = 0;  // [1/Bohr]
};
RshExchangeParams rsh_exchange_params(XCType type);
// hybrid_exchange_fraction(CAMB3LYP) remains the *global* short form only
// if still used; prefer rsh_exchange_params for Fock.

// --- Product ---
cpu::scf::RKS rks(cpu::scf::XCType::CAMB3LYP);
rks.basis = cpu::scf::BasisSet::STO3G;  // or existing BasisSet member path
rks.tol = 1e-10;
rks.run(frame);  // frame.meta.get(kEnergy) [eV]
```

Rules:

- **No public `*Config` bags**, no `make_*` / `*_config` factories.
- Parameters are ctor args / members / `XCType` catalog (same pattern as B3LYP).
- Integrals and AO maps stay **private workspace** inside the method.
- v1 product: **closed-shell RKS energy only** (no UKS CAM, no analytic gradient,
  no TDDFT kernel for RSH).

### Reuse decision

| Candidate | Tag | Decision |
| --------- | --- | -------- |
| `native::rys_eri_block` / `rys_eri_block_deriv` | reuse / generalize | **generalize** — add parallel `rys_eri_erf_block` sharing shell-quartet dispatch, HRR, and root budget; do not fork a second ERI codebase |
| `hybrid_exchange_fraction` + `hamiltonian_builder` hybrid K path | generalize | **generalize** — when `rsh_exchange_params(type).beta != 0`, build \(K[1/r]\) and \(K[\mathrm{erf}/r]\) and mix \(\alpha,\beta\); pure global hybrids keep single-operator path |
| `evaluate_xc` B3LYP semi-local pieces | reuse | **reuse** — CAM semi-local is the same B88+VWN/LYP family as B3LYP with published scales; do not reimplement GGA kernels |
| `XCType::{B3LYP,PBE0}` catalog | pattern | **pattern** — add `CAMB3LYP` enum value + table rows like existing hybrids |
| PySCF offline pin path `tools/pyscf_oracle/` | reuse | **reuse** — same pin workflow as `b2plyp-h2o` / hybrid TDA regressions |

## Files to create or modify

- `include/atomiverse/cpu/scf/rys_eri.h` — declare `rys_eri_erf_block`
- `src/cpu/scf/rys_eri.cpp` — implement attenuated path
- `include/atomiverse/cpu/scf/rys_roots.h` / `src/cpu/scf/rys_roots.cpp` — if Boys/root path needs ω (or document reuse of existing roots with modified \(T\))
- `include/atomiverse/cpu/scf/xc_functional.h`
- `src/cpu/scf/xc_functional.cpp`
- `include/atomiverse/cpu/scf/internal/hamiltonian_builder.h`
- `src/cpu/scf/internal/hamiltonian_builder.cpp`
- `src/cpu/scf/stages/fock_stage.cpp` (if E_K accounting needs RSH split)
- `tests/cpu/scf/test_rys_eri_erf.cpp` (new)
- `tests/cpu/scf/test_xc_cam_b3lyp_params.cpp` (new)
- `regressions/scf/cam-b3lyp-h2-sto3g/` (new) + `REFERENCE.md`
- `tools/pyscf_oracle/` pin script or documented one-liner (extend existing)
- `docs/roadmap/dft.md` — flip CAM-B3LYP ○→✓ when closed
- `docs/developer/pyscf-method-queue.md` — `cam-b3lyp` blocked → done

## Tasks

- [ ] Write failing unit tests for `rys_eri_erf_block`: ω→0 identity vs `rys_eri_block`; known (ss\|ss) erf values for H-like STO-3G primitives (`tests/cpu/scf/test_rys_eri_erf.cpp`)
- [ ] Implement `rys_eri_erf_block` by generalizing the Rys engine (`rys_eri.h` / `rys_eri.cpp`, roots as needed)
- [ ] Write failing unit tests for `rsh_exchange_params(XCType::CAMB3LYP)` = (0.19, 0.46, 0.33) (`tests/cpu/scf/test_xc_cam_b3lyp_params.cpp`)
- [ ] Add `XCType::CAMB3LYP` + semi-local scales + `rsh_exchange_params` in `xc_functional.{h,cpp}`
- [ ] Generalize `hamiltonian_builder` (and Fock energy accounting) to mix bare + erf exchange for RSH types
- [ ] Offline PySCF pin (CAM-B3LYP / STO-3G cart H₂ or H₂O); write `REFERENCE.md` with version + command
- [ ] Add regression example `regressions/scf/cam-b3lyp-h2-sto3g/` (public API only; hard-coded golden; no live PySCF)
- [ ] Wire `atv_add_regression` under `ATV_BUILD_HF`; mark queue row + roadmap; Doxygen on new public symbols
- [ ] Run full check + test suite

## Testing strategy

**Unit (`tests/` only, single-function):**

1. `TestRysEriErf` — ω=0 matches bare ERI (max abs err ≤ 1e-12 on random s/p quartets).
2. `TestRysEriErf` — fixed primitive (ss\|ss) vs hard-coded reference (from offline PySCF/libcint or analytic limit).
3. `TestRshExchangeParams` — CAM constants; pure functionals return zeros; B3LYP stays non-RSH.

**Regression (`regressions/scf/cam-b3lyp-h2-sto3g/`):**

- Molecule: H₂ (or H₂O if H₂ is too soft), STO-3G Cartesian, CAM-B3LYP defaults.
- Observable: \(E_\mathrm{tot}\) [Ha] vs offline PySCF 2.14.x pin, tol **1e-6 Ha** (same class as B2PLYP / BLYP energy pins).
- Public path: `cpu::scf::RKS{XCType::CAMB3LYP}` + `run(frame)` only.

**Perf:** optional ratio vs PySCF; no standing ⚡ suite required; `perf_exempt` allowed if fair bench not wired.

## Out of scope

- Analytic nuclear gradients / Hessian for CAM-B3LYP
- UKS / ROKS CAM, TDA/TDDFT with RSH kernels
- Other RSH functionals (ωB97X, LC-ωPBE, HSE) — catalog can share the engine later
- CUDA / Kokkos ERI attenuation
- Optimal-tuning of \(\mu\)
- Public Config bags / factories
- Live PySCF at ctest time

## Notes for implementers

- Prefer **one** attenuated operator (`erf`) + bare Coulomb; form CAM mix in the Fock builder. Avoid three parallel full ERI engines.
- Keep `rys_eri_block` behavior bitwise-stable for existing hybrids (regression: B3LYP H₂ TDA / RPA cases must not move).
- Deferred: do **not** auto-run `/mol:impl` until an operator schedules this slug.
