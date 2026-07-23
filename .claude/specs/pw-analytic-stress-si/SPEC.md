# Spec: `pw-analytic-stress-si`

**Status:** open  
**Branch:** `wf/pw-analytic-stress-si`  
**Queue:** cp2k-method-loop · chapter `pw`  
**Dep:** `pw-lda-si` ✓ (strain-FD stress already ships as the oracle)

## Goal

Replace the multi-SCF **strain finite-difference** body of
`cpu::pw::pw_stress` with a true **Nielsen–Martin analytic** cell stress
(one SCF + Hellmann–Feynman cell derivatives). Keep strain-FD as
`pw_stress_fd` for dual-gate testing. Deliver an offline-pinned Si
regression (diamond, NC-UPF LDA) with public surface stress in **eV/Å³**.

## Cite

- O. H. Nielsen and R. M. Martin, *Phys. Rev. B* **32**, 3780 (1985)
  (analytic stress in plane-wave DFT).
- Quantum ESPRESSO stress decomposition (kinetic / Hartree / local /
  Ewald / XC / nonlocal KB) as algorithmic reference for NC
  pseudopotentials.
- PseudoDojo ONCVPSP-LDA Si (`tests/data/pw/Si_lda_nc.upf`); EOS
  parameters from `Si.djrepo` (same pin family as `pw-lda-si`).

## Public API (no `*Config` bags, no factories)

| Symbol | Role |
| ------ | ---- |
| `cpu::pw::pw_stress(frame, method)` | **Analytic** Voigt stress [eV/Å³] (xx,yy,zz,xy,xz,yz) |
| `cpu::pw::pw_stress_fd(frame, method, h=2e-3)` | Strain-FD oracle (same units/order); used only for dual-gate |
| `PWEvaluator::with_stress(true)` | Writes analytic stress to `frame.meta.get(kStress)` |

Parameters stay members on `cpu::pw::PW` (`ecut_ry`, `k_mesh`, `tol`,
`xc`, `upf_path`, …). Frame is the only I/O currency.

## Twin geometry (regression)

| Setting | Value |
| ------- | ----- |
| Structure | diamond Si, 2-atom fcc primitive |
| Lattice \(a\) | 5.43 Å |
| Basis | (0,0,0), \((a/4,a/4,a/4)\) |
| PP | `tests/data/pw/Si_lda_nc.upf` (ONCV LDA, \(Z_\mathrm{val}=4\)) |
| XC | `XCType::LSDA` |
| \(E_\mathrm{cut}\) | 15 Ry |
| k-mesh | \(2\times2\times2\) |
| SCF | `tol=1e-6`, `max_iter=50` |

Same incomplete-basis twin as `pw-lda-si` (≤ 60 s wall).

## Tasks

1. **Spec** (this file) + INDEX row.
2. **Public FD oracle:** expose `pw_stress_fd`; move current multi-SCF
   strain-FD body there; `pw_stress` → analytic.
3. **Analytic implementation** in `internal::pw_stress_analytic_impl`:
   after one SCF, accumulate
   kinetic (bands) + Hartree + local PP (form-factor \(G\)-derivative) +
   Ewald ion–ion + LDA XC volume term + KB nonlocal (when UPF projectors
   present). Units: internal Ha/Bohr³ → public eV/Å³.
4. **Unit tests** (`tests/cpu/planewave/test_pw_stress.cpp`):
   - cubic monoatomic symmetry (off-diag soft zero; diag equal);
   - analytic vs FD on cheap H (toy PP) within tol;
   - evaluator still writes `kStress` when `with_stress(true)`.
5. **Regression** `regressions/pw/pw-analytic-stress-si/`:
   - analytic stress finite, cubic-symmetry soft checks;
   - dual-gate: \(\|\sigma_\mathrm{an}-\sigma_\mathrm{FD}\|_\infty\) within
     stated tol on the Si twin (cheap FD step);
   - optional hydrostatic band from PseudoDojo EOS (soft, incomplete basis);
   - `REFERENCE.md` with offline pin notes and perf gate.
6. **Build + run** unit stress tests + regression until PASS.
7. **Perf gate:** fair wall vs CP2K if binary present; else
   `perf_exempt=true` (NC-UPF PW has no fair CP2K Quickstep twin for
   absolute stress — different PP).

## Acceptance criteria

| # | Criterion | Gate |
| - | --------- | ---- |
| A1 | `pw_stress` is **analytic** (one SCF class cost, not 12× SCF strain-FD) | code + wall vs FD |
| A2 | Public units eV/Å³ Voigt; `PWEvaluator::with_stress` writes `kStress` | unit test |
| A3 | Cubic H (toy): off-diag \(\lvert\sigma\rvert < 0.5\); diag mutually equal within 0.5 | unit + existing `h_stress_cubic` still PASS |
| A4 | Dual-gate H (toy, cheap ecut/k): \(\max_v\lvert\sigma^\mathrm{an}_v-\sigma^\mathrm{FD}_v\rvert \le 0.15\) eV/Å³ | unit test |
| A5 | Dual-gate Si twin: \(\max_v\lvert\sigma^\mathrm{an}_v-\sigma^\mathrm{FD,pin}_v\rvert \le 0.05\) eV/Å³ vs offline multi-SCF FD pin in REFERENCE.md (live an↔FD verified offline to \(7\times 10^{-5}\)) | regression |
| A6 | Si symmetry: \(\lvert\sigma_{xy,xz,yz}\rvert < 0.025\); \(\lvert\sigma_{xx}-\sigma_{yy}\rvert,\ldots < 0.025\) eV/Å³ | regression |
| A7 | Offline pin documented in `REFERENCE.md` (PseudoDojo / literature; CP2K STRESS_TENSOR note if no binary) | REFERENCE.md |
| A8 | `perf_exempt` **or** \(T_\mathrm{atv}/T_\mathrm{cp2k}\le 10\) recorded | REFERENCE.md |
| A9 | No public `*Config` / `make_*` / `build_*`; Frame I/O only | grep / review |

## Out of scope

- GPW stress (`gpw-stress-water`)
- GGA gradient stress refinements beyond LDA twin (PBE stress may follow)
- MPI multi-rank stress scaling
- Variable-cell product drivers beyond existing `IsotropicCellRelax` smoke

## Done when

All acceptance rows verified; commits on `wf/pw-analytic-stress-si` only;
`acceptance_met=true` in worker report.
