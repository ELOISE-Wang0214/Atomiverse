# xl-bomd-linearized-kpoints — zero-SCF linearized XL-BOMD on multi-k RKS/UKS GPW

Status: **code-complete on `wf/xl-bomd-linearized-kpoints`, awaiting merge**
(2026-07-31; dual gate measured — delete this spec at close per house rule)
Queue: `docs/developer/cp2k-method-queue.md` (GapCensus wave 14 row; held
wave 1 behind the Γ-UKS pick, freed by its merge).

## Gap

`GPWEvaluator::build` throws `"xl_bomd_linearized is Γ-point only"`
(`src/cpu/gpw/gpw_evaluator.cpp:182`) — the Niklasson shadow-functional
zero-SCF force path is product on Γ RKS (`xl-bomd-linearized-forces`) and
Γ-UKS (`xl-bomd-linearized-uks`), while the k path only ships the
SCF-converged XL flavors (`xl-bomd-kpoints`, `xl-bomd-uks-kpoints`).

Underneath, the multi-k **UKS** solver has no fixed-density KS eval at all:
`run_gpw_scf_kpoints_uks` has no `max_iter == 0` branch and `GpwResult`
has no `H_ks_k_beta` field (the RKS k path got both for multi-k CPMD).

## Method (extends the shipped Γ product; Niklasson shadow functional)

Cold step (per-k XL histories not ready): unchanged — the existing SCF-XL
fall-through does one full multi-k SCF and seeds one `XlBomdCplx` per
sampled k (dual per spin for UKS).

Warm step, RKS (one complex aux density \(n_k\) per k):

1. Fixed-density eval (`max_iter = 0`, `guess_P_k =` \(\{n_k\}\)) →
   \(H_k[n]\), \(S_k\), \(E_{\mathrm{KS}}[n]\) (k-weighted, grid at
   \(n(\mathbf r)=\sum_k w_k n_k(\mathbf r)\)).
2. One generalized diagonalization per k (`internal::KOrth` +
   `internal::diagonalize_k`) → \(C_k\), \(\varepsilon_k\);
   \(D_k = 2\,C^{\mathrm{occ}}_k C^{\mathrm{occ}\dagger}_k\).
3. Shadow energy
   \(\mathcal E = \sum_k w_k\,\mathrm{Re}\,\mathrm{Tr}\!\big[H_k[n]\,(D_k-n_k)\big]
   + E_{\mathrm{KS}}[n]\) — recovers \(E_{\mathrm{KS}}\) at \(D_k=n_k\).
4. Forces: k-path assembly `gpw_forces_from_scf` on the payload
   \(\{P_k{=}D_k,\ C^{\mathrm{occ}}_k,\ \varepsilon_k\}\) (the exact-eigenvector
   case of the off-BO \(W = C(C^\dagger H C)C^\dagger\), matching the Γ
   product's `gpw_forces_at_density` semantics).
5. `advance(D_k)` per k; `iters = 0`, `converged = true` (stationary for the
   shadow functional). No iterative SCF anywhere in the warm step.

Warm step, UKS: same shape spin-diagonally — dual aux \(n^\sigma_k\), one
dual-spin fixed-density eval (**new** `max_iter == 0` branch in
`run_gpw_scf_kpoints_uks`, filling **new** `GpwResult::H_ks_k_beta`), one
diagonalization per (k, σ) with occ 1, shadow
\(\mathcal E=\sum_k w_k \sum_\sigma \mathrm{Re}\,\mathrm{Tr}[H^\sigma_k(D^\sigma_k-n^\sigma_k)]
+E_{\mathrm{KS}}[n^\alpha,n^\beta]\), forces from the dual k payload,
dual advance.

## Touch set

- `include/atomiverse/cpu/gpw/gpw_solver.h` — `GpwResult::H_ks_k_beta`;
  fixed-density payload doc covers multi-k UKS.
- `src/cpu/gpw/gpw_solver.cpp` — hoist `build_uks_k_energy` to function
  scope in `run_gpw_scf_kpoints_uks`; add the `max_iter == 0` dual-spin
  fixed-density branch (requires both spin guesses; refuses `use_ot` and
  metallic smearing — no-fallback).
- `src/cpu/gpw/gpw_evaluator.cpp` — drop the Γ-only throw; add the multi-k
  linearized warm branch (RKS + UKS); shared `density_from_cocc` /
  `flatten_occ` helpers (second use — CPMD multi-k branch already has
  private copies); **new guard**: `xl_bomd_linearized` + `compute_stress`
  throws at `build()` (the linearized warm branches return before the
  stress block — previously a silent drop on Γ too).
- `include/atomiverse/cpu/gpw/gpw_evaluator.h`,
  `include/atomiverse/cpu/electronic/xl_bomd_force.h` — doc updates.
- `tests/cpu/gpw/test_gpw_xl_bomd.cpp` — flip the two multi-k rejection
  pins to acceptance; add: k-weighted complex shadow identity, UKS-k
  fixed-density solver payload, linearized+stress guard.
- `regressions/dynamics/xl-bomd-linearized-kpoints/` — e2e (below).
- `tools/cp2k_oracle/inputs/xl-bomd-linearized-kpoints/` — CP2K twin input.
- Queue row, roadmap MD/DFT prose, `docs/tutorial/methods/` page.

## Oracle

CP2K 2026.1 **is** available on this host (`~/.local/bin/cp2k`, source rev
`85b8a9b49c`; also the EasyBuild module used by wave 1) — the family's
"no binary" `perf_exempt` basis is stale. CP2K has no zero-SCF
shadow-functional XL product, so the *shadow dynamics itself* has no
oracle twin; the **BO surface it must shadow does**: the same periodic-H₂
1×1×2 NVE twin that pinned `gpw-bomd-kpoints` (CP2K `.ener` per-step Pot
trace, CUTOFF 300, documented engine model band 4.07e-4 Ha). Re-run
offline for this slug, pin the trace in REFERENCE.md, and gate the
linearized shadow \(E_{\mathrm{pot}}(t)\) against it. Literature
(Niklasson PRL 100, 123004 (2008); JCP 130, 214109 (2009) Table I;
PRB 86, 174308 (2012); JCP 147, 054103 (2017)) still names the scheme;
the pin numbers are CP2K's.

## Regression (`regressions/dynamics/xl-bomd-linearized-kpoints`)

System: the `xl-bomd-kpoints` / `gpw-bomd-kpoints` twin — periodic H₂,
R = 1.0 Bohr along x, cubic 12 Bohr cell, DZV-GTH-PADE / GTH-PADE-q1,
PADE, CUTOFF 300 Ry, mesh 1×1×2 unreduced (2 k), NVE dt 0.5 fs,
v₀ = ∓0.001 Å/fs, 1 cold seed + 4 production steps.

Three legs in one binary:

- **A (anchor):** SCF-XL multi-k (`XlBomdForce{gpw, 5}`) — the shipped
  product, same trajectory settings.
- **B (product under test):** linearized multi-k
  (`XlBomdForce{gpw, 5, XlBomdLinearized{}}`), RKS.
- **C (UKS path):** linearized multi-k with `n_alpha = n_beta = 1` —
  spin-symmetric, so it must reproduce leg B's energies (exact
  degeneracy identity, dual-spin machinery genuinely exercised).

Gates:

1. B warm steps: `iters == 0`, `converged`, `n_kpts == 2`, `P_k` payload.
2. B per-step \(E_{\mathrm{pot}}\) vs the CP2K `.ener` Pot rows within the
   band (sibling band 2e-3 Ha + measured shadow/trajectory margin;
   exact band recorded in REFERENCE.md after measurement).
3. B NVE drift envelope \(\le 0.05\) (family policy).
4. A/B trace: per-step \(|E^{\mathrm{lin}}_{\mathrm{pot}}
   - E^{\mathrm{SCF-XL}}_{\mathrm{pot}}|\) within a measured band (the
   shadow residual is \(O((D-n)^2)\); band from measurement, expected
   ≪ 1 mHa class).
5. C vs B: per-step energy agreement within a tight spin-degeneracy band;
   C warm `iters == 0`.
6. Geometry sanity: finite coords, soft H–H window.

## Perf gate

Fair wall twin: CP2K 4-step MD (same input class as `gpw-bomd-kpoints`,
33.4 s measured wave 1) vs the ATV linearized 1+4-step leg. Expect
ratio ≪ 10 (sibling plain-BOMD leg was 0.064); measure and record —
`perf_exempt` only if the twin cannot be re-run this wave.

## Acceptance (all met 2026-07-31; numbers in the regression REFERENCE.md)

- [x] Multi-k RKS linearized: warm `iters == 0`; shadow \(E_{\mathrm{pot}}\)
      vs CP2K 2026.1 `.ener` worst **4.066e-4 Ha** (band 2e-3; ≈ the
      documented 4.07e-4 engine model band); A/B vs SCF-XL worst
      **1.099e-4 Ha** (band 1e-3); NVE drift 4.06e-3 (CP2K 4.4e-3).
- [x] Multi-k UKS linearized productized: leg C spin-degeneracy identity
      worst **2.31e-8 Ha** (band 1e-4) + `MultiKUksFixedDensityEvalFillsDualPayload`
      unit pin (E at converged dual densities within 1e-6 of the SCF).
- [x] Γ behavior unchanged: `test_cpu_gpw_xl_bomd` green (0.16 s);
      loop-domain suite at the merge gate (parent).
- [x] Unit: flipped acceptance pins, k-weighted shadow identity, UKS-k
      fixed-density payload + no-guess throw, linearized+stress guard.
- [x] Regression `reg_xl-bomd-linearized-kpoints` PASS, wall 7.7 s ≤ 60 s;
      CP2K trace re-run offline this wave (EasyBuild CP2K/2026.1 module).
- [x] Dual gate: ratio **1.3** (ATV leg B 1.96 s / CP2K 1.49 s, same host,
      OMP=8) — within_10x, `perf_exempt = false` (family's "no binary"
      basis is stale; binary present and used).
