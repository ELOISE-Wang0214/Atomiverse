# Acceptance — `gpw-kpoints-si`

System (regression): Si diamond, 4-atom bct cell
`cell = (a/√2, a/√2, a)` with a = 5.43 Å (10.26120 Bohr), atoms at
fractional (0,0,0), (½,½,½), (½,0,¼), (0,½,¾); DZVP-MOLOPT-GTH +
GTH-PBE-q4; XC = PBE; n_occ = 8 (16 valence e⁻).

All energies Hartree. Every criterion is verified by the single regression
binary `reg_gpw-kpoints-si` (hard-coded numbers, no external process) plus
the unit tests listed in spec.md.

Measured (first green 2026-07-25): A1 |Δ| = 1.4e-11 Ha; A2 |Δ| = 6.9e-14 Ha
(the A2 identity caught a real conjugation-side bug in the phased KB
nonlocal at first run — V_nl at −k, 2.09 Ha off); A3 |Δ| = 0.0;
A4 |Δ| = 0.599 Ha; A5 gap = 3.95 eV; A6 6 DIIS iters;
A7 E = −15.8023941333 Ha frozen; A9 wall ~17 s (OMP=8).

| # | Criterion | Type | Tolerance | Status |
|---|-----------|------|-----------|--------|
| A1 | **Γ-limit identity**: explicit `k_points = {(0,0,0)}, w=1` through the complex k-machinery reproduces the legacy real Γ path. Runs on a compact-basis periodic H₂ toy (DZV-GTH-PADE + GTH-PADE-q1, cubic 8 Bohr, PADE XC) — the legacy Γ path is impractical for diffuse-MOLOPT Si in a 7 Bohr cell (duplicated-triple collocation), which is itself part of what this slug fixes | scientific | \|ΔE\| ≤ 1e-6 Ha | **pass** (2026-07-25, reg first green) |
| A2 | **Supercell folding identity** (Si): E(bct, mesh 1×1×2 ⇒ κ_z = ±¼ TR→(0,0,¼), complex phases) equals E(z-doubled 8-atom cell at explicit κ = (0,0,½))/2, commensurate grids (unit-tested cutoffs) | scientific | \|E_super/2 − E_mesh\| ≤ 5e-5 Ha | **pass** (2026-07-25, reg first green) |
| A3 | **Time-reversal identity** (Si): MP 2×2×2 with `k_time_reversal=false` (8 k) equals `=true` (4 k, w doubled) | scientific | \|ΔE\| ≤ 5e-6 Ha | **pass** (2026-07-25, reg first green) |
| A4 | **k-sampling activation** (Si): \|E(MP 2×2×2) − E(k-path explicit Γ)\| > 1e-3 Ha on the 4-atom cell (k-mesh must matter; guards against the k path silently collapsing to Γ) | scientific | floor 1e-3 Ha | **pass** (2026-07-25, reg first green) |
| A5 | **Insulator band structure**: on the MP 2×2×2 run, max occupied ε < min virtual ε across all sampled k, gap > 0.1 eV; `eps_k`/`k_frac`/`k_weights` populated with matching shapes | scientific | gap floor 0.1 eV | **pass** (2026-07-25, reg first green) |
| A6 | **Convergence + finiteness**: MP 2×2×2 SCF converges (`converged == true`) with finite E within `max_iter` | functional | — | **pass** (2026-07-25, reg first green) |
| A7 | **Absolute e_tot band** (first-green frozen ATV value; replaced by a cited CP2K 2024.1+ number when an offline CP2K exists — twin input committed) | pin | \|ΔE\| ≤ 1e-3 Ha vs frozen value | **pass** (2026-07-25, reg first green) |
| A8 | Γ path regression-neutral: `k_mesh = {1,1,1}` and empty `k_points` routes the pre-existing real Γ code; scoped Γ gpw tests stay green | guard | existing tests pass | **pass** (2026-07-25, reg first green) |
| A9 | Regression ≤ 60 s wall (Release, OMP capped), no external process | perf/law | 60 s | **pass** (2026-07-25, reg first green) |

Performance: `within_10x` n/a — **perf_exempt** expected (no CP2K binary on
the worker host for a fair wall twin; verified `module avail cp2k` empty).
Record T_atv in REFERENCE.md.

Green requires: A1–A9 all pass in the registered regression + unit tests;
REFERENCE.md carries the pinned numbers, twin input echo, CP2K-absence
statement, and T_atv.
