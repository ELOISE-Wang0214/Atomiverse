# MetalWalls method loop — handoff

Last land: **wave 5 — 2026-07-27** (10/10 staffed slugs green, all merged into
`dev` via `merge(wf/<slug>)` commits; harness at `dev` head `ca3b511a`).
All 10 wave-5 slugs → `## done` (done total **46**). Specs closed (10).
Pending residual (**3**) and blocked (**4**) remain for the next
GapCensus/assign. **incomplete_remaining=0**.

**Workflow order (binding):** Recover incomplete first → Census →
Assign (`incomplete ≫ pending`) → Worktrees (**spec-first**) → Merge → Harness →
Tutorial. See `docs/developer/metalwalls-method-loop.md`.

## Wave 5 land (2026-07-27, merged → `dev`)

| Slug | Component / change | Regression | oracle | ratio |
| ---- | ------------------ | ---------- | ------ | ----- |
| `per-electrode-eta` | `ConstPotBO::add_electrode(type, psi, eta)` per-group Gaussian widths, mixed \(\eta_{ij}\) | `reg_per-electrode-eta` | literature (Reed 2007 DOI 10.1063/1.2464084; Scalfi 2020 DOI 10.1063/5.0028232; MW JOSS 2020) | perf_exempt |
| `const-charge-electrode` | `ConstPotBO::add_electrode_charge` const-\(Q\) control (LAMMPS `electrode/conq` twin) | `reg_const-charge-electrode` | lammps | perf_exempt |
| `const-pot-charge-halo-mpi` | post-solve ghost-charge refresh under MPI DD (`ConstPot*` charge halo) | `reg_const-pot-charge-halo-mpi` | literature / self-consistency | perf_exempt |
| `maze-symplectic-charge-dynamics` | `ConstPotMaZe.symplectic` true SHAKE mass-zero path (Coretti) | `reg_maze-symplectic-charge-dynamics` | metalwalls (Coretti 2020 literature pin) | perf_exempt |
| `constant-d-polarization-feedback` | `FixDfield.polarization_feedback` self-consistent \(E = D - 4\pi k_e P\) (Zhang–Sprik) | `reg_constant-d-polarization-feedback` | metalwalls | perf_exempt |
| `extended-lagrangian-dipoles` | `cpu::PimXlDipole` time-reversible XL induced dipoles | `reg_extended-lagrangian-dipoles` | metalwalls (iEL/SCF literature pins — Albaugh/Niklasson; MW MaZe route cited) | perf_exempt |
| `qeq-fluctuating-charge` | `cpu::md::QEq` (LAMMPS `fix qeq/point` twin) | `reg_qeq-fluctuating-charge` | lammps | **0.45** |
| `pim-npt-polarizable-il` | bulk PIM + NPT polarizable-IL product twin | `reg_pim-npt-polarizable-il` | metalwalls (offline numpy Thole-PIM pin; MW JOSS 2020, Thole 1981) | perf_exempt |
| `voronoi-volume-tf` | `VoronoiVolume` TF site-volume companion | `reg_voronoi-volume-tf` | metalwalls | perf_exempt |
| `const-pot-il-edl-box` | planar const-pot + non-polarizable IL EDL production box | `reg_const-pot-il-edl-box` | metalwalls | perf_exempt |

All regressions registered behind `ATV_BUILD_REGRESSIONS`; every `REFERENCE.md`
states `oracle:` + `oracle_version`; ctest never launches LAMMPS/MetalWalls.

## Recover inventory (2026-07-27, post-harness wave 5)

| Source | Result |
| ------ | ------ |
| A queue `## in_progress` `[~]` | **0** (wave 5 fully merged) |
| B open `.claude/specs/` MD/electrochem | **0** (esmm-\*/dftb2-\*/periodic-rhf-gamma/rsh-eri-cam-b3lyp belong to other loops — never touch) |
| C unmerged `wf/<slug>` (this loop) | **0** (all wave-5 tips ancestors of `dev`) |
| D handoff incomplete rows | **none** |
| E half-written `regressions/{md,forcefield}/` | **0** MW slugs (method-ladder orphans unchanged, not this loop) |

**Resume set: empty.** **incomplete_remaining=0**.

## Next (after harness wave 5 — 2026-07-27)

Residual queue: **pending=3**, **blocked=4**, in_progress=0.

1. Pending (next assign staffs these): `nanoporous-carbon-const-pot` (E9
   product geometry), `maze-molten-salt-al` (MW example scenario; MaZe +
   Fumi–Tosi footholds landed), `thermo-integration-capacitance` (E10 TI
   workflow).
2. Blocked (out of loop, unchanged): `cuda-const-pot`, `aim-daim-product`,
   `ml-pinnwall-electrode`, `periodic-hybrid-electrode-dft`.
3. Next wave starts with **GapCensus** (map re-diff; foothold ≠ closed), then
   worktrees `wf/<slug>` off `dev`, **spec-first**, MD-only build.

Parent stays on **`dev`**. Default **staff_count=10**.

## Oracle

- Overlap with LAMMPS → **LAMMPS** offline pin (`tools/metalwalls_oracle/` + LAMMPS inputs)
- MW-only → **MetalWalls** offline pin (or literature when no binary)
- `ctest` never launches LAMMPS or MetalWalls
- Wave-5 oracle tags: `qeq-fluctuating-charge` + `const-charge-electrode` →
  **lammps**; `per-electrode-eta` + `const-pot-charge-halo-mpi` →
  **literature**; the other six → **metalwalls** (literature/numpy pins where
  no MW binary)

## Build (workers)

MD-only flags: `ATV_BUILD_MD=ON`, `ATV_BUILD_KSPACE=ON`, `ATV_BUILD_TESTS=ON`,
`ATV_BUILD_REGRESSIONS=ON`, **HF/PW/GPW/CUDA/Torch OFF**.

## History (compact)

- **Wave 1** (2026-07-23): 10 merged → done (const-pot BO/observer/NVT,
  charging work, Gaussian electrode, slab Ewald, E-field, Buckingham, TIP4P
  foothold, freeze). Ratios: freeze 3.98, slab 6.49, efield 0.60, tip4p 0.338;
  rest perf_exempt.
- **Wave 2** (2026-07-24): 10 merged → done (ConstPotCG, ASPC predictor,
  finite-field single electrode, bulk PIM, piston, Thomas–Fermi, Fumi–Tosi,
  Steele wall, Born, capacitance fluctuation). All perf_exempt.
- **Wave 3** (2026-07-24): 6 merged → done (force pull 0.12, oscillating
  efield 2.12, surface tension 2.88, virtual sites 0.084, \(\Psi(z)\) +
  dipole⊕electrode joint perf_exempt).
- **Wave 4** (2026-07-24 write, manual land re-verified 2026-07-25): 10 → done
  (ConstPotMaZe CG+ASPC foothold, electrode forces, global electroneutrality,
  Gaussian/Born NL bulk paths, `FixDfield` vacuum foothold, TIP4P EDL,
  \(C_{QQ}\), elastance observer, `FixMove` SMD). All perf_exempt.
- **Wave 5** (2026-07-27): table above; done total **46**.

Queue: `docs/developer/metalwalls-method-queue.md` (do not delete done rows).
Map: `.claude/notes/metalwalls-capability-map.md` (ATV marks re-audited
2026-07-27 post-wave-5).
