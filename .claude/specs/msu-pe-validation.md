---
title: MSU amorphous-PE tension — correctness and multi-core performance validation vs LAMMPS
status: in-progress
created: 2026-04-19
last-updated: 2026-04-19
parent: parallel-md-engine.md
---

# MSU amorphous-PE validation

## Status snapshot

| Phase | Sub-item | Status |
|---|---|---|
| 0.1 | Ghost-aware CellList | ✅ landed |
| 0.2 | AtomMigration (atom data via MPI_Alltoallv) | ✅ landed |
| 0.3 | Tag-based bonded topology + cross-rank bonds + LJ exclusions | ✅ landed |
| 0.4 | ReduceEnergy 8-double allreduce | ✅ landed |
| 0.5 | NoseHoover global-KE allreduce | ✅ landed |
| 0.6 | VirialStress + Berendsen Barostat (multi-rank-correct) | ✅ landed |
| 0.7 | Deform uniaxial affine box update | ✅ landed |
| A.2 | C20H42 single-chain force-field check | ✅ landed |
| A.3 | LJ special-bonds exclusions + multi-chain box | ✅ landed |
| A.4 | Single-chain NVE conservation 1000 steps | ✅ landed |
| A.5 | Single-chain NoseHoover NVT stability | ✅ landed |
| B.2 | 8-rank PE bonded NVE conservation | ✅ landed |
| B.3 | 8-rank NoseHoover NVT, real chain crossing rank boundary | ✅ landed |
| End-to-end | 8-rank pipeline (bonded + LJ + halo + NVT + NPT + tension) | ✅ landed |
| C | MSU 50k-atom system construction (data file) | ⏳ external -- needs LAMMPS data file or PSP/moltemplate generation |
| D | LAMMPS reference run + comparison script | ⏳ external -- needs LAMMPS install + scripted run |
| E.2/E.3 | Strong / weak scaling at 1, 2, 4, 8, 16, 32 ranks | ⏳ external -- needs cluster access |
| E.4 | vs LAMMPS head-to-head | ⏳ depends on D |

What can be done **right now** without external tooling:
* Any multi-rank PE simulation (NVE / NVT / NPT / tension) on systems
  the user constructs in code.
* Any qualitative check (energy conservation, no-NaN, rank-consistent
  state, box growth direction).

What needs **external setup**:
* Reading MSU's exact data file (a LAMMPS-format file -- requires the
  `molrs` LAMMPS data reader and a copy of the file).
* Quantitative match against LAMMPS (requires LAMMPS install + a
  scripted run + `tools/validate_polyethylene.py`).
* Multi-rank performance benchmarks at scale (single-machine 8 ranks
  works today; cluster 32-128 ranks needs MPI on a multi-node host).

---


Final acceptance deliverable for the parallel CPU MD engine: run the MSU amorphous-polyethylene uniaxial-tension protocol (<https://icme.hpc.msstate.edu/mediawiki/index.php/Deformation_of_Amorphous_Polyethylene.html>) on Atomiverse, cross-validate every observable against a LAMMPS reference run of the same system, and publish a multi-core performance report (strong + weak scaling, hybrid MPI×OpenMP, profile breakdown) that closes the parallel-md-engine roadmap.

**Primary success criterion**: at 8 MPI ranks on one node, Atomiverse produces a stress-strain curve σ_xx(ε) whose elastic modulus, yield stress, and hardening slope agree with LAMMPS within ±5% relative error, while wall-clock time per MD step is within 2× of LAMMPS at the same rank count.

**Secondary success criterion**: strong-scaling efficiency ≥ 70% at 32 ranks, weak-scaling efficiency ≥ 85% at 128 ranks (if a cluster is available), on the same system.

Everything in this spec is downstream of `parallel-md-engine.md`. Engine phases referenced as **"P2.x"** / **"P3.x"** refer to that spec.

---

## Dependency graph

```
Phase 0 (engine blockers)  →  Phase A (single-rank FF correctness)
                                    ↓
                         Phase B (multi-rank correctness)
                                    ↓
                         Phase C (MSU system setup)
                                    ↓
                  Phase D (LAMMPS correctness comparison)
                                    ↓
                Phase E (multi-core performance benchmarks)
                                    ↓
                      Phase F (reproducibility artifacts)
```

Phase 0 is hard-blocking. A through F are largely sequential but C.1 (system construction) can start in parallel with A.3+ once the bonded force field is correct on one rank.

---

## Phase 0 — Engine blockers (must land in parallel-md-engine before anything in this spec runs)

These are items from `parallel-md-engine.md` that are **not yet implemented** and are prerequisites for the MSU run specifically. Tracked here as gates, not redesigned.

| ID | Blocker | Spec ref | Notes |
|---|---|---|---|
| 0.1 | Ghost-aware `CellList` | P2.7 | Size bookkeeping to `dd->capacity`; bin `n_local + n_ghost`; bin grid covers halo skin. Unblocks multi-rank LJ NVE exit gate. |
| 0.2 | `AtomMigration` component | P3.1 | `MPI_Alltoallv` of atoms crossing `sublo/subhi`; bump `dd->neighbor_generation`; rebuild `tag_to_local`. PE chains will diffuse in 10⁵–10⁶ steps — without migration, owners drift from atoms. |
| 0.3 | Bonded cross-boundary | P3.2 | Ownership rule: bond = min-tag, angle = center, dihedral = 2nd atom. Bonded tuples migrate alongside atoms. Ghost skin widened to `cutoff + δ + bonded_extension`. |
| 0.4 | `ReduceEnergy` batched allreduce | P3.3 | 9-double `MPI_Allreduce` of `{energy, kinetic, virial[6], n_local}` at `end_of_step @ -50`. Required for meaningful thermo output. |
| 0.5 | `NoseHoover` multi-rank KE | P3.4 | Pull reduced KE (not rank-local) for thermostat coupling. Must run *after* `ReduceEnergy`. |
| 0.6 | `Barostat` multi-rank virial | P3.5 | Same: use globally-reduced virial. |
| 0.7 | `Deform` multi-rank affine | P3.6 | Every rank applies identical affine box update; ghost slab bounds recomputed; no desync of `dd->sublo/subhi` across ranks. |

**Exit criterion for Phase 0**: a 32-rank LJ + harmonic-bond run on 10⁵ atoms matches a 1-rank run to 1e-10 relative energy through 10 000 steps (= the existing Phase 3 exit gate of parallel-md-engine). When that lands, Phase A of this spec can proceed.

---

## Phase A — Single-rank force-field correctness

Goal: prove Atomiverse's CPU PE force field reproduces LAMMPS energies and forces on one rank, so any multi-rank discrepancies in Phase B are attributable to parallel bookkeeping, not to kernels.

### A.1 Force-field choice

Use **Dreiding-style united-atom PE** (same potential family the MSU reference uses):

- Bond: `HarmonicBond` — r_eq = 1.53 Å, k = 350 kcal/mol/Å² (convert to eV/Å²)
- Angle: `HarmonicAngle` — θ_eq = 109.47°, k = 60 kcal/mol/rad²
- Dihedral: `MultiHarmonicDihedral` — CH₂–CH₂–CH₂–CH₂ OPLS Fourier series (`K_1..K_4`)
- Non-bonded: `LJ` — σ = 4.01 Å, ε = 0.112 kcal/mol, r_cut = 10.5 Å, 1-4 scaling factor **0.0** (Dreiding convention; MSU uses exclusion, not scaling; revisit if needed)

All parameters stored in `benchmarks/md/reference/polyethylene/force_field.yaml` so LAMMPS and Atomiverse inputs both pull from one source of truth.

### A.2 Single-chain energy decomposition

New test `tests/cpu/md/test_pe_single_chain.cpp`:

- Build a single all-trans C₂₀H₄₂ UA chain (20 atoms, 19 bonds, 18 angles, 17 dihedrals).
- Evaluate `E_bond`, `E_angle`, `E_dihedral`, `E_vdW` with Atomiverse components at the same geometry.
- Compare each term against LAMMPS running the same geometry (reference written to `benchmarks/md/reference/polyethylene/c20h42_energies.json`).
- Tolerance: **1e-6 eV absolute**, per term.

### A.3 Small-box static energy + forces

New test `tests/cpu/md/test_pe_small_box.cpp`:

- 64 chains × 20 CH₂ = 1 280 atoms in a cubic box at ρ ≈ 0.85 g/cm³.
- Compare `E_total` and `max |Δf_i|` vs LAMMPS.
- Tolerance: **1e-5 relative** on total energy, **1e-6 eV/Å absolute** on any component of any atomic force.

### A.4 NVE energy conservation

Extend `tests/cpu/md/test_pe_small_box.cpp`:

- Run 1 000 steps NVE, `Δt = 0.5 fs`, starting from equilibrated geometry.
- Check `|E(t) − E(0)| / |E(0)| < 1e-4` throughout.

### A.5 NVT reference run

New benchmark `benchmarks/md/pe_nvt_small.cpp`:

- 1 000 steps NVT @ 500 K, `Δt = 1.0 fs`.
- Report `⟨T⟩`, `⟨E_pot⟩`, `⟨P⟩` block-averaged over last 500 steps.
- Compare against LAMMPS run of same duration and protocol.
- Tolerance: **0.5% relative** on each block average.

**Exit gate for Phase A**: A.2–A.5 all green on one rank.

---

## Phase B — Multi-rank correctness

Goal: prove that the multi-rank engine reproduces one-rank trajectories bit-for-bit within FP tolerance when the same seed, initial state, and potential are used.

All tests launched via `mpirun -n N`. Label `mpi` in CTest.

### B.1 Multi-rank LJ NVE (unlocks after P2.7)

Already stubbed in `tests/mpi/test_lj_nve.cpp` as `LjNve.EnergyConservationMultiRank`. Uncomment body. Tolerance: **1e-2 relative drift** over 20 steps at 8 ranks (loose — this is the engine-level smoke test, not a science-grade bound).

### B.2 Multi-rank bonded NVE

New test `tests/mpi/test_pe_bonded_nve.cpp`:

- 1 280-atom PE system (same as A.3).
- Run 100 steps NVE at 1 rank and at 8 ranks.
- Compare `frame.atoms.x[tag]`, `vx[tag]`, `fx[tag]` at step 100 after globally gathering to rank 0 and re-sorting by `tag`.
- Tolerance: **1e-10 absolute** on every field (FP-deterministic tolerance).

### B.3 Multi-rank NVT trajectories

New test `tests/mpi/test_pe_nvt_multirank.cpp`:

- 1 280-atom PE, 500 K, `Δt = 1.0 fs`, same seed.
- Run 2 000 steps at 1 rank and at 8 ranks.
- Block-average `T`, `E_pot`, `P_xx/yy/zz` over last 1 000 steps.
- Tolerance: **0.3% relative** on each block average (slightly tighter than A.5 since it's same-engine comparison).

### B.4 Multi-rank stress tensor

Inside `test_pe_nvt_multirank.cpp`:

- Log full virial tensor (`VirialStress`) every 100 steps.
- 1-rank vs 8-rank: each tensor component agrees to **1e-6 relative** at matching steps.

### B.5 Multi-rank RDF

Inside `test_pe_nvt_multirank.cpp`:

- Register `cpu::observer::RDF{r_max=10, bins=200, freq=100}` (async).
- Finalize after run.
- 1-rank and 8-rank RDF curves agree pointwise to **1e-8 absolute** on `g(r)`.

**Exit gate for Phase B**: B.1–B.5 all green on 1 and 8 ranks.

---

## Phase C — MSU amorphous-PE system construction

Goal: reproduce the MSU initial state so LAMMPS and Atomiverse start from the *same* atom coordinates, velocities, and topology.

### C.1 Initial configuration

Two candidate paths:

1. **Direct download** of the MSU LAMMPS data file if available from the MediaWiki link's supplementary material. Preferred — zero ambiguity.
2. **Regenerate** using `moltemplate` or the RadonPy / PSP packaging-tool chain with the MSU recipe:
   - 20 chains × 1 000 UA CH₂ (target — MSU uses different chain counts across examples; 20×1 000 = 20 000 atoms is a tractable size that matches parallel-md-engine Verification §7).
   - FCC seed at 1.53 Å spacing → random coil relaxation.
   - Box dimensions and ρ chosen to match MSU target.

Stored at `benchmarks/md/reference/polyethylene/msu_pe.data` (LAMMPS data-file format). Both LAMMPS and Atomiverse inputs read it (Atomiverse via `molrs` LAMMPS data reader).

### C.2 Equilibration protocol

Matches `parallel-md-engine.md` Verification §7:

| Stage | Ensemble | Duration | Δt | Purpose |
|---|---|---|---|---|
| 1 | NVT @ 500 K | 10 000 steps | 1.0 fs | initial randomization |
| 2 | NPT @ 500 K, 0 bar | 50 000 steps | 0.5 fs | density relaxation |
| 3 | NVT cooling 500 → 100 K (linear) | 50 000 steps | 0.5 fs | glass formation |
| 4 | NPT @ 100 K, 0 bar | 50 000 steps | 0.5 fs | glassy equilibration |
| 5 | Uniaxial tension along x | 200 000 steps | 0.5 fs | ε̇ = 10⁹ s⁻¹ to ε ≈ 1.0 |

Output: `Thermo`, `VirialStress`, `XYZTrajectory @ 10 000`, `RDF @ 1 000`.

### C.3 Atomiverse input

New file `benchmarks/md/polyethylene_tension.cpp` (already listed in parallel-md-engine §Files-to-create). Builds the Driver with `VelocityVerlet + NoseHoover + Barostat + Deform + VirialStress + LJ + HarmonicBond + HarmonicAngle + MultiHarmonicDihedral`, reads `msu_pe.data`, runs the 5-stage protocol above, writes thermo CSV + stress-strain CSV.

**Exit gate for Phase C**: single-rank Atomiverse completes all 5 stages on 20 000 atoms in ≤ 24 h wall time, writes correctly-formatted thermo and stress-strain CSVs, and its density after stage 4 is within 5% of the MSU reported amorphous-PE density at 100 K (~ 0.94 g/cm³).

---

## Phase D — LAMMPS reference comparison

Goal: demonstrate that Atomiverse and LAMMPS, starting from the same data file and running the same protocol, produce quantitatively matching observables.

### D.1 LAMMPS reference run

`benchmarks/md/reference/polyethylene/lammps.in` — LAMMPS input file implementing stages 1–5 of Phase C. Run at 1 and 8 ranks; collect:

- `log.lammps` → thermo CSV (`step, T, P, E_pot, E_kin, Pxx, Pyy, Pzz`)
- stress-strain CSV (`strain, sigma_xx`)
- trajectory XYZ every 10 000 steps
- final RDF

### D.2 Per-step thermo comparison (first 1 000 steps)

`tools/validate_polyethylene.py`:

- Reads LAMMPS `log.lammps` and Atomiverse `thermo.csv`.
- First 1 000 steps (deterministic regime before dissipative stochastics dominate): agreement to **1e-4 relative** on `T`, `E_pot`, `E_kin`; **1e-3 relative** on `Pxx/Pyy/Pzz`.
- Failure threshold: any step exceeding tolerance flagged with step number and max deviation.

### D.3 Equilibrium-state observables

Compare block averages over stage-4 (100 K NPT, last 25 000 steps):

| Observable | Source | Tolerance (relative) |
|---|---|---|
| ρ (density) | `Thermo` volume | 1% |
| ⟨E_pot⟩ | `Thermo` | 0.5% |
| RDF first-peak position r* | `RDF::finalize()` | 1 bin (< 0.05 Å) |
| RDF first-peak height g(r*) | `RDF::finalize()` | 3% |
| ⟨R²⟩ (end-to-end) | post-process `traj.xyz` | 10% (chain-statistics noise) |

### D.4 Stress-strain curve

`tools/validate_polyethylene.py`:

- Overlay Atomiverse σ_xx(ε) on LAMMPS σ_xx(ε).
- Extract:
  - elastic modulus E = dσ/dε at ε < 0.02
  - yield stress σ_y = max σ_xx for ε ∈ [0, 0.1]
  - strain-hardening slope dσ/dε at ε ∈ [0.5, 1.0]
- Tolerance: **±5% relative** on E, σ_y, hardening slope.

**Exit gate for Phase D**: D.2, D.3, D.4 all pass at 1 rank and at 8 ranks.

---

## Phase E — Multi-core performance benchmarks

**Primary target of this spec.** Goal: quantify Atomiverse parallel efficiency and compare to LAMMPS head-to-head.

### E.1 Single-rank baseline breakdown

`benchmarks/md/pe_profile_single.cpp` — runs MSU system stage 5 for 1 000 steps at 1 rank, with internal timers around:

- `HaloExchangeForward` (should be ~0 ns at 1 rank)
- `NeighborList::update + build_half`
- `LJ::launch`
- `HarmonicBond/Angle/MultiHarmonicDihedral::launch`
- `VelocityVerlet::launch`
- `HaloExchangeReverse`
- `ReduceEnergy`
- `SnapshotPublisher + RDF submission`

Report absolute and relative breakdown. Acceptance: no single component > 60% of step time (= no pathological hotspot); LJ + bonded combined ≥ 50% (forces dominating, as expected).

### E.2 Strong scaling (fixed problem size)

`benchmarks/md/scaling.cpp` at the MSU 20 000-atom system, 1 000 steps per rank-count:

| Ranks | Target efficiency | Notes |
|---|---|---|
| 1 | reference (1.00×) | — |
| 2 | ≥ 90% | — |
| 4 | ≥ 85% | — |
| 8 | ≥ 75% | matches MSU validation rank count |
| 16 | ≥ 60% | atoms-per-rank drops to 1 250 |
| 32 | ≥ 40% | comm-bound regime starts |

Report: time/step vs ranks (log-log plot), speedup, efficiency %. Output CSV for inclusion in `benchmarks/md/reference/polyethylene/scaling_report.csv`.

### E.3 Weak scaling (fixed atoms per rank)

Same benchmark, scaled system:

| Ranks | Atoms | Target efficiency |
|---|---|---|
| 1 | 20 000 | reference |
| 2 | 40 000 | ≥ 95% |
| 4 | 80 000 | ≥ 92% |
| 8 | 160 000 | ≥ 90% |
| 16 | 320 000 | ≥ 88% |
| 32 | 640 000 | ≥ 85% |
| 128 | 2 560 000 | ≥ 85% (if cluster) |

### E.4 Head-to-head vs LAMMPS

Same MSU 20 000-atom system, 1 000 steps, at matching rank counts (1, 2, 4, 8, 16, 32). LAMMPS command: `lmp -in lammps.in -var steps 1000 -pk ...`. Report:

| Ranks | Atomiverse time/step | LAMMPS time/step | Ratio (AV/LAMMPS) |
|---|---|---|---|
| 1 | — | — | target ≤ 2.0× |
| 8 | — | — | target ≤ 2.0× |
| 32 | — | — | target ≤ 2.5× |

Failure = ratio > stated target. Goal is *competitive*, not ≥ LAMMPS. Stretch goal: ratio ≤ 1.0× at any rank count.

### E.5 Profile decomposition vs ranks

At 1, 8, 32 ranks, repeat E.1 breakdown. Expected trends:

- `HaloExchangeForward + Reverse`: grows from 0% → 20–30% at 32 ranks.
- `ReduceEnergy`: grows from ~0% → 5–10%.
- Force components: drop proportionally as comm grows.
- Neighbor-list build: drops to near-zero (infrequent rebuilds).

Publish stacked-bar plot in `scaling_report.csv` companion.

### E.6 Hybrid MPI × OpenMP (optional, if time)

At 32 total cores, compare:

- 32 × 1 (pure MPI)
- 16 × 2
- 8 × 4
- 4 × 8
- 1 × 32 (pure OpenMP)

Report best configuration. Expected: MPI-heavy wins on non-NUMA single-socket; OpenMP-heavy may win on large-atom-per-rank regime.

### E.7 Multi-node (stretch, requires cluster access)

If a compute cluster is available: 128-rank run across 4 nodes × 32 ranks. Report weak-scaling efficiency at 128. Target: ≥ 70% efficiency vs single-node 32 ranks.

**Exit gate for Phase E**: E.1, E.2, E.3, E.4, E.5 all pass their target efficiencies on the local machine. E.6 and E.7 are stretch.

---

## Phase F — Reproducibility artifacts

Everything needed for an external reviewer to reproduce the results.

### F.1 Input files under version control

```
benchmarks/md/reference/polyethylene/
    force_field.yaml             # single source of truth for FF params
    msu_pe.data                  # LAMMPS data file (20 000 atoms)
    lammps.in                    # LAMMPS input for full protocol
    atomiverse.cpp               # symlink or note → benchmarks/md/polyethylene_tension.cpp
    c20h42_geometry.xyz          # Phase A.2 single-chain test input
    c20h42_energies.json         # Phase A.2 LAMMPS reference energies
    small_box_reference.json     # Phase A.3 LAMMPS reference (E, max|Δf|)
    thermo_reference.csv         # Phase D.2 reference thermo (1 rank and 8 ranks)
    stress_strain_reference.csv  # Phase D.4 reference stress-strain
    scaling_report.csv           # Phase E.2/E.3 scaling numbers
    scaling_report.md            # human-readable scaling writeup + plots
```

### F.2 Scripts

```
tools/
    run_msu_atomiverse.sh        # standardized launch at given rank count
    run_msu_lammps.sh            # same, for LAMMPS
    validate_polyethylene.py     # D.2 + D.3 + D.4 automated comparison
    plot_scaling.py              # E.2/E.3 plots from scaling_report.csv
    plot_profile.py              # E.5 stacked-bar plot
```

### F.3 README

`benchmarks/md/reference/polyethylene/README.md` containing:

1. One-paragraph summary and link to MSU source.
2. Exact commands to reproduce every phase (A–E) with expected outputs.
3. Expected-results table (modulus / yield / hardening / timings at 1/8/32 ranks) with ±tolerance.
4. "If this fails, check…" troubleshooting list.
5. Citation of the MSU reference and LAMMPS version used.

**Exit gate for Phase F**: a fresh clone + `./tools/run_msu_atomiverse.sh 8 && ./tools/validate_polyethylene.py` reproduces all green on a clean machine.

---

## Files to create

```
tests/cpu/md/test_pe_single_chain.cpp                # Phase A.2
tests/cpu/md/test_pe_small_box.cpp                   # Phase A.3, A.4
tests/mpi/test_pe_bonded_nve.cpp                     # Phase B.2
tests/mpi/test_pe_nvt_multirank.cpp                  # Phase B.3, B.4, B.5

benchmarks/md/pe_nvt_small.cpp                       # Phase A.5
benchmarks/md/polyethylene_tension.cpp               # Phase C.3 (listed in parallel-md-engine §Files)
benchmarks/md/pe_profile_single.cpp                  # Phase E.1
benchmarks/md/scaling.cpp                            # Phase E.2, E.3 (listed in parallel-md-engine §Files)

benchmarks/md/reference/polyethylene/force_field.yaml
benchmarks/md/reference/polyethylene/msu_pe.data
benchmarks/md/reference/polyethylene/lammps.in       # listed in parallel-md-engine §Files
benchmarks/md/reference/polyethylene/c20h42_geometry.xyz
benchmarks/md/reference/polyethylene/c20h42_energies.json
benchmarks/md/reference/polyethylene/small_box_reference.json
benchmarks/md/reference/polyethylene/thermo_reference.csv  # listed in parallel-md-engine §Files
benchmarks/md/reference/polyethylene/stress_strain_reference.csv  # listed in parallel-md-engine §Files
benchmarks/md/reference/polyethylene/scaling_report.csv
benchmarks/md/reference/polyethylene/scaling_report.md
benchmarks/md/reference/polyethylene/README.md

tools/validate_polyethylene.py                       # listed in parallel-md-engine §Files
tools/run_msu_atomiverse.sh
tools/run_msu_lammps.sh
tools/plot_scaling.py
tools/plot_profile.py
```

## Files to modify

```
src/cpu/molrec.cpp                    # confirm molrs LAMMPS data-file reader handles MSU data
CMakeLists.txt                        # register new benchmarks under ATV_BUILD_BENCHMARKS
.claude/specs/INDEX.md                # add this spec entry
```

---

## Acceptance checklist (at a glance)

- [ ] **0.1–0.7** engine blockers from parallel-md-engine.md all green
- [ ] **A.2** C₂₀H₄₂ energy decomposition vs LAMMPS, 1e-6 eV absolute
- [ ] **A.3** 1 280-atom static E + forces vs LAMMPS, 1e-5 rel E, 1e-6 force
- [ ] **A.4** 1 000-step NVE drift < 1e-4
- [ ] **A.5** NVT block averages vs LAMMPS, 0.5% rel
- [ ] **B.1** 8-rank LJ NVE drift < 1e-2 over 20 steps
- [ ] **B.2** 8-rank bonded NVE atoms/velocities match 1-rank at 1e-10
- [ ] **B.3** 8-rank NVT block averages match 1-rank at 0.3% rel
- [ ] **B.4** 8-rank stress tensor matches 1-rank at 1e-6 rel
- [ ] **B.5** 8-rank RDF matches 1-rank at 1e-8 absolute
- [ ] **C.1** MSU data file in repo; identical to LAMMPS input
- [ ] **C.2/C.3** full 5-stage protocol runs to completion, ρ within 5% of reference
- [ ] **D.2** per-step thermo vs LAMMPS first 1 000 steps, 1e-4 rel
- [ ] **D.3** equilibrium observables vs LAMMPS within listed tolerances
- [ ] **D.4** σ_xx(ε) modulus / yield / hardening within ±5% of LAMMPS
- [ ] **E.1** single-rank profile: no hotspot > 60%, force share ≥ 50%
- [ ] **E.2** strong-scaling efficiency targets met (1–32 ranks)
- [ ] **E.3** weak-scaling efficiency targets met
- [ ] **E.4** Atomiverse time/step ≤ 2× LAMMPS at 8 ranks
- [ ] **E.5** profile breakdown published at 1/8/32 ranks
- [ ] **F.1–F.3** all artifacts in `benchmarks/md/reference/polyethylene/` + reproducible from README

---

## Rationale

- **Why Dreiding-style UA, not exact MSU parameters?** MSU page does not publish a complete parameter file with every cross-term. Dreiding UA is the nearest widely-reproducible PE force field that Atomiverse already supports. Any science-grade parameter match is out of scope — this is an engine-validation spec, not a polymer-physics paper.
- **Why LAMMPS head-to-head instead of "match MSU numbers"?** LAMMPS is a reproducible, open reference; MSU's paper numbers depend on their specific FF parameter values and random seeds, which are not fully exposed. Using LAMMPS as the reference under *identical* inputs turns "correctness" into a deterministic FP comparison.
- **Why 2× LAMMPS as perf target, not ≥ LAMMPS?** LAMMPS has 20+ years of optimization. Matching it on first delivery is unrealistic; 2× is the industry norm for new engines at parity features. The scaling-efficiency targets (70–85%) are absolute performance-engineering criteria independent of the LAMMPS comparison.
- **Why Phase E gets the most detailed structure?** User's stated primary goal is multi-core performance. Correctness (A–D) is necessary groundwork; scaling + head-to-head comparison (E) is the deliverable.
- **Why hold off on MSU paper-number matching until later?** Two different unknowns (our FF parameters vs theirs, plus engine bugs) confounded in one test = impossible to debug. Phase D fixes the FF/protocol alignment via LAMMPS as intermediate reference; matching MSU paper numbers becomes a downstream Phase 5 addendum once both engines agree.
