# Method handbook pages

**Style (current generation): technical report.**  
Dense, sectioned, validation-first. **Not** the narrative Theory→Practice chapters
in `docs/tutorial/ch*.md`. A house narrative style will be applied in a later
rewrite pass.

Produced/updated by workflow **Tutorial** phases after each successful merge
wave:

- `pyscf-method-loop` — molecular post-HF / response (arg `write_tutorials`)
- `cp2k-method-loop` — GPW / PW / periodic DFT vs offline CP2K pins
- `metalwalls-method-loop` — classical electrochemical MD + LAMMPS-overlap MD
  (arg `write_tutorials`; LAMMPS oracle when feature overlaps; MW for MW-only)

## Required sections

| Section | Content |
| ------- | ------- |
| Summary | What / when / API / regression path |
| Method | Equations, assumptions, exclusions |
| Implementation | Minimal public C++; how to run `reg_<slug>` |
| Validation | Hard-coded pins from `REFERENCE.md` (PySCF, CP2K, or LAMMPS/MetalWalls) |
| Performance | Fair wall \(T_{\mathrm{ATV}}/T_{\mathrm{oracle}}\) (or `perf_exempt`) |
| Limitations | Scope + known debt |
| References | Paths + literature |

Template: [`_TEMPLATE.md`](_TEMPLATE.md).

## Naming

```text
docs/tutorial/methods/<slug>.md
```

`slug` = queue / regression directory name.

## Data sources (do not invent)

**PySCF campaign (molecular):**

- Accuracy: `regressions/scf/<slug>/REFERENCE.md`
- Speed: `docs/developer/pyscf-method-perf.md` or `tools/pyscf_oracle/fair_bench.py`
- Oracle rules: `docs/tutorial/appendix-pyscf-oracle.md` (PySCF never in ctest)

**CP2K capability (GPW / PW / periodic):**

- Accuracy: `regressions/{gpw,pw}/<slug>/REFERENCE.md`
- Speed / ratios: `docs/developer/cp2k-method-handoff.md` (fair wall or `perf_exempt`)
- Oracle: `tools/cp2k_oracle/README.md` (offline pins only; CP2K never in ctest)
- Queue / design: `docs/developer/cp2k-method-queue.md`, `cp2k-method-loop.md`

**MetalWalls / classical electrochem + LAMMPS-overlap MD:**

- Accuracy: `regressions/md/<slug>/REFERENCE.md` or
  `regressions/forcefield/<slug>/REFERENCE.md` (`oracle: lammps|metalwalls|literature`)
- Speed / ratios: `docs/developer/metalwalls-method-handoff.md` (fair wall or `perf_exempt`)
- Oracle: `tools/metalwalls_oracle/README.md` (offline pins only; never live in ctest)
- Queue / design: `docs/developer/metalwalls-method-queue.md`, `metalwalls-method-loop.md`
- **Oracle priority:** if the feature also exists in LAMMPS, pin and fair-bench
  vs **LAMMPS**. MetalWalls only for MW-only physics.

## Index — MetalWalls capability wave 1 (2026-07-23)

Classical electrochemical MD + LAMMPS-overlap MD (map E1–E3, E8–E9). Offline
pins only — see [metalwalls-method-handoff](../../developer/metalwalls-method-handoff.md)
and `tools/metalwalls_oracle/`. **Oracle priority:** LAMMPS when the feature
overlaps; MetalWalls / literature for MW-only physics.

| Method | Page | Regression | Oracle |
| ------ | ---- | ---------- | ------ |
| Const-pot BO parallel plate (madelung, E1) | [const-pot-bo-parallel-plate](const-pot-bo-parallel-plate.md) | `regressions/md/const-pot-bo-parallel-plate` | lammps ELECTRODE (`perf_exempt`) |
| Freeze electrode atoms (`setforce 0`) | [freeze-electrode-atoms](freeze-electrode-atoms.md) | `regressions/md/freeze-electrode-atoms` | lammps (ratio **3.98**) |
| Electrode total-charge observer | [const-pot-charge-observer](const-pot-charge-observer.md) | `regressions/md/const-pot-charge-observer` | lammps (`perf_exempt`) |
| Charging work \(-\sum_i \Psi_i Q_i\) | [charging-work-energy](charging-work-energy.md) | `regressions/md/charging-work-energy` | literature MW/Scalfi (`perf_exempt`) |
| Const-pot MD NVT short trajectory | [const-pot-md-nvt-short](const-pot-md-nvt-short.md) | `regressions/md/const-pot-md-nvt-short` | metalwalls (`perf_exempt`) |
| Gaussian electrode charge kernels (E2) | [gaussian-electrode-charge](gaussian-electrode-charge.md) | `regressions/md/gaussian-electrode-charge` | literature Reed/MW (`perf_exempt`) |
| Slab 2D Ewald foothold (Yeh–Berkowitz) | [slab-2d-ewald-foothold](slab-2d-ewald-foothold.md) | `regressions/md/slab-2d-ewald-foothold` | lammps (ratio **6.49**) |
| External uniform E-field MD | [external-efield-md](external-efield-md.md) | `regressions/md/external-efield-md` | lammps (ratio **0.60**) |
| Buckingham pair (LAMMPS `buck`) | [buckingham-pair-lammps](buckingham-pair-lammps.md) | `regressions/forcefield/buckingham-pair-lammps` | lammps (`perf_exempt`) |
| Four-site TIP4P foothold | [four-site-tip4p-foothold](four-site-tip4p-foothold.md) | `regressions/md/four-site-tip4p-foothold` | lammps (ratio **0.338**) |

## Index — MetalWalls capability wave 2 (2026-07-24)

Classical electrochemical MD + LAMMPS-overlap force fields (map E2–E4, E6–E8,
E10). Offline pins only — see
[metalwalls-method-handoff](../../developer/metalwalls-method-handoff.md) and
`tools/metalwalls_oracle/`. **Oracle priority:** LAMMPS when the feature
overlaps; MetalWalls / literature for MW-only physics. Fair walls from
REFERENCE / handoff (`ratio` or `perf_exempt`).

| Method | Page | Regression | Oracle |
| ------ | ---- | ---------- | ------ |
| Const-pot CG large multi-site electrode (E2) | [const-pot-cg-large-electrode](const-pot-cg-large-electrode.md) | `regressions/md/const-pot-cg-large-electrode` | metalwalls / dense BO dual (`perf_exempt`) |
| Electrode charge predictor (Kolafa ASPC, E2) | [electrode-charge-predictor](electrode-charge-predictor.md) | `regressions/md/electrode-charge-predictor` | metalwalls (`perf_exempt`) |
| Finite-field + single electrode (Dufils 3D, E3) | [finite-field-single-electrode](finite-field-single-electrode.md) | `regressions/md/finite-field-single-electrode` | metalwalls / literature (`perf_exempt`) |
| Bulk PIM self-consistent induced dipoles (E4) | [pim-induced-dipole-bulk](pim-induced-dipole-bulk.md) | `regressions/md/pim-induced-dipole-bulk` | metalwalls (`perf_exempt`) |
| Rigid piston electrode constant \(P_z\) (E6) | [piston-electrode-pz](piston-electrode-pz.md) | `regressions/md/piston-electrode-pz` | metalwalls (`perf_exempt`) |
| Thomas–Fermi finite metallicity electrode (E7) | [thomas-fermi-electrode](thomas-fermi-electrode.md) | `regressions/md/thomas-fermi-electrode` | metalwalls (`perf_exempt`) |
| Born–Mayer–Huggins pair foothold (E8) | [born-repulsion-foothold](born-repulsion-foothold.md) | `regressions/forcefield/born-repulsion-foothold` | lammps `pair_style born` (`perf_exempt`) |
| Fumi–Tosi multi-type Born + Coulomb NaCl (E8) | [fumi-tosi-nacl](fumi-tosi-nacl.md) | `regressions/forcefield/fumi-tosi-nacl` | lammps `hybrid/overlay born coul/cut` (`perf_exempt`) |
| Steele 10-4-3 wall (Ar single-point, E8) | [steele-wall-argon](steele-wall-argon.md) | `regressions/forcefield/steele-wall-argon` | lammps `wall/lj1043` (`perf_exempt`) |
| Capacitance from electrode charge fluctuations (E10) | [capacitance-charge-fluctuation](capacitance-charge-fluctuation.md) | `regressions/md/capacitance-charge-fluctuation` | literature (`perf_exempt`) |

## Index — MetalWalls capability wave 3 (2026-07-24)

Classical electrochemical MD + LAMMPS-overlap MD (map E3, E5, E9–E10 + shared
MD). Offline pins only — see
[metalwalls-method-handoff](../../developer/metalwalls-method-handoff.md) and
`tools/metalwalls_oracle/`. **Oracle priority:** LAMMPS when the feature
overlaps; MetalWalls / literature for MW-only physics. Fair walls from
REFERENCE / handoff (`ratio` or `perf_exempt`).

| Method | Page | Regression | Oracle |
| ------ | ---- | ---------- | ------ |
| Constant-force ion pull (`fix addforce`) | [constant-force-pull-ion](constant-force-pull-ion.md) | `regressions/md/constant-force-pull-ion` | lammps (ratio **0.12**) |
| Oscillating external E-field MD (E3) | [oscillating-efield-md](oscillating-efield-md.md) | `regressions/md/oscillating-efield-md` | lammps (ratio **2.12**) |
| Slab \(\Psi(z)\) planar Poisson observer (E10) | [potential-profile-observer](potential-profile-observer.md) | `regressions/md/potential-profile-observer` | literature (`perf_exempt`) |
| Surface-tension / NPγT coupling | [surface-tension-coupling](surface-tension-coupling.md) | `regressions/md/surface-tension-coupling` | lammps (ratio **2.88**) |
| Multi-topology virtual sites TIP4P+TIP5P (E9) | [virtual-site-multi-topology](virtual-site-multi-topology.md) | `regressions/md/virtual-site-multi-topology` | lammps (ratio **0.084**) |
| Coupled PIM dipoles ⊕ electrode \(Q\) (E5) | [dipole-electrode-joint](dipole-electrode-joint.md) | `regressions/md/dipole-electrode-joint` | metalwalls (`perf_exempt`) |

## Index — CP2K capability wave 1 (2026-07-23)

| Method | Page | Regression |
| ------ | ---- | ---------- |
| GPW Γ-RKS H₂O PBE (MT) | [gpw-water-pbe-mt](gpw-water-pbe-mt.md) | `regressions/gpw/gpw-water-pbe-mt` |
| GPW analytic forces H₂O | [gpw-forces-water](gpw-forces-water.md) | `regressions/gpw/gpw-forces-water` |
| Γ-UKS GPW O₂ | [gpw-uks-o2-gamma](gpw-uks-o2-gamma.md) | `regressions/gpw/gpw-uks-o2-gamma` |
| GPW PBE-D3(BJ) ice dimer | [gpw-pbe-d3bj-ice](gpw-pbe-d3bj-ice.md) | `regressions/gpw/gpw-pbe-d3bj-ice` |
| PW LDA Si (k-mesh) | [pw-lda-si](pw-lda-si.md) | `regressions/pw/pw-lda-si` |
| PW true PBE GGA Si | [pw-pbe-gga-si](pw-pbe-gga-si.md) | `regressions/pw/pw-pbe-gga-si` |
| Metallic Fermi–Dirac Al | [pw-smearing-al](pw-smearing-al.md) | `regressions/pw/pw-smearing-al` |
| PW collinear spin O₂ | [pw-spin-polarized-o2](pw-spin-polarized-o2.md) | `regressions/pw/pw-spin-polarized-o2` |
| PW BOMD H (NVE foothold) | [pw-bomd-h](pw-bomd-h.md) | `regressions/pw/pw-bomd-h` |
| OT SCF water (GPW PBE) | [ot-minimizer-water](ot-minimizer-water.md) | `regressions/gpw/ot-minimizer-water` |

## Index — CP2K capability wave 2 (2026-07-24)

Offline CP2K pins only — see [cp2k-method-handoff](../../developer/cp2k-method-handoff.md)
and `tools/cp2k_oracle/`. Fair walls from harness (`ratio` or `perf_exempt`).

| Method | Page | Regression | ratio \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}}\) |
| ------ | ---- | ---------- | --------------------------------------------- |
| GPW cell stress strain-FD H₂O (NUMERICAL pin) | [gpw-stress-water](gpw-stress-water.md) | `regressions/gpw/gpw-stress-water` | **4.75** |
| GPW BOMD H₂O NVE energy-drift foothold | [gpw-bomd-water](gpw-bomd-water.md) | `regressions/gpw/gpw-bomd-water` | **perf_exempt** |
| PW Nielsen–Martin analytic stress Si | [pw-analytic-stress-si](pw-analytic-stress-si.md) | `regressions/pw/pw-analytic-stress-si` | **perf_exempt** |
| FD supercell phonons diamond Si (Γ) | [fd-phonon-supercell](fd-phonon-supercell.md) | `regressions/pw/fd-phonon-supercell` | **perf_exempt** |
| Dudarev DFT+U NiO mol (Γ-UKS GPW) | [gpw-dftu-nio](gpw-dftu-nio.md) | `regressions/gpw/gpw-dftu-nio` | **0.21** |
| PW bands + DOS bulk Si | [bands-dos-si](bands-dos-si.md) | `regressions/pw/bands-dos-si` | **perf_exempt** |
| PW k-point collinear UKS O₂ | [pw-uks-kmesh-fe](pw-uks-kmesh-fe.md) | `regressions/pw/pw-uks-kmesh-fe` | **perf_exempt** |
| GPW Γ-RKS PBEsol Si₂ (MT) | [gpw-pbesol-si](gpw-pbesol-si.md) | `regressions/gpw/gpw-pbesol-si` | **perf_exempt** |
| GPW PBE-D4 periodic water dimer | [gpw-d4-water](gpw-d4-water.md) | `regressions/gpw/gpw-d4-water` | **perf_exempt** |

## Index — PySCF campaign methods

| Method | Page | Regression |
| ------ | ---- | ---------- |
| UMP2 | [ump2-o2-sto3g](ump2-o2-sto3g.md) | `regressions/scf/ump2-o2-sto3g` |
| Hybrid TDDFT-TDA (B3LYP) | [h2-b3lyp-hybrid-tda](h2-b3lyp-hybrid-tda.md) | `…/h2-b3lyp-hybrid-tda` |
| Frozen-core RMP2 | [frozen-core-mp2-h2o](frozen-core-mp2-h2o.md) | `…/frozen-core-mp2-h2o` |
| CIS(D) | [cisd-h2co-sto3g](cisd-h2co-sto3g.md) | `…/cisd-h2co-sto3g` |
| BLYP GGA | [blyp-water-sto3g](blyp-water-sto3g.md) | `…/blyp-water-sto3g` |
| GHF analytic gradient | [ghf-grad-h3-sto3g](ghf-grad-h3-sto3g.md) | `…/ghf-grad-h3-sto3g` |
| C-PCM free-energy forces | [pcm-grad-h2o](pcm-grad-h2o.md) | `…/pcm-grad-h2o` |
| SCS-MP2 | [scs-mp2-h2o-sto3g](scs-mp2-h2o-sto3g.md) | `…/scs-mp2-h2o-sto3g` |
| Full LR-TDDFT (RKS) | [full-lr-tddft-h2-pbe](full-lr-tddft-h2-pbe.md) | `…/full-lr-tddft-h2-pbe` |
| ROKS | [roks-o2-sto3g](roks-o2-sto3g.md) | `…/roks-o2-sto3g` |
| RHF static polarizability | [static-polariz-h2o](static-polariz-h2o.md) | `…/static-polariz-h2o` |
| IEF-PCM | [ief-pcm-h2o](ief-pcm-h2o.md) | `…/ief-pcm-h2o` |
| ROMP2 | [romp2-o2-sto3g](romp2-o2-sto3g.md) | `…/romp2-o2-sto3g` |
| ROKS analytic gradient | [roks-grad-o2-sto3g](roks-grad-o2-sto3g.md) | `…/roks-grad-o2-sto3g` |
| Hybrid full-LR TDDFT (B3LYP RPA) | [h2-b3lyp-rpa-tddft](h2-b3lyp-rpa-tddft.md) | `…/h2-b3lyp-rpa-tddft` |
| Spin-flip TDA (SF-CIS) | [spin-flip-tda-o2](spin-flip-tda-o2.md) | `…/spin-flip-tda-o2` |
| Spin-flip full TDDFT (SF-RPA) | [sf-tddft-o2](sf-tddft-o2.md) | `…/sf-tddft-o2` |
| ROHF nuclear Hessian | [rohf-hessian-o2](rohf-hessian-o2.md) | `…/rohf-hessian-o2` |
| ADC(2) | [adc2-h2co-sto3g](adc2-h2co-sto3g.md) | `…/adc2-h2co-sto3g` |
| ADC(2)-x | [adc2x-h2co-sto3g](adc2x-h2co-sto3g.md) | `…/adc2x-h2co-sto3g` |
| IP-ADC(2) | [ip-adc2-h2o-sto3g](ip-adc2-h2o-sto3g.md) | `…/ip-adc2-h2o-sto3g` |
| EA-ADC(2) | [ea-adc2-h2-sto3g](ea-adc2-h2-sto3g.md) | `…/ea-adc2-h2-sto3g` |
| RMP3 | [mp3-h2-sto3g](mp3-h2-sto3g.md) | `…/mp3-h2-sto3g` |
| RMP4(SDTQ) | [mp4-h2-sto3g](mp4-h2-sto3g.md) | `…/mp4-h2-sto3g` |
| CCSD | [ccsd-h2-sto3g](ccsd-h2-sto3g.md) | `…/ccsd-h2-sto3g` |
| CCSD nuclear gradient | [ccsd-grad-h2o-sto3g](ccsd-grad-h2o-sto3g.md) | `…/ccsd-grad-h2o-sto3g` |
| CC2 | [cc2-h2-sto3g](cc2-h2-sto3g.md) | `…/cc2-h2-sto3g` |
| CCSD(T) | [ccsd-t-h2o-sto3g](ccsd-t-h2o-sto3g.md) | `…/ccsd-t-h2o-sto3g` |
| G0W0 | [g0w0-h2-sto3g](g0w0-h2-sto3g.md) | `…/g0w0-h2-sto3g` |
| BSE@G0W0 | [bse-h2-sto3g](bse-h2-sto3g.md) | `…/bse-h2-sto3g` |
| B2PLYP double hybrid | [b2plyp-h2o](b2plyp-h2o.md) | `…/b2plyp-h2o` |
| CASCI | [casci-h2-sto3g](casci-h2-sto3g.md) | `…/casci-h2-sto3g` |
| CASSCF | [casscf-h2-sto3g](casscf-h2-sto3g.md) | `…/casscf-h2-sto3g` |
| SC-NEVPT2 | [nevpt2-h2-sto3g](nevpt2-h2-sto3g.md) | `…/nevpt2-h2-sto3g` |
| ΔSCF (IMOM) | [delta-scf-n2](delta-scf-n2.md) | `…/delta-scf-n2` |
| EOM-EE-CCSD | [eom-ee-ccsd](eom-ee-ccsd.md) | `…/eom-ee-ccsd` |
| EOM-SF-CCSD | [eom-sf-ccsd-o2](eom-sf-ccsd-o2.md) | `…/eom-sf-ccsd-o2` |
| EOM-CC2 | [eom-cc2-h2](eom-cc2-h2.md) | `…/eom-cc2-h2` |
| DFT-D4 (PBE params) | [dft-d4-h2o-pbe](dft-d4-h2o-pbe.md) | `…/dft-d4-h2o-pbe` |
| EOM-EA-CCSD | [eom-ea-ccsd](eom-ea-ccsd.md) | `…/eom-ea-ccsd` |
| EOM-IP-CCSD | [eom-ip-ccsd](eom-ip-ccsd.md) | `…/eom-ip-ccsd` |
| One-shot GF2 (G0F2) | [gf2-h2-sto3g](gf2-h2-sto3g.md) | `…/gf2-h2-sto3g` |
| Direct RPA | [direct-rpa-h2-sto3g](direct-rpa-h2-sto3g.md) | `…/direct-rpa-h2-sto3g` |
| RKS nuclear Hessian (FD of forces) | [rks-hessian-h2o-sto3g](rks-hessian-h2o-sto3g.md) | `…/rks-hessian-h2o-sto3g` |
| UKS nuclear Hessian (FD of forces) | [uks-hessian-o2-sto3g](uks-hessian-o2-sto3g.md) | `…/uks-hessian-o2-sto3g` |
| UKS full-LR TDDFT (LSDA O₂) | [uks-full-lr-o2-lsda](uks-full-lr-o2-lsda.md) | `…/uks-full-lr-o2-lsda` |
| UMP2 analytic gradient | [ump2-grad-o2-sto3g](ump2-grad-o2-sto3g.md) | `…/ump2-grad-o2-sto3g` |

Existing pages may still use an older Theory/Practice layout; next workflow run
or a bulk rewrite should move them to the technical-report template above.
