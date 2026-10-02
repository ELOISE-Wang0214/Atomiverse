<!-- mol:note:topic:es-architecture-debt -->
## [2026-09-24] ES public-shape debt

Science factorization is sound: shared SCF/energy loop + concrete
Hamiltonian / representation; PW / GPW / SE domain-split; method-as-nested
Driver for BOMD. **Engineering does not yet obey one public shape** across
method families. Review verdict **BLOCK** until the 🔴 rows below close
or are explicitly carved as named exceptions in CLAUDE.md.

**Rule**: Do not extend any dialect below (new Frame-bypass geometry,
new free-function force entry, new binding path, deeper inheritance).
Prefer closing toward: concrete method = `Driver<*::Phase>` (or one chosen
Component shape) + `run(Frame&)` / `forces(Frame&)` + one outer bind rule.
Detail lives here; CLAUDE.md keeps only a thin Remaining-debt pointer.

**Supersedes**: ad-hoc “Remaining debt” bullets that implied only CUDA RHF /
PW free-loop / GPW OT — those are rows in this ledger, not the whole story.

### Target principle (claimed)

| Axis | Claim |
| --- | --- |
| Ownership | Concrete method owns Hamiltonian + parameters |
| Loop | Self-consistency is a typed-stage `Driver` |
| I/O | Public currency is `Frame` (`run` / `forces`) |
| Kind | `HF` / `Semiempirical` / `PlaneWave` / `GPW` are energy-source tags only |
| Depth | Inheritance ≤ 2 |

### Live dialects (do not grow)

Gaussian post-HF, continuum solvent, the shared nuclear FD Hessian, and
CPHF response are **not** members of `MeanFieldMethod`. `Stability<RHF>`,
`TDA<RHF>`, `TDHF<RHF>`, `CISD<RHF>`, `G0W0<RHF>`, and `CCSD<RHF>::eom_ee_ccsd` are not members of
`RHF` (`pyscf.scf.stability.rhf_stability`, `pyscf.tdscf.rhf.TDA` /
`TDHF`, `gw.gw_exact.GWExact`, `cc.CCSD.eomee_ccsd_singlet`). CIS(D) is `CISD<RHF>` (PySCF has no CIS(D) class). Nuclear Hessians are
`Hessian<RHF|UHF|RKS|UKS|ROHF>` (`pyscf.hessian.rhf.Hessian` /
`hessian.uhf.Hessian`; KS and ROHF differentiate `forces`). KS linear
response is `TDA<RKS|UKS>` / `TDDFT<RKS|UKS>` (`pyscf.tdscf.rks` /
`uks`). Spin-flip is `SfTDA<UHF>` / `SfTDHF<UHF>`. `MP2` / `CCSD` /
`CASCI` / … take the mean-field reference (`pyscf.mp.mp2.MP2(mf)`,
`cc.CCSD.ccsd_t`). Solvent is `WithSolvent<Method>`
(`pyscf.solvent._attach_solvent.SCFWithSolvent`). The FD skeleton is
`central_difference_nuclear_hessian`. Nuclear CPHF is
`Hessian<RHF>` (`pyscf.hessian.rhf.Hessian(mf).kernel`); static
polarizability is `Polarizability<RHF>`. `pyscf.scf.cphf.solve` is not an
SCF method. Density fitting stays an SCF JK selector (`ri`);
`pyscf.df.df_jk.density_fit` returns another SCF.


1. **Method-as-Driver** — CPU Gaussian mean-field (`MeanFieldMethod`).
2. **Method + `*Evaluator`** — PW / AM1 / RM1 (Frame write often only on the
   wrapper; evaluator duplicates method params). GPW is the reference shape:
   Frame-in method + evaluator that only owns it.
3. **Loop inside Solver / Component** — NDDO SCC; some DFTB as bare `Semiempirical`.

### Debt ledger

Priority: 🔴 bind/forces/Driver alignment before any new ES public API;
🟡 are science/docs/harness hygiene.

| Sev | ID | Where | Debt |
| --- | --- | --- | --- |
| 🔴 | `es-cuda-gpw-frame-io` | `cuda/gpw/gpw_solver.cuh` | CUDA product still `build/run(vector<internal::GpwAtom>)` + public `cell[3]`; CPU twin is Frame-in (`build/run/forces(Frame&)`, `Kind` ctor arg). Mirror it ([cuda-gpw-public-shape.md](cuda-gpw-public-shape.md)); until then CUDA regs/tests reach the CPU oracle via `cpu::gpw::internal::GpwDftAccess`. |
| 🔴 | `es-bind-dialects` | `*Evaluator` vs direct nest vs DFTB-as-Component | GPW done: `GPWEvaluator(GpwDft, compute_forces)` owns the method, no parameter copy. PW / AM1 / RM1 evaluators still duplicate method params. Pick one outer bind story; document matrix; delete redundant shells. |
| 🔴 | `es-forces-api` | PW free functions; AM1 `compute_forces`; GPW stress | GPW has `forces(Frame&)` (eV/Å). PW forces are free functions, AM1 `compute_forces`; GPW stress is free-fn `gpw_stress(frame, …)` returning Ha/Bohr³ (not public units). Unify on Frame-writing members; keep FD oracles (`gpw_forces_fd`, `gpw_stress_fd`) as free functions. |
| 🔴 | `es-kind-model` | `cpu/component/hf.h` vs SE/PW/GPW Kinds | `HF` is pure tag; others are Component kinds. One model for energy sources. |
| 🔴 | `es-se-driver` | `am1_solver.h`, `rm1_solver.h`, NDDO SCC | Self-consistency not a typed Driver; RM1 lacks `QmEngine` / Kind / Frame write-back. |
| 🔴 | `es-aimd-depth` | `cpu/electronic/aimd_method.h` | `CarParrinelloMethod` inheritance depth ≥ 3. Compose, do not publicly derive. |
| 🔴 | `es-cuda-pw-driver` | `cuda/planewave/pw_solver.cuh` | No Driver; CPU/CUDA mirror broken (new vs CLAUDE). |
| 🟡 | `es-gpw-ot-freeloop` | `GpwDft` `use_ot` | OT still free-loop; DIIS paths are Driver stages — do not add a third topology. |
| 🟡 | `es-hf-cast` | `cpu/driver/driver.h` | `dynamic_cast<const HF*>` for energy sources vs protocol-probe law. |
| 🟡 | `es-ghost-atoms` | `MeanFieldMethod::ghost_atoms` | Selection on method, not Frame (named exception or move to Frame). |
| 🟡 | `es-dftb-qmengine` | DFTB0/3 / GFN0 vs DFTB2 / GFN1/2 | Uneven `QmEngine` coverage for QM/MM `set_engine`. |
| 🟡 | `es-kind-name-hf` | `hf.h` | Kind name `HF` covers pure KS; rename or document as Gaussian mean-field. |
| 🟡 | `es-docs-drift` | tutorials / `*Evaluator` comments | “mirrors Gaussian” without Evaluator caveat; ch09 `role()/capabilities()` prose (no Capability set exists); stale UKS / multi-k “throws” guard lists in `periodic-{mp2,ri-mp2,rpa,gw}-si.md`, `periodic-gw-kpoints.md`, `ot-smearing.md`; `second-gen-cpmd-uks.md` presents private `second_gen_cp*` flags as members. |
| 🟡 | `es-reg-internal` | `regressions/gpw/*` standalone dispersion pins | ~25 GPW cases still include `internal/{dftd3,dftd3_params,ts_dispersion,mbd_dispersion,xdm_dispersion}.h` to pin the standalone kernel — out-of-library consumer law. Pin through the public method (E_disp on `GpwResult`, F(disp)−F(plain)) or move the kernel pin to `tests/`. |
| 🟡 | `es-cuda-gpw-tests-e2e` | `tests/cuda/gpw/test_gpw_full_product.cu`, `test_gpw_kpoints.cu`, `test_gpw_stress.cu` | Run full SCFs (and CPU-oracle SCFs) under ctest — iron law 1. Move to `regressions/gpw/cuda-*`. |

### Conformance snapshot (at note date)

| Method family | Driver? | Frame `run`/`forces`? | Conforms? |
| --- | --- | --- | --- |
| CPU RHF/UHF/RKS/UKS/ROHF/ROKS/GHF | Yes | Yes | Pass (ghost_atoms / HF cast aside) |
| CUDA RHF | Yes | Yes | Pass |
| AM1 / RM1 | No | Partial / no | Fail (`es-se-driver`, `es-forces-api`) |
| CPU PW | Yes | `run` yes; forces free-fn | Partial |
| CUDA PW | No | `run` yes | Fail (`es-cuda-pw-driver`) |
| CPU GpwDft | Yes (OT residual) | Yes (`run`/`forces`; stress free-fn) | Pass (OT / stress units aside) |
| CUDA GpwDft | No | No (`internal::GpwAtom`) | Fail (`es-cuda-gpw-frame-io`) |
| `*Evaluator` wrappers | Component | launch → Frame | Pass as bind shell only |
| `CarParrinelloMethod` | via base | via base | Fail (`es-aimd-depth`) |

### Close order (recommended)

1. `es-forces-api` + `es-bind-dialects` (PW / AM1 / RM1 onto the GPW shape)
2. `es-se-driver` / `es-cuda-pw-driver` / `es-cuda-gpw-frame-io`
   (Driver + Frame shape on CUDA; CUDA RHF is already `cuda::Driver<scf::Phase>`)
3. `es-kind-model` + QmEngine coverage
4. Hygiene: `es-posthf-ks-gate`, docs, ghost_atoms, OT residual,
   `es-reg-internal`, `es-cuda-gpw-tests-e2e`
