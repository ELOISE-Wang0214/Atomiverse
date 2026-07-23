# Specs index

## CPU breadth — next methods (2026-07-23)

CPU-only breadth-first queue from `docs/roadmap/index.md`. CUDA is a separate
ownership lane and is not part of this campaign.

- ~~fire-minimizer~~ — **SHIPPED 2026-07-23** Bitzek FIRE 1.0 as `cpu::minimizer::FIRE`; spec deleted, history in git log

## `backend-mirror` chain — CPU/CUDA layout rewrite (2026-07-22)

Discuss decisions: pillar=`frame` data only; no facade/shared/runtime/loop;
full-word dirs; per-backend implementations; Driver/Criterion; dissolve Backend;
**no backward compatibility; rewrite tests**. Implement in order:

- [backend-mirror-01-core](backend-mirror-01-core.md) — move contracts facade→core; delete facade [done]
- [backend-mirror-02-absorb](backend-mirror-02-absorb.md) — absorb top-level longrange/potentials; cuda kspace→longrange [done]
- [backend-mirror-03-paths](backend-mirror-03-paths.md) — full-word renames (dynamics, communication, …) [done]
- [backend-mirror-04-driver](backend-mirror-04-driver.md) — Driver/DriverNode rename; Criterion on CUDA [done]
- [backend-mirror-05-dissolve](backend-mirror-05-dissolve.md) — dissolve Backend into communication/device/asynchronous/driver [done]
- [backend-mirror-06-gates](backend-mirror-06-gates.md) — architecture scrubs + test layout freeze [done]
- [backend-mirror-07-docs](backend-mirror-07-docs.md) — CLAUDE/AGENTS/architecture.md narrative sync [done]

Start: `/mol:impl backend-mirror-01-core` (or `/mol:impl-all backend-mirror`).

**COMPLETE 2026-07-22** — all 7 phases shipped + verified (CPU 236/236); spec
files deleted, history in git log (`dev..backend-mirror`, commits ba03fac0…18b8a6ee).
Follow-ons split out: namespace→directory alignment; CUDA resurrection (below).

## `cuda-resurrect` — revive the dead CUDA backend (2026-07-22)

The CUDA backend does not compile on `dev`/`backend-mirror` (~180 pre-existing
errors: stale `Frame::atoms/box/scalars`, deleted `Real`/`core/real.h`, renamed
units constants) — orthogonal to backend-mirror, surfaced by its GPU validation.
Revive it purely by mirroring the current CPU Frame/precision/units API; **CPU
untouched (stays 236/236)**; GPU-node build+ctest is the acceptance.

- [cuda-resurrect](cuda-resurrect.md) — single spec, staged tasks (precision/units → frame boundary → domains → GPU-green) [done]

**COMPLETE 2026-07-23** — the CUDA backend configures, **builds and links** on
GH200 (`BUILD EXIT 0`) and passes **276/278** ctest; login-node CPU suite
**217/217 green** (zero CPU code touched). 14 GPU build rounds. Beyond the
Frame/precision/units mirror the operator also approved finishing the
**kind-base validation port** (6 CUDA kind tags + 15 components reparented + the
3 presets switched from name-substring to `dynamic_cast`, matching CPU's
contract and throw messages) — that retired the "CUDA kind bases pending" debt
and took ctest from 8 failures to 0 validation failures.
Gate: `tests/architecture/cuda_no_stale_api_scrub.sh` (ctest #212, login-node
runnable) locks the stale host-Frame / `Real` / units API out.

Two residual failures are **pre-existing and orthogonal** (documented, split out):

- [cuda-pme-accuracy](cuda-pme-accuracy.md) — CUDA PME reciprocal energy 3% off the
  NaCl Madelung reference (CPU hits 5e-4 on identical params). Ruled out: unit
  renames, float precision (fails under `-DATV_CUDA_REAL=double` too), missing
  terms, the Frame/scalars bridge. Needs on-device numerical bisect. [approved]
- `test_cpu_scf_ghf` — **CPU** UHF-pin test; PASSES on x86_64, fails only on
  aarch64 (FP-sensitive GHF/UHF basin). Not CUDA; routed to the SCF owner —
  see `.claude/notes/notes.md` debt entry.

---

Completed HF campaign work lives in the git log (`hf-capability-gaps-01`…`21`,
`scf-driver-refactor-06-criterion-unify` — closed 2026-07-21 by operator
manual attestation after review). Residual HF/Jensen/PySCF gaps are parked
in the inventory below; **do not** `/mol:impl` the umbrella — split a child
spec when work resumes.

- [hf-deferred-backlog](hf-deferred-backlog.md) — residual Jensen/PySCF HF gaps (forces, hard gates, response, solvent/x2c depth, basis, backends); parking lot only [approved]

## `hf-parity` chain — P0–P5 backlog decomposed (2026-07-22)

The whole `hf-deferred-backlog` inventory is expanded into the 21-link
`hf-parity` chain below, **all `[draft]`** (roadmap). Each is a real theme with
scope / design direction / files / provisional tasks / testing, but **no
acceptance ledger yet** — most need `/mol:litrev` + a pinned PySCF/literature
oracle before "done" is definable. To work a link: `/mol:spec <slug>` (adds the
researched design + acceptance, promotes to `approved`) → `/mol:impl <slug>`.
Priority = numeric order (P0 forces first; P5 debt last). Speculative links
(14–16, 20) promote only when a concrete consumer appears.

- ~~hf-parity-01-ecp-forces~~ — **SHIPPED 2026-07-22** (semi-numerical dU/dR; I₂ vs PySCF 1e-9 Ha/Bohr); spec deleted, history in git log
- ~~hf-parity-02-openshell-forces~~ — **ROHF forces SHIPPED 2026-07-22** (spin-Fock W; O2 vs PySCF 5.8e-11); GHF split out (no PySCF oracle); spec deleted, history in git log
- [hf-parity-03-environment-forces](hf-parity-03-environment-forces.md) — environment forces (umbrella): **embed + external-field forces SHIPPED 2026-07-22** (qmmm 5.4e-11; field FD 1.5e-7); RI-JK/C-PCM/x2c/CFMM-J members remain [draft, in progress]
- ~◐ hf-parity-04-ir-intensities — **APT (∂μ/∂R) tensor SHIPPED 2026-07-22** (H2O vs PySCF FD 4e-7 au); per-mode km/mol intensity projection deferred (no pyscf.prop oracle); spec deleted, history in git log
- [hf-parity-05-cfmm-j-gate](hf-parity-05-cfmm-j-gate.md) — CFMM-J hard gate (energy identity + crossover + honest ch19) [draft]
- [hf-parity-06-x2c-gate](hf-parity-06-x2c-gate.md) — x2c-1e hard gate (pVp vs PySCF, decoupling, Ar energy) [draft]
- [hf-parity-07-cpcm-gate](hf-parity-07-cpcm-gate.md) — C-PCM hard gate (ΔG_solv vs PySCF + honest docs) [draft]
- [hf-parity-08-embed-gate](hf-parity-08-embed-gate.md) — embed/field hard gate (finite-field dipole, MM bookkeeping) [draft]
- [hf-parity-09-docs-honesty](hf-parity-09-docs-honesty.md) — docs/Doxygen honesty sweep (Handler/SCFSolver/Results) [draft]
- [hf-parity-10-polarizability](hf-parity-10-polarizability.md) — static polarizability (E-field CPHF); optional β [draft]
- [hf-parity-11-uhf-ks-hessian](hf-parity-11-uhf-ks-hessian.md) — UHF CPHF/Hessian + RKS/UKS second derivatives (XC kernel) [draft]
- [hf-parity-12-tddft](hf-parity-12-tddft.md) — UHF-TDHF / spin-flip TDA; TDDFT (RKS/UKS A/B) [draft]
- [hf-parity-13-delta-scf-occupation](hf-parity-13-delta-scf-occupation.md) — MOM/IMOM/non-Aufbau ΔSCF + fractional occupation / Fermi smearing [draft]
- [hf-parity-14-advanced-solvent](hf-parity-14-advanced-solvent.md) — IEF-PCM/SS(V)PE/SMD, ddCOSMO, surface switching (speculative) [draft]
- [hf-parity-15-mm-embedding](hf-parity-15-mm-embedding.md) — polarizable/periodic MM embedding; Ewald-corrected charges (speculative) [draft]
- [hf-parity-16-relativity-depth](hf-parity-16-relativity-depth.md) — x2c picture-change props; DKH; SO-ECP; SOC/complex GHF (speculative) [draft]
- [hf-parity-17-population-analysis](hf-parity-17-population-analysis.md) — NPA/NBO, IAO/IBO, RESP [draft]
- [hf-parity-18-basis-expansion](hf-parity-18-basis-expansion.md) — static-table basis expansion + g angular momentum (no parsers) [draft]
- [hf-parity-19-ri-mp2](hf-parity-19-ri-mp2.md) — RI-MP2 product path / DF-ao2mo completeness [draft]
- [hf-parity-20-cuda-scf-gradients](hf-parity-20-cuda-scf-gradients.md) — CUDA SCF gradients (energy/density-only today; GPU node) [draft]
- [hf-parity-21-retire-scf-debt](hf-parity-21-retire-scf-debt.md) — retire SCFSolver/SCFConfig/IntegralPack/AtomList debt (do last; refactor) [draft]
