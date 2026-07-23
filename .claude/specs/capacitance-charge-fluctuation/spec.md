# Spec: `capacitance-charge-fluctuation`

| Field | Value |
| ----- | ----- |
| Slug | `capacitance-charge-fluctuation` |
| Branch | `wf/capacitance-charge-fluctuation` |
| Domain | classical electrochemistry (E10 ecosystem / analysis) |
| Method | Capacitance / charge-fluctuation analysis |
| Chapter | E10 |
| Geometry | Parallel-plate EDL foothold (two electrode type groups) |
| Observable | \(C\) estimate [e/V] (differential from \(\langle(\delta Q)^2\rangle\); optional integral \(\langle Q\rangle/\Delta\Psi\)) |
| Public API | `cpu::CapacitanceChargeFluctuation` (Observer Component on `md::MD`) |
| Regression | `regressions/md/capacitance-charge-fluctuation/` |
| Oracle | **literature** (MetalWalls / Scalfi–Salanne charge-fluctuation formula; no MW binary on agent host) |
| Queue row | `docs/developer/metalwalls-method-queue.md` `## in_progress` `[~] capacitance-charge-fluctuation` |
| Map | MetalWalls capability map E10 — capacitance / charge-fluctuation tools |
| Deps | E1 ✓ (`ElectrodeChargeObserver` pattern); charge on Frame already written |

## Goal

Land a **Frame-immutable Observer** that accumulates electrode total charge
samples \(Q(t)\) (sum of `kAtom_charge` over owned atoms of a chosen electrode
type) and estimates capacitance:

\[
C_{\mathrm{diff}}
  = \frac{\langle Q^2\rangle - \langle Q\rangle^2}{k_{\mathrm{B}}T}
  = \frac{\langle(\delta Q)^2\rangle}{k_{\mathrm{B}}T}
  \qquad [\mathrm{e}/\mathrm{V}]
\]

(optional integral capacitance when a fixed \(\Delta\Psi\) [V] is supplied):

\[
C_{\mathrm{int}} = \frac{\langle Q\rangle}{\Delta\Psi}
  \qquad [\mathrm{e}/\mathrm{V}]
\]

with \(Q\) in elementary charge [e], \(T\) in K, and
\(k_{\mathrm{B}} = \texttt{units::kB\_ev\_per\_k}\) so that \(k_{\mathrm{B}}T\) is
in eV and \(1\,e\cdot 1\,\mathrm{V} = 1\,\mathrm{eV}\).

Does **not** solve electrode charges, does **not** run MD, and is **not** a
new Driver category — pure analysis on whatever \(q_i\) already sit on the
Frame (fixed charges, const-pot BO, etc.).

## Non-goals

- Constant-potential BO / CG / MaZe charge solve.
- Thomas–Fermi metallicity product (separate slug).
- Poisson potential-profile \(\Psi(z)\) observer (separate pending slug).
- Live MetalWalls or LAMMPS at ctest time.
- Welding capacitance into `KSpace*` or into a Driver preset.
- Public `*Config` bag, `make_*` / `build_*` factories.
- CUDA path.

## Architecture (binding)

- **Kind:** `cpu::Observer` (`FrameImmutable | Deferred`).
- **Bind:** default `md::Phase::end_of_step` (after force / charge update),
  order after a typical charge source / `ElectrodeChargeObserver` (e.g. +50).
- **I/O:** read `kAtom_charge` + `kAtom_type` on owned atoms only
  (`owned_count`); never mutates x/v/f/box/topology/charges.
- **Primary electrode:** constructor electrode type list; **sample uses
  electrode index 0** for \(Q\) (parallel-plate left wall). Remaining types
  are accepted for API symmetry with `ElectrodeChargeObserver` and may be
  inspected via `total_charge(e)` from the last sample only (not used in \(C\)).
- **Running moments:** Welford / sum-\(Q\), sum-\(Q^2\), sample count; no
  trajectory storage.
- **Publish:** members `mean_charge()`, `charge_variance()`,
  `differential_capacitance()`, `integral_capacitance()`, `n_samples()`;
  Frame meta `kCapacitanceDiff` / `kCapacitanceIntegral` [e/V] when
  \(n_{\mathrm{samples}}\ge 1\) (integral only if \(\Delta\Psi \ne 0\)).
- **MPI:** local sum of \(Q\) on owned atoms, then `MPI_Allreduce(SUM)` when a
  domain communicator is present (same as electrode charge observer); moments
  accumulate the global \(Q\) once per step on every rank.

## Oracle decision

| Check | Result |
| ----- | ------ |
| Same algorithm in LAMMPS? | **No** — core LAMMPS has no built-in electrode charge-fluctuation \(C\) compute (ELECTRODE package not installed on oracle `lmp` either). |
| MetalWalls-only analysis? | **Yes** — MW post-processing / JOSS ecosystem tools estimate \(C\) from \(Q(t)\). |
| MW binary on agent host? | **No** |

Therefore offline pin is the **closed-form literature formula** (Scalfi /
Limmer / MetalWalls JOSS lineage) evaluated on a **synthetic** \(Q\) sample
sequence that the regression drives by rewriting wall charges between
observer launches. `oracle: literature` with full citation;
`oracle_version` names the formula source. `perf_exempt` (no MW/LAMMPS twin
binary for this analysis micro-kernel).

## Tasks

1. Write this spec + `acceptance.md` on `wf/capacitance-charge-fluctuation`.
2. Offline pin: literature formula on a documented \(Q\) sample series + \(T\),
   \(\Delta\Psi\); record in
   `regressions/md/capacitance-charge-fluctuation/REFERENCE.md` with
   `oracle: literature` and `oracle_version`.
3. Implement `cpu::CapacitanceChargeFluctuation` under `cpu/observe/`; add meta
   keys; wire CMake (MD sources) + unit test under `tests/cpu/observe/`.
4. E2E regression `reg_capacitance-charge-fluctuation` (public headers only):
   parallel-plate Frame; drive known \(Q\) samples; assert \(C_{\mathrm{diff}}\)
   (and \(C_{\mathrm{int}}\)) vs pin within tol.
5. MD-only configure/build (HF/PW/GPW/CUDA/Torch **OFF**); run unit test +
   regression until PASS.
6. Perf gate: no twin binary → `perf_exempt` with reason.
7. Commit on `wf/capacitance-charge-fluctuation` only (never `dev`, never push).

## Acceptance criteria

See `acceptance.md` (binding).
