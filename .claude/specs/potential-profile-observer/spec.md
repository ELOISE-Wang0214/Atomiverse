# Spec: `potential-profile-observer`

| Field | Value |
| ----- | ----- |
| Slug | `potential-profile-observer` |
| Branch | `wf/potential-profile-observer` |
| Domain | classical electrochemistry (E10 ecosystem / analysis) |
| Method | Slab electrostatic potential profile \(\Psi(z)\) (planar Poisson) |
| Chapter | E10 |
| Geometry | Parallel-plate capacitor / EDL slab foothold |
| Observable | \(\Psi(z)\) curve [V], \(\rho(z)\) [e/Å³], \(E_z(z)\) [V/Å] |
| Public API | `cpu::PotentialProfile` (Observer Component on `md::MD`) |
| Regression | `regressions/md/potential-profile-observer/` |
| Oracle | **literature** Poisson planar (MetalWalls potential-profile analysis lineage; no MW binary on agent host) |
| Queue row | `docs/developer/metalwalls-method-queue.md` `## in_progress` `[~] potential-profile-observer` |
| Map | MetalWalls capability map E10 — potential-profile post-processing |
| Deps | E1 charges on Frame; long-range not required for this analysis observer |

## Goal

Land a **Frame-immutable Observer** that bins atomic charge into a 1-D
slab profile \(\rho(z)\) and solves the **planar Poisson equation** in MD
metal units to produce the electrostatic potential profile \(\Psi(z)\):

\[
\frac{\mathrm{d} E_z}{\mathrm{d} z}
  = 4\pi\,k_{\mathrm{e}}\,\rho(z)
  = \frac{\rho(z)}{\varepsilon_0},
\qquad
E_z = -\frac{\mathrm{d}\Psi}{\mathrm{d} z},
\]

with \(k_{\mathrm{e}} = \texttt{units::coulomb}\) (14.399645 eV·Å/e²),
\(\rho\) in e/Å³, \(E_z\) in V/Å, \(\Psi\) in V (eV/e). Gauge:
\(\Psi(z_{\min})=0\) at the left edge of the binning window; left-edge
field \(E_z(z_{\min})=0\) (ideal neutralized capacitor / vacuum left of
the left wall).

Does **not** solve constant-potential charges, does **not** run MD, and
is **not** a new Driver category — pure analysis on whatever \(q_i\)
already sit on the Frame.

## Non-goals

- Constant-potential BO / CG / MaZe charge solve.
- Capacitance fluctuation (separate done slug).
- 3-D / reciprocal-space potential from PME (this is the **planar slab
  post-process** MetalWalls tools provide).
- Live MetalWalls or LAMMPS at ctest time.
- Public `*Config` bag, `make_*` / `build_*` factories.
- CUDA path.
- New Driver category; charge solve inside `KSpace*`.

## Architecture (binding)

- **Kind:** `cpu::Observer` (`FrameImmutable | Deferred`).
- **Bind:** default `md::Phase::end_of_step` (after force / charge update),
  order +45 (near other electrode observers).
- **I/O:** read `kAtom_charge` + `kAtom_z` (and box for \(A=L_x L_y\)) on
  **owned** atoms only (`owned_count`); never mutates x/v/f/box/topology/charges.
- **Ctor args (no Config bag):** `n_bins` (\(>0\)), `z_min` [Å], `z_max` [Å]
  with \(z_{\max}>z_{\min}\). Axis is **z** only in this foothold.
- **Algorithm (each launch):**
  1. Zero bin charge sums; for each owned atom with
     \(z\in[z_{\min},z_{\max})\), add \(q_i\) to the matching bin.
  2. MPI_Allreduce(SUM) of bin charges when a domain communicator is present.
  3. \(\rho_i = Q_i / (A\,\Delta z)\) with \(A = L_x L_y\),
     \(\Delta z = (z_{\max}-z_{\min})/n_{\mathrm{bins}}\).
  4. Piecewise-constant \(\rho\): edge field
     \(E_{i+1}=E_i + 4\pi k_{\mathrm{e}}\rho_i\Delta z\), \(E_0=0\);
     edge potential
     \(\Psi_{i+1}=\Psi_i - \bigl(E_i\Delta z + \tfrac12 4\pi k_{\mathrm{e}}\rho_i\Delta z^2\bigr)\),
     \(\Psi_0=0\); publish **bin-center** \(E\) and \(\Psi\).
- **Publish:** members `n_bins()`, `z_min()`, `z_max()`, `dz()`,
  `rho()`, `field()`, `psi()`, `z_center(i)`; optional Frame meta
  `kPotentialDrop` = \(\Psi(z_{\max})-\Psi(z_{\min})\) [V] (right-edge
  minus left-edge = \(\Psi_{n}\)).
- **MPI:** bin locally, then allreduce bin charge sums; every rank
  integrates identically.

## Oracle decision

| Check | Result |
| ----- | ------ |
| Same algorithm in LAMMPS? | **No** built-in planar Poisson \(\Psi(z)\) compute for electrode slabs (chunk/ave + external post only). |
| MetalWalls-only analysis? | **Yes** — MW ecosystem potential-profile tools. |
| MW binary on agent host? | **No** |

Therefore offline pin is the **closed-form planar Poisson / parallel-plate
capacitor** (literature) evaluated on a synthetic two-sheet charge
geometry, with discrete bin values recorded from the same integration
scheme. `oracle: literature`; `perf_exempt` (no MW/LAMMPS twin binary for
this analysis micro-kernel).

## Tasks

1. Write this spec + `acceptance.md` on `wf/potential-profile-observer`.
2. Offline pin: parallel-plate sheets \(Q=\pm 1\,e\), \(A=100\) Å², bin
   window; record \(E_{\mathrm{gap}}\), \(\Psi\) key bins, continuum
   cross-check in `regressions/md/potential-profile-observer/REFERENCE.md`
   with `oracle: literature` and `oracle_version`.
3. Implement `cpu::PotentialProfile` under `cpu/observe/`; add meta key;
   wire CMake (MD sources) + unit test under `tests/cpu/observe/`.
4. E2E regression `reg_potential-profile-observer` (public headers only):
   assert mid-gap \(E\), potential drop, key \(\Psi\) bins vs pin.
5. MD-only configure/build; run unit test + regression until PASS.
6. Perf gate: no twin binary → `perf_exempt` with reason.
7. Commit on `wf/potential-profile-observer` only (never `dev`, never push).

## Acceptance criteria

See `acceptance.md` (binding).
