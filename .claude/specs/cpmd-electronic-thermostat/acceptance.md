# Acceptance — `cpmd-electronic-thermostat`

## Unit (`tests/cpu/gpw/test_gpw_cpmd.cpp` + thermostat cases)

| ID | Gate |
| -- | ---- |
| A1 | Existing Cpmd / RATTLE / GPWEvaluator guards still green |
| A2 | `set_thermostat` / members: non-positive \(K_0\) or \(\tau_e\) throws when thermostat on |
| A3 | `seed(C, Cdot)` with nonzero \(\dot C\): `kinetic_ha()` matches \(\tfrac12\mu\lVert\dot C\rVert_F^2\) |
| A4 | With thermostat + target \(K_0\): after several advances from hot seed, \(T_e\) moves toward \(K_0\) (final closer than initial, or within \(5\times K_0\)) |
| A5 | After thermostatted advance: RATTLE \(\lVert C^\top S C-I\rVert_\infty\le 1\mathrm{e}{-9}\), \(\mathrm{sym}(C^\top S\dot C)\le 1\mathrm{e}{-9}\) |
| A6 | `thermostat_energy_ha()` finite; thermostat off ⇒ energy 0 and \(\xi=\eta=0\) after rest seed |
| A7 | GPWEvaluator `build()` rejects `cpmd_thermostat` without `cpmd`, or bad \(K_0/\tau_e\) |

## Regression (`regressions/dynamics/cpmd-electronic-thermostat`)

| ID | Gate |
| -- | ---- |
| A10 | Public API: `MD` + `VelocityVerlet` + `GPWEvaluator{cpmd=true, cpmd_thermostat=true, …}` |
| A11 | Warm-up: finite PE \(<-100\) eV, nuclear KE ≥ 0 |
| A12 | Production: finite coords; O–H \(< 1.8\) Å |
| A13 | \(E_{\mathrm{ext}}=KE_{\mathrm{n}}+PE+T_{\mathrm{e}}+E_{\mathrm{therm}}\): relative drift \(\le 0.10\) |
| A14 | Production force steps: `last_result().iters==0` |
| A15 | \(T_e\) stays \(\le 0.05\) eV over production (does not run away) |
| A16 | Wall ideally ≤60 s; else `perf_exempt` in REFERENCE |
| A17 | REFERENCE cites NH electronic thermostat literature; offline twin under `tools/cp2k_oracle/inputs/` |

## Docs / queue

| ID | Gate |
| -- | ---- |
| A18 | On merge: roadmap electronic thermostat ✓; queue row done; spec deleted |

## Green criteria (worker schema)

- Dual gate: unit tests + `run_reg_cpmd-electronic-thermostat`
- `within_10x` or `perf_exempt` (literature envelope pin when no CP2K wall twin)
- `acceptance_met=true` only if A1–A17 hold
