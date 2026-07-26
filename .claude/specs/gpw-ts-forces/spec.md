# Spec: `gpw-ts-forces`

## Goal

Ship a **dedicated offline force pin** for free-atom **Tkatchenko–Scheffler
(TS)** ionic forces on GPW. Product energy (`gpw-ts-ice`) already reports
`e_disp` via `GpwDft::ts`; force assembly in `gpw_forces.cpp` already **folds**
geometry-only TS forces when `cfg.ts`. This slug closes the residual gap:
unit dual-gates + regression e2e + libMBD-class oracle notes (same policy as
`gpw-mbd-forces` / `gpw-d4-forces`).

## Scope

| In | Out |
| -- | --- |
| Free-atom TS forces \(F_A=-\partial E_{\mathrm{TS}}/\partial R_A\) [Ha/Bohr] | Density-chain Hirshfeld \(\partial(V/V_{\mathrm{free}})/\partial R\) |
| Isolated + periodic pure-RS (rcut 50 Å) dual-gates | Stock CP2K `TYPE TS` (CP2K 2026.1 has no TS pair type) |
| Wire verification: `assemble` / `gpw_forces_from_scf` with `ts=true` | Live oracle in ctest |
| Unit tests + `regressions/gpw/gpw-ts-forces` | MBD / XDM / D4 (sibling slugs) |
| Offline pin notes under `tools/cp2k_oracle/inputs/gpw-ts-forces/` | Molecular post-HF / classical MW / CUDA |

## API (no public Config bags / factories)

- Existing internal: `TsResult { energy, forces }`; `ts_pbe(...)` fills both
  (already on tree from `gpw-ts-ice`).
- Product: `GpwDft::ts = true` / `GPWEvaluator::ts` + `compute_forces` /
  `gpw_forces_*` — **no new public symbols**.
- Forces are the **negative gradient of the same free-atom TS energy** already
  matched to libMBD (`get_ts_energy` real-space path; PBE \(s_R=0.94\), \(d=20\)).

## Science

- Tkatchenko & Scheffler, *Phys. Rev. Lett.* **102**, 073005 (2009).
- libMBD: Hermann et al., *J. Chem. Phys.* **159**, 174802 (2023) — C API
  `cmbd_ts_energy` for isolated E_disp; pure-RS lattice sum for periodic.
- Free-atom \(\alpha_0\), \(C_6\), \(R_{\mathrm{vdW}}\): `basis_data/ts_tables.h`.
- Public Frame energy still eV; force pins Ha/Bohr (internal).

## Oracle pin (offline)

| Quantity | Source | Notes |
| -------- | ------ | ----- |
| E_disp isolated | libMBD `cmbd_ts_energy` (git 88069c6, ATV Å→Bohr) | −5.9764872295816794×10⁻⁴ Ha water dimer |
| E_disp periodic | pure-RS rcut 50 Å (same as `gpw-ts-ice`) | −6.099630773386013×10⁻⁴ Ha, 12 Å cell |
| F_TS | Independent central FD of ATV `ts_pbe` energy | Dual-gate max\|ΔF\| ≤ 1e-6 Ha/Bohr (h=1e-5 Bohr) |
| Wiring | \(F(\mathrm{ts})-F(\mathrm{plain})\) vs standalone `ts_pbe` forces | max\|Δ\| ≤ 1e-7 Ha/Bohr (geometry-only add) |

CP2K 2026.1 stock `PAIR_POTENTIAL TYPE` is only DFTD2/D3/D3(BJ)/D4 — **no TS**.
Oracle is **libMBD** + pure-RS dual formula (same as `gpw-ts-ice`). Twin notes:
`tools/cp2k_oracle/inputs/gpw-ts-forces/`.

## Acceptance

See `acceptance.md`. Summary:

1. Isolated E_disp still within **1e-12 Ha** of libMBD (unit).
2. Periodic E_disp within **1e-12 Ha** of pure-RS pin (unit) / **1e-8 Ha** (reg).
3. Product forces vs independent FD: max\|ΔF\| ≤ **1e-6** Ha/Bohr.
4. \(\|\sum_A F_A\| \le 1e-8\) Ha/Bohr.
5. Force assembly with `ts=true` adds TS forces (not energy-only).
6. Regression `reg_gpw-ts-forces` green ≤ 60 s.
7. No public `*Config` / `make_*` / `build_*`.
8. `perf_exempt` allowed (geometry-only kernel vs full SCF wall unfair).

## Files (expected)

- (existing) `include/.../ts_dispersion.h`, `src/.../ts_dispersion.cpp`
- (existing) `src/.../gpw_forces.cpp` TS branch — verify only
- `tests/cpu/gpw/test_cpu_gpw_ts.cpp` — periodic FD + product wiring
- `regressions/gpw/gpw-ts-forces/`
- `tools/cp2k_oracle/inputs/gpw-ts-forces/`
- `regressions/CMakeLists.txt` registration
