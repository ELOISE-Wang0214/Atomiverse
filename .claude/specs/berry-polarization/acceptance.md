# berry-polarization — acceptance

All criteria concrete: observable + tolerance vs pin. Gate = scoped build/test
per loop law (`atomiverse` + touched `test_*` + `reg_berry-polarization`).

## Criteria

1. **Unit — string geometry** (`test_cpu_pw_berry`):
   `BerryPolarization::string_kpoints(dir)` returns
   `k_perp² · k_string` fractional points, string-major, spacing `1/k_string`
   along `dir`, transverse offsets `i/k_perp`, all in `[0,1)`; throws on
   `dir ∉ [0,2]`, `k_perp < 1`, `k_string < 2`. Exact (1e-12).
2. **Unit — phase math**: `phase_from_link_dets` reproduces
   `arg ∏ det` for analytic complex links incl. modulus ≠ 1 and phase
   wrap-around (1e-12); throws on an (near-)zero det (|det| ≤ 1e-12).
3. **Unit — quantum wrap**: `wrap_centered(p, q)` maps into `(-q/2, q/2]`
   for positive/negative/boundary inputs (1e-12); throws on `q ≤ 0`.
4. **Regression `regressions/pw/berry-polarization`** (public API only,
   hard-coded pins, ≤ 60 s):
   a. SCF converged for both geometries; `p_ion_frac` of ideal Si =
      (1,1,1) within 1e-9 (Z_val=4 at (¼,¼,¼)).
   b. **Centrosymmetric pin** (Spaldin 2012): every computed total `p_i` of
      ideal diamond Si sits on the symmetry lattice — distance to the nearest
      multiple of 0.5 < **0.05** (units of e).
   c. **Born-charge pin** (ASR, Pick–Cohen–Martin 1970: Z*(elemental) = 0):
      sublattice displacement δ = 0.02 frac along a₃ →
      `|Z*| = |Δp₃|/δ < 0.5` — i.e. the electronic string phase cancels the
      bare ionic 4δ = 0.08 to better than 12.5 %. A sign error or a dead
      electronic term shows up as |Z*| ≈ 4 or 8.
   d. Consistency: displaced `Δp_ion,3 = 0.08` within 1e-9 (ionic bookkeeping).
5. **Architecture**: no public `*Config`/factory; class + members mirrors
   `Bands`/`Dos`; `run_restricted_scf` extracted (shared, not duplicated);
   `bands-dos-si` regression still green (same behaviour after extraction).
6. **Perf**: fair CP2K wall twin or `perf_exempt` with reason (no CP2K binary
   on build host — precedent `bands-dos-si`).

## Status (2026-07-25)

- [x] 1–3 unit tests green (`test_cpu_pw_berry`: strings, phase math, wrap,
      plus `jacobi_eigh_hermitian` residual/order/unitarity tests)
- [x] 4 regression green (`reg_berry-polarization`): ideal phases quantize to
      1e-6, p_ion=(1,1,1) exact, p on the half-quantum lattice, Z* = 0.0000
- [x] 5 architecture: `run_restricted_scf` extracted (shared by bands/DOS/
      Berry); no new public bags/factories; full PW suite green — 17
      `test_cpu_pw_*` unit tests + all 14 `regressions/pw/*` (two frozen
      offline method pins re-pinned with provenance after the solver
      root-cause fixes; their literature gates unchanged)
- [x] 6 perf: `perf_exempt=true` in REFERENCE.md (no CP2K binary on host);
      ATV twin wall ~18 s

Root-cause fixes shipped (iron law): `enumerate_g_vectors` index-bound
clipping; Davidson complex-Hermitian subspace (was real-truncated). CUDA
band solver has the analogous truncation — flagged for the CUDA backlog,
out of scope here.
