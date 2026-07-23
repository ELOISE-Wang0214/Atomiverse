# Spec: `oscillating-efield-md`

| Field | Value |
| ----- | ----- |
| Slug | `oscillating-efield-md` |
| Branch | `wf/oscillating-efield-md` |
| Domain | classical MD / electrochemistry (map **E3** residual) |
| Method | Oscillating / AC uniform external electric field on MD charges |
| Observable | \(F_i(t)=q_i\mathbf{E}(t)\), field energy \(U(t)=-\sum_i q_i\mathbf{E}(t)\cdot\mathbf{r}_i^{\mathrm{u}}\) vs time |
| Geometry | 3 point charges (ions) in a cubic metal-unit box; no pair force for the pin |
| Oracle | **LAMMPS** `fix efield` with equal-style variable (mandatory; feature exists in LAMMPS) |
| Regression | `regressions/md/oscillating-efield-md/` |
| Queue row | `docs/developer/metalwalls-method-queue.md` `## in_progress` `[~] oscillating-efield-md` |
| Depends | constant `FixEfield` foothold (`external-efield-md` ✓) |

## Goal

Ship a public CPU component that applies a **time-dependent uniform** electric
field

\[
\mathbf{E}(t) = \mathbf{E}_0\cos\!\bigl(2\pi\, t / T + \varphi\bigr),\qquad
t = n_{\mathrm{step}}\,\Delta t
\]

to charged atoms in classical MD, matching LAMMPS `fix efield` with an
equal-style variable for the field components (metal units) on forces and on
the instantaneous field energy (field frozen over the step; same \(-\mathbf{F}\cdot\mathbf{r}^{\mathrm{u}}\) as the constant path).

This is the **E3 residual** after constant efield ✓ and finite-field
single-electrode ✓.

## Non-goals

- Atom-style (spatially varying) variable fields / region keyword.
- Dipole / torque path.
- Full general equal-style expression engine (arbitrary `v_` formulas beyond AC).
- CUDA path; const-pot electrode solve (E1).
- Public `*Config` bags or `make_*` / `build_*` factories.
- New Driver category (component on existing MD / Evaluate only).

## Public API (binding shape)

```cpp
// Amplitudes [V/Å], period T [fs], timestep dt [fs], phase φ [rad] (default 0)
cpu::md::FixOscillatingEfield osc(
    /*Ex0=*/0.0, /*Ey0=*/0.0, /*Ez0=*/0.2,
    /*period=*/1000.0, /*dt=*/1.0, /*phase=*/0.0);

// MD (default_bindings → post_force):
cpu::md::MD sim;
sim.add(cpu::md::FixOscillatingEfield{0.0, 0.0, 0.2, 1000.0, 1.0});

// Instantaneous field helpers (read-only):
double Ez = osc.field_z(/*step=*/125);  // Ez0 * cos(2π step dt / T + φ)
```

- **No** `OscillatingEfieldConfig` / nested `::Config`.
- Parameters are **constructor args** + public members.
- Kind: `cpu::Bias` (external energy/force).
- Units (public surface): \(\mathbf{E}\) **V/Å**, \(q\) in \(e\), \(\mathbf{F}\)
  **eV/Å**, \(U\) **eV**, \(T\) and \(\Delta t\) **fs**, \(\varphi\) **rad**.
  Metal `qe2f = 1`: \(\mathbf{F}=q\mathbf{E}(t)\),
  \(U=-\sum q\,\mathbf{E}(t)\cdot\mathbf{r}^{\mathrm{u}}\).

## Physics (LAMMPS `fix efield` equal-style variable)

At MD step \(n\) with \(t = n\Delta t\):

\[
\mathbf{E}(t)=\mathbf{E}_0\cos(2\pi t/T+\varphi),\qquad
\mathbf{F}_i \mathrel{+}= q_i\mathbf{E}(t),\qquad
U \mathrel{+}= -q_i\mathbf{E}(t)\cdot\mathbf{r}_i^{\mathrm{u}}.
\]

Orthorhombic image unwrap only (same foothold as `FixEfield`). Writes
`kAtom_f*` and `meta.get(kEnergy)`.

LAMMPS twin (metal units; \(T=1\,\mathrm{ps}=1000\,\mathrm{fs}\),
\(\Delta t=0.001\,\mathrm{ps}=1\,\mathrm{fs}\)):

```text
variable Ez equal "0.2*cos(2.0*PI*time/1.0)"
fix efield 0.0 0.0 v_Ez potential v_phi   # phi atom = -v_Ez*z for U pin
```

(Force is always applied for variable efield; energy pin uses the
`potential` atom-style path so \(U=\sum q\phi\) matches \(-\sum q E_z z\).)

## Tasks

1. Write this spec + `acceptance.md` on `wf/oscillating-efield-md` **before** production code.
2. Offline LAMMPS pin → `regressions/md/oscillating-efield-md/REFERENCE.md`
   (`oracle: lammps`, `oracle_version: …`) + input under
   `tools/metalwalls_oracle/inputs/lammps/oscillating-efield-md/`.
3. Implement `cpu::md::FixOscillatingEfield`
   (`include/…/dynamics/fix_oscillating_efield.h`,
   `src/cpu/dynamics/fix_oscillating_efield.cpp`); wire CMake under `ATV_BUILD_MD`.
4. Unit test `tests/cpu/dynamics/test_fix_oscillating_efield.cpp`
   (\(E(t)\) schedule, \(F=qE(t)\), \(U\) formula).
5. E2E regression `reg_oscillating-efield-md` vs offline multi-step pin
   (dual force + energy at several \(t\)).
6. MD-only build + run until PASS.
7. Perf gate vs same LAMMPS binary (or `perf_exempt` with reason).
8. Commit on `wf/oscillating-efield-md` only.

## Acceptance criteria

See `acceptance.md` (binding).
