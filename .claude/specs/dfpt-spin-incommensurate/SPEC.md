# Spec — `dfpt-spin-incommensurate`

| Field | Value |
| ----- | ----- |
| Slug | `dfpt-spin-incommensurate` |
| Loop | CP2K capability (`cp2k-method-loop`) |
| Oracle | offline CP2K / same-Hamiltonian ATV twin (no live oracle in ctest) |
| Domain | response / PW DFPT |
| Deps | `dfpt-spin-si` ✓ + `dfpt-incommensurate-si` ✓ |
| Branch | `wf/dfpt-spin-incommensurate` |

## Goal

Ship **collinear-spin** DFPT phonons at **incommensurate** (non-rational /
oversize-supercell) \(\mathbf{q}\) on the plane-wave method.

Today:

- Γ + commensurate finite \(\mathbf{q}\): collinear spin works
  (`dfpt_gamma_spin_impl` / supercell Γ path).
- Incommensurate \(\mathbf{q}\): restricted only
  (`dfpt_q_force_constants` throws on `method.spin_polarized`).

Target: remove the spin+incommensurate throw; direct Sternheimer
`dfpt_q_force_constants` accepts collinear fixed-occupation pure LSDA /
pure GGA with \(f=1\) per spin-orbital and a spin XC kernel (same as Γ spin).

## Public API (no `*Config`, no factories)

```cpp
cpu::pw::PW method;
method.upf_path = "...";
method.ecut_ry = 10.0;
method.k_mesh = 1;
method.xc = cpu::scf::XCType::LSDA;
method.spin_polarized = true;
method.multiplicity = 1;   // M = 0 closed-shell collinear framework

cpu::pw::DfptPhonon dfpt;
dfpt.q_frac = {0.3, 0.0, 0.0};   // incommensurate under n≤6
auto r = dfpt.compute(method, frame);  // ω [cm⁻¹], C(q) [eV/Å²]
```

Parameters remain members on `DfptPhonon` / `PW`. Frame is sole I/O currency;
spectra live on `DfptPhononResult` only.

## Physics

Finite-\(\mathbf{q}\) Sternheimer (Baroni et al., RMP **73**, 515 (2001))
with collinear spin:

\[
P_c^{k+q}(H_\sigma^{k+q}-\varepsilon_{k\sigma i})\,P_c^{k+q}\,|\Delta\psi_{k\sigma i}\rangle
  = -P_c^{k+q}\,(\Delta V_{\mathrm{bare}}+\Delta V_{\mathrm{H}}[\Delta n]
     +K_{\mathrm{xc}}^{\sigma\sigma'}[\Delta n_{\sigma'}])\,|\psi_{k\sigma i}\rangle
\]

- Occupation \(f=1\) per spin-orbital.
- Hartree from total \(\Delta n=\Delta n_\uparrow+\Delta n_\downarrow\) with
  \(|G+q|\) Poisson kernel.
- Pure LSDA: local \(f_{\mathrm{xc}}^{\sigma\sigma'}\) on complex lattice-periodic
  \(\Delta\tilde n_{q\sigma}\); pure GGA: FD of `XCBuilder::step_spin` on
  Re/Im of each channel (mirrors restricted finite-q GGA path).
- Hybrids + metallic smearing remain rejected.

Closed-shell collinear (\(M=0\)) on diamond Si must match the restricted
incommensurate spectrum within numerical noise of the twin settings
(sanity gate vs `reg_dfpt-incommensurate-si` windows).

## Acceptance

1. **Unit:** `ThrowsOnSpinIncommensurateQ` removed / inverted — spin +
   `q_frac={0.3,0,0}` no longer throws at the spin gate (empty UPF still
   hits UPF guard).
2. **Regression** `regressions/response/dfpt-spin-incommensurate/`:
   diamond Si, a=5.43 Å, LSDA, NC-UPF, ecut=10 Ry, k=1×1×1,
   `spin_polarized=true`, multiplicity=1, `q_frac=(0.3,0,0)`;
   6 real modes; |\(ω\)| windows aligned with restricted incommensurate
   twin; SCF + response converge; wall ≤ 60 s.
3. **Offline pin:** `tools/cp2k_oracle/inputs/dfpt-spin-incommensurate.note`
   + `REFERENCE.md` (CP2K binary optional; same-Hamiltonian / restricted
   twin + literature-biased windows when CP2K absent).
4. **Perf:** `T_atv/T_oracle ≤ 10` **or** `perf_exempt=true` when no CP2K
   DFPT wall twin on host.
5. No public `*Config` / `make_*` / `build_*`; public energy still eV on Frame
   when energy is written (spectra stay on result object).

## Out of scope

- Metallic occupation smearing DFPT (`dfpt-metallic-si`).
- GPW DFPT; CUDA ports; molecular post-HF; classical MW MD.
- Non-collinear / SOC phonons.

## Files (expected)

| Path | Role |
| ---- | ---- |
| `src/cpu/planewave/internal/dfpt.cpp` | `dfpt_q_spin_impl` + dispatch |
| `include/atomiverse/cpu/planewave/internal/dfpt.h` | docs |
| `include/atomiverse/cpu/planewave/dfpt_phonon.h` | public docs |
| `tests/cpu/planewave/test_dfpt.cpp` | spin+incommensurate gate |
| `regressions/response/dfpt-spin-incommensurate/` | e2e + REFERENCE |
| `tools/cp2k_oracle/inputs/dfpt-spin-incommensurate.note` | offline note |
| `regressions/CMakeLists.txt` | reg target |

## Done when

`acceptance_met=true`: scoped unit + `run_reg_dfpt-spin-incommensurate` green;
commit on `wf/dfpt-spin-incommensurate`; build wiped.
