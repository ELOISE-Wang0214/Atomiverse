# Spec — `dfpt-finite-q-si`

## Intent

Ship **finite-q** density-functional perturbation theory (DFPT) phonons for the
plane-wave method: self-consistent Sternheimer response at a user-chosen
fractional wavevector \(\mathbf{q}\), producing \(\omega(\mathbf{q})\) [cm⁻¹]
and the complex Hermitian dynamical matrix. Closes the follow-on declared by
`dfpt-phonon-si` (Γ-only foothold) and the roadmap row “finite-q DFPT
Sternheimer ○”.

Orthogonal to `phonon-dispersion-si` (IFC Fourier of classical SW supercell
Φ) — this path is **analytic DFPT** on the PW LSDA Hamiltonian.

## Public API (no `*Config`, no factories)

```cpp
cpu::pw::DfptPhonon dfpt;
dfpt.q_frac = {0.5, 0.0, 0.0};   // fractional crystal coords of the cell
// response_tol / cg_tol / mix_beta — same members as Γ
auto r = dfpt.compute(pw, frame); // ω [cm⁻¹] ascending; D(q) → λ → ω
```

| Symbol | Role |
| ------ | ---- |
| `cpu::pw::DfptPhonon::q_frac` | `std::array<double,3>` fractional \(\mathbf{q}\) (default `{0,0,0}` = Γ → existing path) |
| `cpu::pw::DfptPhononResult` | `frequencies` [cm⁻¹]; `force_constants` stores **Re** \(C(\mathbf{q})\) in eV/Å² for diagnostics; new optional imag diagnostic if needed |
| Internal | `internal::dfpt_q_force_constants`, `ewald_dynamical_matrix_q` |

Frame is the only I/O; workspace stays private to the method.

## Physics (Baroni et al., RMP 73, 515 (2001))

**Γ (`q_frac = 0`):** existing `dfpt_gamma_force_constants` (Sternheimer PCG,
Hartree + LDA kernel, analytic Ewald/KB).

**Commensurate finite q (foothold product path):** a rational
\(\mathbf{q}=(h/n_a,k/n_b,l/n_c)\) is the **Γ point of the
\(n_a\times n_b\times n_c\) supercell**. Evaluate:

1. `expand_supercell(frame, na, nb, nc)`
2. `dfpt_gamma_force_constants` on the supercell (same Sternheimer engine)
3. IFC Fourier \(D(\mathbf{q})\) via `PhononDispersion::frequencies_at_q`
   (identical folding as the classical SW path, Φ from DFPT)

This is analytically equivalent to finite-q DFPT for commensurate \(\mathbf{q}\)
(Born–von Karman). Direct incommensurate Sternheimer (`dfpt_q_force_constants`
scaffolding + `ewald_dynamical_matrix_q`) remains internal follow-on.

Scope (same as Γ): NC-UPF, `XCType::LSDA`, fixed occupations, restricted,
no metallic smearing. Throws `std::invalid_argument` otherwise. Supercell
product \(n_a n_b n_c\le 8\).

## System / regression twin

| Setting | Value |
| ------- | ----- |
| Structure | diamond Si, 2-atom fcc primitive, \(a=5.43\) Å |
| PP | `tests/data/pw/Si_lda_nc.upf` |
| XC | LSDA |
| Cutoff | 10 Ry (minimal twin; matches `dfpt-phonon-si`) |
| k-mesh | \(2\times2\times2\) MP |
| \(\mathbf{q}\) | \(X=(1/2,0,0)\) fractional of the **primitive** reciprocal cell |

## Acceptance gates (offline pins only)

1. **Γ recovery:** `q_frac={0,0,0}` reproduces `reg_dfpt-phonon-si` optical
   window / C00 energy-FD pin (no regression of Γ path).
2. **X-point spectrum:** 6 real finite frequencies; literature LDA Si X-point
   windows (Giannozzi et al., PRB 43, 7231 (1991)) with the same generous
   minimal-settings bias as Γ (+~12% class) — see `REFERENCE.md`.
3. **Hermiticity:** max \(|C-C^\dagger|/2\) below pin.
4. **Offline energy-FD (frozen phonon):** second derivative of SCF energy in a
   \(2\times1\times1\) supercell for one symmetry-allowed X pattern, hard-coded
   pin — never live CP2K/QE in ctest.
5. **Perf:** `perf_exempt` if no CP2K DFPT wall; else \(T_\mathrm{atv}/T_\mathrm{oracle}\le 10\).

## Out of scope

- Metallic smearing / UKS DFPT
- GGA response kernel
- LO–TO non-analytic correction (Si is non-polar; not needed)
- CUDA port
- Full Brillouin-zone dispersion product (caller loops `q_frac`)

## Files

| Path | Change |
| ---- | ------ |
| `include/.../dfpt_phonon.h` | `q_frac` member; docs |
| `include/.../internal/dfpt.h` | finite-q API |
| `src/.../internal/dfpt.cpp` | Ewald-q + Sternheimer-q |
| `src/.../dfpt_phonon.cpp` | dispatch Γ vs q; Hermitian mass-weight |
| `tests/cpu/planewave/test_dfpt.cpp` | Ewald-q ASR/FD; guards; q=0 dispatch |
| `regressions/response/dfpt-finite-q-si/` | e2e + `REFERENCE.md` |
| `regressions/CMakeLists.txt` | register reg |
| `tools/cp2k_oracle/inputs/` | optional note/input when CP2K available |
