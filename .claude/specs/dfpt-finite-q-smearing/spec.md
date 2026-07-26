# Spec — `dfpt-finite-q-smearing`

## Goal

Productize **finite-q DFPT phonons** for the plane-wave method when the
ground state uses **restricted metallic occupation smearing**
(Fermi–Dirac / Methfessel–Paxton / Gaussian).

Today:

- Γ DFPT + restricted smear is product (`reg_dfpt-metallic-si`).
- Finite-q fixed-occ is product (`reg_dfpt-finite-q-si` commensurate X;
  `reg_dfpt-incommensurate-si` direct Sternheimer).
- `dfpt_q_force_constants` **throws** on `method.smearing != none`
  ("use Γ restricted smear; finite-q follow-on").

Close the residual by accepting restricted smear on:

1. **Commensurate** finite \(\mathbf{q}\) (supercell Γ-DFPT + IFC Fourier —
   reuses the metallic Γ engine already shipped).
2. **Direct** finite-\(q\) Sternheimer (`dfpt_q_force_constants`) with
   fractional \(f_{nk}\) weights (no \(\delta\mu\) at \(q\neq 0\)).

Scope: **restricted** (closed-shell) PW + pure LSDA/GGA + UPF.
Collinear **spin + smear** remains rejected (`dfpt-spin-smearing`).

## Public API

No new public types. Members only:

```cpp
cpu::pw::PW pw;
pw.upf_path = "Si.upf";
pw.ecut_ry = 10.0;
pw.k_mesh = 1;
pw.xc = cpu::scf::XCType::LSDA;
pw.smearing = cpu::pw::Smearing::fermi_dirac;
pw.degauss = 0.02;   // Ry
pw.extra_bands = 0;  // gapped Si twin

cpu::pw::DfptPhonon dfpt;
dfpt.q_frac = {0.5, 0.0, 0.0}; // X → 2×1×1 supercell path
auto r = dfpt.compute(pw, frame);
```

No `*Config`, no factories. Spectra on `DfptPhononResult` (Frame purity).

## Algorithm

### Commensurate path (product twin)

Rational \(\mathbf{q}=(h/n_a,\ldots)\) with \(n_a n_b n_c\le 8\): expand
supercell → `dfpt_gamma_force_constants` (already smear-aware) →
`PhononDispersion::frequencies_at_q`.

### Direct Sternheimer path

1. Ground-state SCF with `run_scf_smear_restricted` when smear is set
   (\(n_b=\lceil N_e/2\rceil+\mathrm{extra\_bands}\)).
2. Freeze \(V_{\mathrm{eff}}\); eigenstates at every \(\mathbf{k}\) and
   \(\mathbf{k}+\mathbf{q}\).
3. Shared \(\mu\) from the frozen \(\mathbf{k}\) spectrum; fill
   \(f_{n\mathbf{k}}\in[0,2]\) and \(f_{n,\mathbf{k}+\mathbf{q}}\).
4. Response only for bands with \(f_{n\mathbf{k}}\ge 10^{-8}\); Sternheimer
   projector at \(\mathbf{k}+\mathbf{q}\) removes bands with
   \(f_{n,\mathbf{k}+\mathbf{q}}\ge 10^{-8}\).
5. Weight \(\Delta n\) and force-constant terms by \(f_{n\mathbf{k}}\).
6. **No** de Gironcoli \(\delta\mu\) at \(q\neq 0\) (charge neutrality is
   automatic for the lattice-periodic response at finite \(q\)).

Cite: Baroni et al., *Rev. Mod. Phys.* **73**, 515 (2001);
de Gironcoli, *Phys. Rev. B* **51**, 6773 (1995).

## Acceptance

| Gate | Criterion |
| ---- | --------- |
| Unit | Commensurate / direct finite-q + FD smear passes XC/smear gates (UPF I/O may still throw on missing file); spin+smear still throws; hybrid / empty UPF still throw |
| Regression | `regressions/response/dfpt-finite-q-smearing` — diamond Si at X with FD smear: \(\lvert\omega\rvert\) windows + asymmetry (sibling of `dfpt-finite-q-si`) |
| Dual gate | Correctness + `within_10x` or `perf_exempt` |
| Arch | No public Config/factory; Frame I/O only |

## Out of scope

- Collinear spin + metallic smear DFPT
- True-metal Al finite-q DFPT with large `extra_bands` / dense k
- Live CP2K in ctest (offline pins only)

## Deps

- `dfpt-metallic-si` ✓
- `dfpt-finite-q-si` ✓
- `dfpt-incommensurate-si` ✓
