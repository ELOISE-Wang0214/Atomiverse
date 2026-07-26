# Spec — `periodic-wannier-uks`

| Field | Value |
| ----- | ----- |
| Slug | `periodic-wannier-uks` |
| Loop | CP2K capability (`cp2k-method-loop`) |
| Oracle | offline CP2K / literature geometry (no live oracle in ctest) |
| Domain | response / PW MLWF |
| Deps | `periodic-wannier-si` ✓ |
| Branch | `wf/periodic-wannier-uks` |
| Method | response |

## Goal

Ship **collinear spin-polarized** periodic multi-k MLWF (Marzari–Vanderbilt)
on the plane-wave method.

Today:

- Closed-shell restricted MLWF works (`PeriodicWannier` + `reg_periodic-wannier-si`).
- `PeriodicWannier::compute` **throws** on `method.spin_polarized`.

Target: remove the spin throw; localize **α and β occupied manifolds
independently** on frozen collinear \(V_{\mathrm{eff},\sigma}\), fixed
occupations, insulator only (no metallic smearing).

## Public API (no `*Config`, no factories)

```cpp
cpu::pw::PW method;
method.upf_path = "Si.upf";
method.ecut_ry = 15.0;
method.k_mesh = 2;                 // SCF mesh
method.xc = cpu::scf::XCType::LSDA;
method.spin_polarized = true;
method.multiplicity = 1;           // M = N↑−N↓ = 0 (closed-shell collinear)

cpu::pw::PeriodicWannier mlwf;
mlwf.k_mesh = 2;                   // MLWF MP mesh (may differ from SCF)
auto r = mlwf.compute(method, frame);
// r.n_wann_up / r.n_wann_dn; r.centers [Å] = α then β;
// r.spin[i] ∈ {0,1}; r.spreads [Å²]
```

Parameters remain members on `PeriodicWannier` / `PW`. Frame is sole geometry
I/O; spectra/centers live on `PeriodicWannierResult` only.

## Physics

Collinear Marzari–Vanderbilt (Marzari & Vanderbilt, *Phys. Rev. B* **56**,
12847 (1997); Marzari et al., *Rev. Mod. Phys.* **84**, 1419 (2012)) per spin:

1. Fixed-occupation collinear SCF → \(n_\uparrow,n_\downarrow\) and
   \(V_{\mathrm{eff},\sigma}\).
2. NSCF Bloch states on the MLWF Monkhorst–Pack mesh at frozen
   \(V_{\mathrm{eff},\sigma}\) for each spin with \(n_\sigma>0\).
3. Neighbor overlaps \(M^{(\mathbf k,\mathbf b)}\), projection init (optional),
   MV steepest-descent, Berry-phase centers — **identical pure-math path** as
   the restricted foothold, with occupation weight \(f=1\) (spin-orbital).

Diamond Si with multiplicity 1 (\(N_\uparrow=N_\downarrow=4\)) must place **both**
spin manifolds on the four tetrahedral Si–Si bond midpoints (same geometric
pin as `periodic-wannier-si`; spin is a framework productization, not a new
geometry).

## Acceptance

1. **Unit** (`test_cpu_pw_periodic_wannier`): spin-polarized gate no longer
   throws the foothold message; invalid `k_mesh`/`tol` still throw; pure math
   helpers unchanged.
2. **Regression** `regressions/pw/periodic-wannier-uks/`:
   diamond Si, \(a=5.43\) Å, LSDA, in-tree Si NC-UPF, ecut 15 Ry, SCF k=2,
   MLWF k=2, `spin_polarized=true`, multiplicity=1;
   `n_wann_up==n_wann_dn==4`; each spin’s four centers assign to bond midpoints
   within **0.90 Å** min-image; per-spin \(\Omega_\mathrm{final}<0.25\,\Omega_\mathrm{init}\);
   SCF converged; wall ≤ 60 s.
3. **Offline pin:** `tools/cp2k_oracle/inputs/si_periodic_wannier_uks.note`
   + `REFERENCE.md` (CP2K `&DFT` spin + `&LOCALIZE` note; absolute spreads
   differ GTH/MOLOPT vs ONCV; bond geometry is the shared content).
4. **Perf:** `T_atv/T_oracle ≤ 10` **or** `perf_exempt=true` when no fair
   CP2K PW/NC-UPF wall twin on host.
5. No public `*Config` / `make_*` / `build_*`; public Frame energy still eV
   when written (MLWF energies stay on result as Ry/Ha when cited).

## Out of scope

- Metallic smearing MLWF / disentanglement.
- Non-collinear / SOC Wannier.
- GPW-grid MLWF; CUDA; molecular post-HF; classical MW MD.
- Open-shell geometric pins beyond Si \(M=0\) (O₂ product pin is follow-on).

## Files (expected)

| Path | Role |
| ---- | ---- |
| `include/atomiverse/cpu/planewave/periodic_wannier.h` | spin result fields + docs |
| `src/cpu/planewave/periodic_wannier.cpp` | collinear SCF + per-spin MV |
| `tests/cpu/planewave/test_periodic_wannier.cpp` | unit gates |
| `regressions/pw/periodic-wannier-uks/` | e2e + REFERENCE |
| `tools/cp2k_oracle/inputs/si_periodic_wannier_uks.note` | offline note |
| `regressions/CMakeLists.txt` | reg target + UPF define |

## Done when

`acceptance_met=true`: scoped unit + `run_reg_periodic-wannier-uks` green;
commit on `wf/periodic-wannier-uks`; build wiped.
