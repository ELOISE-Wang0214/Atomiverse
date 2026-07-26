# Spec — dos-spin-si

**Slug:** `dos-spin-si`  
**Method:** Collinear-spin PW DOS (diamond Si)  
**Domain:** `pw`  
**Branch:** `wf/dos-spin-si`  
**Oracle:** offline CP2K note (`tools/cp2k_oracle/inputs/dos-spin-si.note`)  
**Status:** in progress (wave 11)

## Gap

`cpu::pw::Dos::compute` throws on `method.spin_polarized` after
`bands-spin-si` closed collinear spin **bands**. Residual: collinear-spin
**DOS** (Gaussian-broadened total + optional per-spin channels).

Roadmap: `docs/roadmap/dft.md` § Band structure / DOS residual
`dos-spin-si`. Queue: `docs/developer/cp2k-method-queue.md`.

## Public API (no `*Config`, no factories)

```cpp
cpu::pw::PW method;
method.spin_polarized = true;
method.multiplicity = 1;   // M = 0 → N↑ = N↓ = 4 for Si
method.upf_path = "..."; method.ecut_ry = 15; method.k_mesh = 2;
method.xc = cpu::scf::XCType::LSDA;

cpu::pw::Dos dos;
dos.k_mesh = 4;
dos.n_bands = 8;
dos.sigma_ry = 0.05;
const DosResult dr = dos.compute(method, frame);
// dr.nspin == 2; ∫^{VBM} D(E) dE ≈ n_e; gap open; |M| ~ 0
```

## Algorithm

1. Fixed-occupation collinear SCF on `method.k_mesh` (reuse
   `run_collinear_scf` / same path as `Bands` spin).
2. Freeze \(V_{\mathrm{eff},\sigma}\) from mixed \(\rho_\sigma\).
3. Rebuild denser Monkhorst–Pack mesh (`Dos::k_mesh`).
4. NSCF diagonalization per spin on frozen \(V_{\mathrm{eff},\sigma}\).
5. Joint α/β gap from fixed occupations (\(N_\uparrow,N_\downarrow\)).
6. Gaussian DOS: spin-orbital degeneracy **1** (not closed-shell 2);
   total \(D(E)\) = \(D_\alpha + D_\beta\); optional per-spin curves.
7. Occupied integral to joint VBM ≈ \(n_e\).

**Out of scope this slug:** spin + metallic smearing (reject with clear
message; residual for a later metallic-spin DOS if needed). Metallic
restricted DOS remains `dos-metallic-si`.

## Acceptance

See `acceptance.md`. Regression: `regressions/pw/dos-spin-si/`.

## Perf

Prefer fair CP2K wall when binary available; else `perf_exempt` with
reason (host has no CP2K module — same policy as `bands-spin-si` /
`dos-metallic-si`).
