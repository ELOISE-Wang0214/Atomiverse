# Reference — plane-wave true PBE GGA on bulk Si

## Method

Spin-unpolarized plane-wave Kohn–Sham total energy and ionic forces with the
**true GGA PBE** exchange–correlation functional (Perdew, Burke, Ernzerhof,
*Phys. Rev. Lett.* **77**, 3865 (1996)). Density gradients enter the XC energy
and potential spectrally on the FFT grid (not the historical `σ = 0`
LDA-limit call of `evaluate_xc(PBE, …)`).

Public API: concrete `cpu::pw::PW` with members (`upf_path`, `ecut_ry`,
`k_mesh`, `xc = XCType::PBE`, …); Frame I/O; `pw_forces_analytic` for forces
[eV/Å]. No public `*Config` bag.

## System

Silicon diamond (FCC primitive cell, 2 atoms):

| | |
| --- | --- |
| Lattice constant | \(a = 5.43\) Å (standard experimental Si) |
| Basis | \((0,0,0)\), \((a/4,a/4,a/4)\) |
| Pseudopotential | in-tree ONCVPSP NC UPF `tests/data/pw/Si_lda_nc.upf` (scalar-relativistic, LDA-generated; used here only as a fixed NC valence PP for the GGA path smoke) |
| Cutoff | \(E_{\mathrm{cut}} = 15\) Ry |
| k-mesh | \(2\times 2\times 2\) Monkhorst–Pack |
| Mixing | Pulay, `mix_beta = 0.5`, `tol = 1e-5`, `max_iter = 40` |

## Offline oracle status

**No CP2K or Quantum ESPRESSO binary** was available on the worker host at pin
time (`module avail cp2k` / `qe` empty; no `pw.x` / `cp2k.psmp` on `PATH`).
Absolute total-energy digits are therefore **not** invented from memory and
**not** taken from a live third-party process at ctest time.

Pin strategy (binding for this regression):

1. **Crystal-symmetry force pin (literature physics).** On ideal diamond
   sites the Hellmann–Feynman force on every atom vanishes by lattice
   symmetry. Residual after a modest SCF is bounded by
   \[
   |F_i| < 0.5~\mathrm{eV/Å}
   \]
   for both Si atoms (same bound style as `regressions/pw/h_force_analytic`).

2. **GGA activation pin.** The same geometry is also run with
   `xc = XCType::LSDA`. True PBE GGA must shift the total energy by at least
   \[
   |E_{\mathrm{PBE}} - E_{\mathrm{LSDA}}| > 10^{-4}~\mathrm{Ry}
   \]
   (documents that the production path is not stuck in a pure density-only
   LDA limit of a single functional call). Unit tests separately prove
   PBE@\(σ=0\) vs true GGA on a nonuniform density.

3. **Energy finiteness / convergence.** \(E_{\mathrm{PBE}}\) is finite;
   preferred: SCF converges within `max_iter`.

When a CP2K Quickstep or QE twin becomes available offline, replace/extend
this file with the twin total energy (cite version, cutoff, PP) without
changing the public API.

### Observed Atomiverse numbers (informational, not a self-golden pin)

Worker run 2026-07-23 (`wf/pw-pbe-gga-si`, Release, single rank):

| Quantity | Value |
| --- | --- |
| \(E_{\mathrm{PBE}}\) | −15.70481672 Ry (converged, 6 iters) |
| \(E_{\mathrm{LSDA}}\) | −15.67175235 Ry (converged, 5 iters) |
| \(\|E_{\mathrm{PBE}}-E_{\mathrm{LSDA}}\|\) | 0.0330644 Ry (> 1e−4 floor) |
| \(\|F_0\|\), \(\|F_1\|\) | ~0.0015 eV/Å (< 0.5 eV/Å) |
| Wall (PBE+forces+LSDA) | ~27 s |

## Performance gate

**perf_exempt = true**

Reason: no offline CP2K/QE binary on host for a fair wall twin
(\(T_{\mathrm{atv}}/T_{\mathrm{cp2k}}\) undefined). Correctness pins above still
apply. Re-run fair bench when oracle is installed; target ratio ≤ 10.

## What the binary checks

`reg_pw-pbe-gga-si`:

1. `PW{xc=PBE}` SCF on Si diamond → finite energy [Ry].
2. `pw_forces_analytic` → \(|F_i| < 0.5\) eV/Å for both atoms.
3. Second run `xc=LSDA` → finite; \(|E_{\mathrm{PBE}}-E_{\mathrm{LSDA}}| > 1e-4\) Ry.
4. No external process spawn.

## Provenance

- Spec: `.claude/specs/pw-pbe-gga-si/spec.md`
- Functional: PBE (1996); shared `cpu::scf::evaluate_xc` (libxc-parity tree)
- Geometry: textbook Si diamond \(a=5.43\) Å
- PP: Hamann ONCVPSP NC UPF shipped under `tests/data/pw/`
