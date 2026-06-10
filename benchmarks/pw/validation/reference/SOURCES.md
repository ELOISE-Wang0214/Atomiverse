# Reference data provenance — pw / Si

## Si_lda_nc.upf

ONCVPSP norm-conserving scalar-relativistic LDA (PZ) pseudopotential
for Si, **without** core correction, from the PseudoDojo project
(v0.4 tables):

- Source: <https://github.com/PseudoDojo/ONCVPSP-LDA-SR-noCCv0.4>
  (file `Si/Si.upf`, retrieved 2026-06-10)
- Generator: ONCVPSP 3.3.0 (D. R. Hamann), lmax = 2, 2 projectors per
  channel, Z_val = 4.
- Cite: D. R. Hamann, *Phys. Rev. B* **88**, 085117 (2013);
  M. J. van Setten et al., *Comput. Phys. Commun.* **226**, 39 (2018).
- The noCC variant is used because the pw backend applies no nonlinear
  core correction in the XC term.

## cp2k/ decks

`si_bulk.inp.template` inlines:

- Basis `Si DZVP-GTH-q4` from cp2k `data/GTH_BASIS_SETS`,
- Pseudopotential `Si GTH-PADE-q4` from cp2k `data/GTH_POTENTIALS`,

both retrieved from <https://github.com/cp2k/cp2k> (master,
2026-06-10). Run `cp2k/run_eos.sh /path/to/cp2k` to produce
`cp2k/eos_reference.csv`; official static binaries are at
<https://github.com/cp2k/cp2k/releases>.

## Literature LDA fallback (used until eos_reference.csv exists)

LDA equilibrium properties of Si (diamond structure):
a0 ≈ 5.40 Å, B0 ≈ 96 GPa — e.g. R. J. Needs & A. Mujica,
*Phys. Rev. B* **51**, 9652 (1995) (a0 = 5.396 Å, B0 = 95.9 GPa);
experimental: a0 = 5.431 Å, B0 = 98.8 GPa.

## Comparability notes

- CP2K Quickstep is GPW (Gaussian orbitals + plane-wave density),
  Atomiverse pw is pure plane waves: **absolute energies are not
  comparable**; the benchmark compares EOS-derived observables
  (a0, B0) fitted with the same quadratic-in-V code path.
- Pseudopotentials differ by construction (ONCV-PZ vs GTH-PADE), and
  the LDA parametrizations differ slightly (PW92-limit vs PADE);
  both effects are ≪ the PASS tolerances (0.08 Å / 25 GPa).
