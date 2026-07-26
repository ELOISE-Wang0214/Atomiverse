# gpw-xdm-ice — Exchange-hole dipole moment (XDM) free-atom dispersion

## Goal

Ship **free-atom XDM** (Becke–Johnson exchange-hole dipole moment model) as an
additive geometry-only dispersion correction on Γ-point GPW, exposed as
`GpwDft::xdm` / `GPWEvaluator::xdm`, reporting `GpwResult::e_disp` [Hartree].

Ice-relevant H-bonded **water dimer** in a periodic 12 Å cell (same geometry as
`gpw-ts-ice` / `gpw-pbe-d3bj-ice`).

## Scope

**In**

- Free-atom multipole moments \(M_1,M_2,M_3\) + CRC free-atom polarizabilities
  \(\alpha\) for H and O (in-tree tables; no file I/O).
- Pairwise XDM \(C_6,C_8,C_{10}\) combination rules + BJ damping with
  **QE NC PBE** \((a_1,a_2)=(0.4283,\,2.4690\,\text{Å})\).
- Periodic pure real-space lattice sum (orthorhombic, `pair_image_counts`,
  rcut 50 Å) and isolated path.
- Analytic forces on the free-atom energy (FD-gated in unit tests).
- Mutual exclusion vs `d3` / `d3bj` / `d4` / `ts` / `mbd` / `vdw_df` /
  `vdw_df2` / `rvv10`.
- Offline pin (dual pure-Python lattice sum) → regression hard gate.

**Out**

- Density-dependent Hirshfeld volumes / BR hole on the SCF density (follow-on;
  free-atom foothold only — same policy as free-atom TS before Hirshfeld).
- C9 three-body term.
- Stock CP2K `TYPE XDM` (does not exist; PAIR_POTENTIAL is only
  DFTD2|DFTD3|DFTD3(BJ)|DFTD4).
- CUDA port; molecular Gaussian-path XDM.

## API

```cpp
cpu::gpw::GpwDft gpw;
gpw.xdm = true;   // exclusive vs other dispersion flags
// ... build/run → res.e_disp [Ha], folded into e_total
```

## Science

Becke & Johnson, *J. Chem. Phys.* **127**, 154108 (2007); Otero-de-la-Roza &
Johnson, *J. Chem. Phys.* **136**, 174109 (2012). Free-atom \(\alpha\) from CRC
(postg/QE); free-atom moments inverted from literature homo \(C_n\) so
identical-atom \(C_6,C_8,C_{10}\) match Yan/Thakkar (H) and Chu–Dalgarno (O)
with the XDM combination rules.

## Acceptance

See `acceptance.md`.
