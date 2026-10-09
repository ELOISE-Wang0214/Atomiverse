# Reference — DFTB0 analytic nuclear forces vs DFTB+ non-SCC

```text
oracle: dftb+
oracle_version: 24.1
```

## Method

Non-SCC DFTB / **DFTB0**. Energy is the band energy of \(H^0\) plus the
pair repulsive,

\[
E = \sum_i 2\varepsilon_i + E_\mathrm{rep}.
\]

Forces are the closed-form gradient of that energy at the non-SCC density
(Hellmann–Feynman + Pulay + repulsive). There is no Hubbard \(\gamma\), no
\(\Delta q\), and no SCC shift on the overlap. The pair walk is one
minimum image, the same geometry the energy uses. Not a \(6N\) finite
difference.

Reference: Elstner et al., *Phys. Rev. B* **58**, 7260 (1998),
DOI 10.1103/PhysRevB.58.7260, the zeroth-order band + repulsive piece,
before the second-order charge fluctuation. Force-module boundary:
Hourahine et al., *J. Chem. Phys.* **152**, 124101 (2020),
DOI 10.1063/1.5143190. **Not** OpenMOPAC. **Not** xtb.

## Offline oracle run

DFTB+ **24.1**, mio-1-1, `SCC = No`, `PrintForces = Yes`. Nothing is
executed at ctest time. Inputs and archived `detailed.out` transcripts:

| Case | Path |
| ---- | ---- |
| H₂, \(R = 0.74\) Å | `tools/tb_oracle/inputs/dftbplus/dftb0-analytic-forces/h2_074/` |
| H₂, \(R = 0.90\) Å | `tools/tb_oracle/inputs/dftbplus/dftb0-analytic-forces/h2_090/` |
| off-eq H₂O | `tools/tb_oracle/inputs/dftbplus/dftb0-analytic-forces/h2o/` |
| tilted CO | `tools/tb_oracle/inputs/dftbplus/dftb0-analytic-forces/co/` |

H₂ at 0.74 Å is also the homonuclear reduction of the SCC transcript in
`h2_mio_scc/` (Δq = 0). The numbers below are from a fresh **non-SCC**
run, not from that reduction. Water and CO are heteronuclear: an SCC run
is a different Hamiltonian and is not a DFTB0 pin.

Conversion, from `core/units.h`: 1 Ha = 27.211386245988 eV,
1 Ha/Bohr = 51.42206747467393 eV/Å. DFTB+ prints forces in atomic units
(`Maximal derivative component: … au`).

### H₂, \(R = 0.74\) Å (cubic isolation 30 Å)

| Quantity | Transcript | Unit |
| -------- | ---------- | ---- |
| Total energy | −0.6749509324 | Ha |
| \(F_{1x}\) | −0.001976354621 | Ha/Bohr |
| \(F_{2x}\) | +0.001976354621 | Ha/Bohr |

The published eV figure **−18.36635052 eV** (`reg_dftb0-h2`) must not move.

### H₂, \(R = 0.90\) Å

| Quantity | Transcript | Unit |
| -------- | ---------- | ---- |
| Total energy | −0.6630166170 | Ha |
| \(F_{1x}\) | +0.065975227977 | Ha/Bohr |
| \(F_{2x}\) | −0.065975227977 | Ha/Bohr |

### Off-equilibrium water (Å)

| Atom | Z | x | y | z |
| ---- | - | - | - | - |
| O | 8 | 0.00 | 0.00 | 0.00 |
| H | 1 | 1.05 | 0.00 | 0.00 |
| H | 1 | −0.26 | 1.00 | 0.00 |

| Quantity | Transcript | Unit |
| -------- | ---------- | ---- |
| Total energy | −4.0964033746 | Ha |
| \(F_\mathrm{O}\) | (0.039206443065, 0.036588999897, 0) | Ha/Bohr |
| \(F_\mathrm{H1}\) | (−0.048418613962, −0.000229747051, 0) | Ha/Bohr |
| \(F_\mathrm{H2}\) | (0.009212170898, −0.036359252846, 0) | Ha/Bohr |

### Tilted CO (ppσ / ppπ and a transverse component)

| Atom | Z | x | y | z |
| ---- | - | - | - | - |
| C | 6 | 0.00 | 0.00 | 0.00 |
| O | 8 | 1.20 | 0.15 | 0.00 |

| Quantity | Transcript | Unit |
| -------- | ---------- | ---- |
| Total energy | −5.0270725023 | Ha |
| \(F_\mathrm{C}\) | (0.183901804021, 0.022987724205, 0) | Ha/Bohr |
| \(F_\mathrm{O}\) | (−0.183901804021, −0.022987724205, 0) | Ha/Bohr |

## Pins (hard-coded gates)

| Gate | Tol |
| ---- | --- |
| Total energy vs the Ha transcript × `kHartreeInEv` | \(1\times10^{-6}\) eV |
| Each force component vs the Ha/Bohr transcript × `kHartreePerBohrInEvPerAngstrom` | \(1\times10^{-5}\) eV/Å |
| Analytic \(F\) vs central FD of the same energy (\(h = 10^{-4}\) Å) | \(1\times10^{-5}\) eV/Å |
| \(\lvert\sum_i F_i\rvert\) | \(1\times10^{-8}\) eV/Å |
| Published H₂ 0.74 eV pin −18.36635052 | \(1\times10^{-6}\) eV |

**Tolerance rationale.** The transcript prints energy to \(10^{-10}\) Ha
(\(\approx 3\times10^{-9}\) eV) and forces to \(\sim10^{-12}\) Ha/Bohr
(\(\approx 5\times10^{-11}\) eV/Å). The \(10^{-6}\) eV / \(10^{-5}\) eV/Å
gates sit above that print precision and above the central-difference
floor of `MioSkTable::grid_derivative` (\(h = 10^{-5}\) Å on the hopping
interpolant). They are not fitted to an ATV residual. A factor-of-two or
a sign error in the Pulay contraction is eV/Å, not \(10^{-5}\).

A wrapped H₂ (atoms at 0.20 Å and 7.60 Å in an 8 Å cube) is an ATV
min-image identity only: DFTB+ `C` geometry does not wrap. It gates
\(F = -\nabla E\) and \(\sum F = 0\), not an absolute oracle number.

## Gaps

- No cell stress, no k-mesh, no lattice-sum images beyond the single
  minimum image the energy already uses.
- Closed shell, integer occupation, H/C/N/O only.
- The SK radial derivative is still a finite difference of the interpolant
  (`dftb-sk-analytic-derivative`); this slug does not change that.

## Perf

Warm DFTB+ 24.1 single point, `OMP_NUM_THREADS=1`, process wall
(five repeats, 2026-10-09): H₂ \(R=0.90\) Å median **0.03 s**, water
median **0.04 s**. Cold start of the binary was 0.92 s and is not the
comparison.

ATV `Evaluate::run`, same thread count, Release, in-process chrono in
`reg_dftb0-analytic-forces`: first H₂ launch **0.014 s** (runtime init),
later H₂ **< 0.1 ms**, water **0.0010 s**, CO **0.0003 s**.

\(T_\mathrm{atv}/T_\mathrm{dftb+} = 0.014/0.03 = 0.45\) using the cold
ATV launch against the warm DFTB+ process (the slower reading). Warm
repeats are ~0.03. **Within 10×.** Not `perf_exempt`.

## Observed (2026-10-09, Release, 1 thread)

Not a retune of the gates. Max |ΔE| vs the Ha transcript was
\(8.3\times10^{-7}\) eV (CO). Max |ΔF| stayed under the
\(1\times10^{-5}\) eV/Å gate (H₂ 0.74 Å, ~2×10⁻⁶ eV/Å). FD identity
max residual \(1.0\times10^{-6}\) eV/Å (wrapped H₂). H₂ 0.74 Å energy
−18.3663505230 eV, |Δ| 3×10⁻⁹ eV from the published pin.
