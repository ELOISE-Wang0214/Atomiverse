# H2O SCC-DFTB3 (static mio-style SK)

```text
oracle: dftb+
oracle_version: formula-Gaus2011 (DFTB+ module boundary; no host DFTB+ binary)
```

## Method

Self-consistent-charge **DFTB3** (third-order charge fluctuation):

- Gaus et al., *J. Chem. Theory Comput.* **7**, 931 (2011),
  DOI [10.1021/ct100684s](https://doi.org/10.1021/ct100684s)
- Second-order SCC backbone: Elstner et al., *Phys. Rev. B* **58**, 7260
  (1998), DOI [10.1103/PhysRevB.58.7260](https://doi.org/10.1103/PhysRevB.58.7260)

Family-locked numerical oracle package is **DFTB+** (never OpenMOPAC / xtb).

## Third-order convention (what the code computes)

With \(\gamma_{AB}\) the Ohno–Klopman-like charge kernel and
\(\Gamma_{AB}=\partial\gamma_{AB}/\partial q_A\) (on-site \(\gamma_{AA}=U_A\),
\(\Gamma_{AA}=U^{\mathrm{d}}_A\); note \(\Gamma_{AB}\neq\Gamma_{BA}\)):

\[
E^{(2)}=\tfrac12\sum_{AB}\gamma_{AB}\,\Delta q_A\Delta q_B,\qquad
E^{(3)}=\tfrac13\sum_{AB}\Delta q_A^2\,\Delta q_B\,\Gamma_{AB}
\]

which is the Gaus/Cui/Elstner (2011) form. The SCC potential is its exact
derivative,

\[
V^{(3)}_C=\tfrac23\,\Delta q_C\sum_B\Gamma_{CB}\Delta q_B
         +\tfrac13\sum_B\Delta q_B^2\,\Gamma_{BC},
\]

so \(\sum_A\Delta q_A V_A=2E^{(2)}+3E^{(3)}\) and the band-structure
double-counting correction is \(E^{(2)}+2E^{(3)}\). Pinned by
`tests/cpu/semiempirical/test_dftb3_scc.cpp` (Euler identity, finite
differences, \(\Delta q\equiv 0\) reduction).

## Hamiltonian twin note

In-tree static tables `basis_data/dftb_mio_tables.h` for H/C/N/O only.
**Not** a full DFTB.org mio-1-1 `.skf` dump — hoppings / \(U^{\mathrm{d}}\)
are mio-style product-path sketches (no runtime parser). Absolute energy
identity with DFTB+ + production mio-1-1 is **deferred** to SK-fidelity
slugs (`dftb2-3ob-organic`, energy-pin). This regression pins the
**ATV DFTB3 Hamiltonian twin** energy at the stated geometry (deterministic
offline value measured under SE-only build), with formula lineage from
DFTB+/Gaus.

## Geometry

Water, slightly off-equilibrium (Å):

| Atom | x    | y    | z   |
|------|------|------|-----|
| O    | 0.00 | 0.00 | 0.0 |
| H    | 0.96 | 0.00 | 0.0 |
| H    |−0.24 | 0.93 | 0.0 |

Cubic box 20 Å (isolated).

## Parameters

- `scc_tol = 1e-10`, `max_scc = 50`, forces off (energy gate)
- Units: energy **eV**, positions **Å**

## Pin (provisional; DFTB+ oracle owed)

| Quantity | Value | Unit | Tol |
| -------- | ----- | ---- | --- |
| \(E_{\mathrm{tot}}\) | **−107.618793920072** | eV | \(10^{-6}\) |

Measured under Release build, `scc_tol=1e-10`, `max_scc=50`, forces off
(2026-07-28). Components at pin geometry (informational):
\(E^{(2)}\approx 22.537\) eV, \(E^{(3)}\approx 5.015\) eV
(sketch SK + \(U^{\mathrm{d}}\)).

**What this number is.** A *provisional self-comparison pin* of the ATV DFTB3
Hamiltonian twin — **not** an oracle value and **not** a converged
self-consistent energy (see the caveat below). Its only jobs are to freeze the
current third-order arithmetic and to fire the moment anything in the SCC
changes.

**What invalidated the previous pin (−84.561719022489 eV).** Two independent
defects, both now fixed in `src/cpu/semiempirical/dftb3.cpp`:

1. The third-order SCC potential implemented only one of the two contractions
   of \(\partial E^{(3)}/\partial\Delta q_C\). Because
   \(\Gamma_{AB}\neq\Gamma_{BA}\), the missing \(\tfrac13\sum_B\Delta
   q_B^2\Gamma_{BC}\) term left \(\sum_A\Delta q_A V^{(3)}_A\approx 2E^{(3)}\)
   instead of \(3E^{(3)}\), i.e. the potential was not the derivative of the
   reported energy (0.457 eV error on O at \(\Delta q=(-0.6,0.3,0.3)\)).
2. The band-structure double counting subtracted \(E^{(2)}+E^{(3)}\) where the
   Euler identity requires \(E^{(2)}+2E^{(3)}\).

The pin was additionally recorded under the branch's pre-upgrade mio H row;
that staleness is resolved by the `dev` merge, not here, and will move this
number again.

**Caveat — the SCC is not converged (known defect, tracked).** The undamped
Picard iteration in `DFTB3::launch` (and identically in `DFTB2::launch`) does
not converge for this system: it settles into a period-2 limit cycle with a
max\(|\Delta\Delta q|\) residual of ≈3.6 e, so `scc_tol` is never met and the
loop always exhausts `max_scc`. The pinned value is therefore the branch of
that cycle selected by an **even** `max_scc`; with `max_scc = 49` or `51` the
same build reports **−337.034306211295 eV**. Nothing in the third-order fix
causes this — the pre-fix code had the same cycle (−84.561719022489 eV even,
−339.982849250505 eV odd). Fixing it requires charge mixing in **both**
`dftb2.cpp` and `dftb3.cpp` and is out of scope for this slug.

**DFTB+ oracle owed.** No DFTB+ binary exists on this host (`module spider
dftb` → nothing), so no fresh oracle number could be generated. Replace this
pin with a DFTB+ 24.1 + mio-1-1 value once the offline oracle infrastructure
lands under `tools/tb_oracle/inputs/dftbplus/`.

## Reduction check (a real oracle, not a snapshot)

For a homonuclear system the Mulliken fluctuations vanish by symmetry
(\(\Delta q\equiv 0\)), so \(E^{(2)}=E^{(3)}=0\) identically and DFTB3 must
reduce **exactly** to DFTB2/DFTB0. Verified on H2 at \(r=0.74\) Å: DFTB3 and
DFTB2 agree bit-for-bit at **−18.137603579048 eV** (difference 0.0e+00), and
\(E^{(2)}=4.6\times10^{-31}\) eV, \(E^{(3)}=0\). That absolute value carries
this branch's pre-upgrade mio H row; the DFTB+ 24.1 + mio-1-1 pin
**−18.36635052 eV** (`regressions/semiempirical/dftb2-dftbplus-energy-pin/`)
becomes reachable after the `dev` merge. Gated by
`tests/cpu/semiempirical/test_dftb3.cpp::ReducesToDftb2ForHomonuclearH2`.

## Gates

- Finite, non-zero total energy [eV] on `frame.meta[kEnergy]`
- Hard-coded energy pin vs table above (tol \(10^{-6}\) eV)
- \(|E_{\mathrm{DFTB3}} - E_{\mathrm{DFTB2}}| > 10^{-6}\) eV (third-order active)

### Open: the DFTB3−DFTB2 magnitude is still not commensurate with \(E^{(3)}\)

At this geometry \(|E_{\mathrm{DFTB3}}-E_{\mathrm{DFTB2}}| = 69.07\) eV against
\(E^{(3)} = 5.01\) eV. The third-order fix does **not** cure this, and the
residual is not a third-order defect. Two second-order-level causes, both
shared with `dftb2.cpp` and therefore out of scope here:

1. **Neither method's SCC converges** (above), so the two energies are sampled
   from two *different* limit cycles at \(|\Delta q_{\mathrm O}|\approx 1.6\)–
   \(1.9\ e\); the comparison is not between two self-consistent states.
2. **The double-counting route itself.** Both methods assemble
   \(E = E_{\mathrm{band}} - E_{\mathrm{dc}} + E_{\mathrm{rep}}\), but the exact
   identity is \(E_{\mathrm{band}} = \mathrm{Tr}(\rho H^0) + \sum_A V_A q_A\)
   with the **gross** Mulliken population \(q_A\), not \(\Delta q_A\). The
   band-energy route therefore carries a spurious \(\sum_A V_A q^0_A\) term that
   is *linear* in \(\Delta q\); at a properly converged \(\Delta q\) for this
   case it is \(\approx 7.7\) eV (DFTB2) / \(7.5\) eV (DFTB3).

With both fixed (charge mixing + the \(\mathrm{Tr}(\rho H^0)+E^{(2)}+E^{(3)}\)
route) the same tables give \(|E_{\mathrm{DFTB3}}-E_{\mathrm{DFTB2}}| =
1.98\times10^{-3}\) eV against \(E^{(3)} = -2.04\times10^{-3}\) eV — i.e. the
smell test only clears when DFTB2 is repaired alongside DFTB3.

## Perf

`perf_exempt: true` — host has no DFTB+ binary for fair wall ratio; formula
oracle only. Revisit when offline DFTB+ pin infrastructure lands under
`tools/tb_oracle/inputs/dftbplus/`.
