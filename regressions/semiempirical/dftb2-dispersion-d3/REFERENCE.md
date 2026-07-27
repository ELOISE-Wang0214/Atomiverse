# Reference — DFTB2 + D3(BJ) water dimer dispersion

## Method

SCC-DFTB2 electronic + repulsive (ATV static mio-style SK sketch) with
optional **Grimme DFT-D3 Becke–Johnson** two-body correction (`DFTB2::d3`).

Damping parameters = published **DFTB3-D3(BJ)** set (DFTB+ / Grimme list;
`s-dftd3 --bj dftb3`):

| param | value |
| ----- | ----- |
| \(s_6\) | 1.0 |
| \(s_8\) | 0.5883 |
| \(a_1\) | 0.5719 |
| \(a_2\) | 3.6017 Bohr |

CN form: molecular Grimme (simple-dftd3 / PySCF path), **not** CP2K tanh.
ATM off. Non-periodic molecular evaluation.

## Oracle

```text
oracle: literature
oracle_version: Simple-DFTD3/s-dftd3 1.2.1 (offline)
```

Primary numerical pin is **pure \(E_{\mathrm{disp}}\)** from offline
**s-dftd3 1.2.1** (Simple-DFTD3). DFTB+ is the family oracle for the electronic
DFTB2 ladder; for this geometry-only D3 correction the literature/simple-dftd3
engine is the correct pin (same model DFTB+ uses for D3). ctest never launches
s-dftd3 or DFTB+.

## Hard gate — E_disp

**E_disp = −0.0016892886213669401 Ha = −0.04596788515696838 eV**

Source: offline s-dftd3 1.2.1 on the geometry below:

```text
s-dftd3 --bj dftb3 --json pin.json dimer.xyz
# Dispersion energy:      -1.6892886213669E-03 Eh
```

Conversion: CODATA / ATV `kHartreeInEv = 27.211386245988`.

Absolute total DFTB2 energy is **not** pinned to DFTB+ (ATV SK is a mio-style
sketch). Regression gates:

1. `|E_disp_ATV − E_disp_pin| ≤ 1e-4 eV`
2. `E_tot(d3) − E_tot(off) = E_disp` within 1e-8 eV

## Geometry (Å)

| Atom | x | y | z |
| ---- | - | - | - |
| O | 0.000 | 0.000 | 0.000 |
| H | 0.757 | 0.000 | 0.586 |
| H | −0.757 | 0.000 | 0.586 |
| O | 0.000 | 0.000 | 2.900 |
| H | 0.757 | 0.000 | 3.486 |
| H | −0.757 | 0.000 | 3.486 |

Cubic box 20 Å (isolated; D3 path ignores PBC).

## Offline reproduce

```bash
export PATH=/software/sse2/el9_epyc9005/easybuild/prefix/software/Simple-DFTD3/1.2.1-gfbf-2025b/bin:$PATH
# twin: tools/tb_oracle/inputs/dftbplus/dftb2-dispersion-d3/dimer.xyz
s-dftd3 --bj dftb3 --json pin.json dimer.xyz
```

Oracle inputs: `tools/tb_oracle/inputs/dftbplus/dftb2-dispersion-d3/`.

## Perf

`perf_exempt: true` — pure D3 pin is microseconds; fair \(T_{\mathrm{atv}}/T_{\mathrm{oracle}}\)
vs s-dftd3 on a 6-atom cluster is not a product gate. Full DFTB2 SCC wall is
electronics-dominated.
