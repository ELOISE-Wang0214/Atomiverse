# Reference — Monkhorst-Pack k-point GPW RKS on diamond Si (bct 4-atom cell)

## Method

Closed-shell GPW Kohn–Sham DFT (PBE) with **Monkhorst–Pack k-point
sampling** — the CP2K Quickstep `&KPOINTS SCHEME MONKHORST-PACK` loop shape:
complex Bloch-phased analytic matrices per k
(S_k / T_k lattice sums, phased Kleinman–Bylander GTH nonlocal), per-k
pre-summed complex basis-on-grid values for the k-summed density
ρ(r) = Σ_k w_k Σ_μν P^k_μν χ^k_μ conj(χ^k_ν), per-k complex Hermitian
canonical-orthogonalization + zheevd diagonalization, and k-summed Pulay
DIIS. CP2K MP convention κ_d = (2r − N − 1)/(2N) (N = 2 → ±1/4, Γ-avoiding);
time-reversal folding of ±k (exact for RKS).

Public API: concrete `cpu::gpw::GpwDft` with members
(`k_mesh`, `k_time_reversal`, explicit `k_points`/`k_point_weights`,
`cutoff_ry`, …); lookup constructors
`GaussianBasis(Element::Si, "DZVP-MOLOPT-GTH")` /
`GthPotential(Element::Si, "GTH-PBE-q4")`. No public `*Config`, no factories.
Result payload: `GpwResult::{k_frac, k_weights, eps_k}`.

## System

Silicon diamond in the **body-centred-tetragonal 4-atom cell**
(the orthorhombic-compatible setting of the fcc diamond lattice):

| | |
| --- | --- |
| Lattice | a = 5.43 Å (10.26120 Bohr); cell = (a/√2, a/√2, a) |
| Atoms (fractional) | (0,0,0), (½,½,½), (½,0,¼), (0,½,¾) |
| Basis | DZVP-MOLOPT-GTH (Si) — in-tree table (`gpw_tables.h`) |
| Pseudopotential | GTH-PBE-q4 |
| Valence | 16 e⁻ → n_occ = 8 |
| XC | PBE |
| Density cutoff | 150 Ry (main run); 95 Ry (identity gates — commensurate supercell grids, unit-test-pinned) |
| k mesh | MONKHORST-PACK 2×2×2, time-reversal folded → 4 k, w = ¼ |
| SCF | DIIS depth 6, tol 1e-8 Ha, EPS_DEFAULT 1e-7 |

## Offline oracle status

**No CP2K binary / module on the worker host at pin time** (2026-07-25:
`module avail cp2k` empty; no `cp2k`, `cp2k.psmp`, `cp2k.sopt`, `cp2k.ssmp`
on PATH) — same status as the `gpw-pbesol-si` / `gpw-d4-water` pins.
Absolute CP2K digits are therefore **not invented** and **no live process**
runs at regression time. Pin strategy (binding, mirrors precedent):

1. **Internal exactness identities** — the decisive k-point correctness
   checks; each holds to SCF tolerance *independent of any oracle*:
   - **Supercell folding** (gate A2): E(bct, mesh 1×1×2 ⇒ κ_z = ±¼,
     genuinely complex phases {±1, ±i}) ≡ E(z-doubled 8-atom cell at the
     folded zone-boundary point κ = (0,0,½))/2. Measured **|Δ| = 6.9e-14 Ha**
     (tol 5e-5).
   - **Time reversal** (A3): full 8-k mesh vs TR-folded 4-k mesh —
     **|Δ| = 0.0** (tol 5e-6).
   - **Γ limit** (A1): explicit k = (0,0,0) through the complex machinery vs
     the pre-existing real Γ path (compact-basis periodic H₂ toy,
     DZV-GTH-PADE) — **|Δ| = 1.4e-11 Ha** (tol 1e-6).
   The A2 identity found a real bug at first run: the phased KB nonlocal
   assembly conjugated the wrong side (V_nl evaluated at −k; Hermitian,
   Γ-exact, TR-invisible — but 2.09 Ha off at κ = ¼). Fixed in
   `internal/kpoints.cpp`; the identity is now machine-exact.
2. **Frozen first-green ATV absolute band** (A7; setup integrity, not a
   third-party oracle). Replace with a cited CP2K number by running
   `tools/cp2k_oracle/inputs/gpw-kpoints-si.inp` when an offline CP2K
   exists, then tighten the tolerance.

### Atomiverse numbers (first green, 2026-07-25)

| Quantity | Value |
| --- | --- |
| **E(MP 2×2×2 TR, 150 Ry)** | **−15.8023941333 Ha** (SCF 6 iters) |
| E(MP 2×2×2 TR, 95 Ry) | −15.8024022460 Ha |
| E(explicit Γ, 150 Ry) | −15.2031525100 Ha |
| k-sampling shift \|E_MP − E_Γ\| | 5.99e-01 Ha (gate floor 1e-3) |
| E(mesh 1×1×2, 95 Ry) | −15.3483059090 Ha |
| E(8-atom supercell, κ_z = ½, 95 Ry)/2 | −15.3483059090 Ha |
| Sampled-k gap (VBM 0.137858 Ha, CBM 0.283137 Ha) | 3.953 eV (floor 0.1 eV) |

## Tolerances

- A1 ≤ 1e-6 Ha; A2 ≤ 5e-5 Ha; A3 ≤ 5e-6 Ha (identities).
- A4 ≥ 1e-3 Ha activation floor; A5 gap ≥ 0.1 eV.
- A7 absolute: 1e-3 Ha vs the frozen first-green value until an offline
  CP2K twin replaces it.

## What is checked

1. Γ-limit, supercell-folding, and time-reversal identities (A1–A3).
2. k-sampling activation and insulating gap + k payload shapes (A4–A5).
3. SCF convergence, finiteness (A6); frozen absolute band (A7).
4. No external process spawn; ≤ 60 s.

## Performance

| Item | Value |
| ---- | ----- |
| `perf_exempt` | **true** |
| reason | No CP2K binary/module on the pin host for a fair wall twin. |
| r = T_atv/T_cp2k | n/a (exempt) |
| `within_10x` | n/a (exempt) |
| T_atv (whole regression: 7 SCF solves, OMP=8, Release) | ~17 s |
| T_atv (main MP 2×2×2 run alone) | ~4 s |

When CP2K is available, run the twin input single-rank (same OMP), record
both walls and replace this section with measured r (target ≤ 10).

## Provenance

- Spec: `.claude/specs/gpw-kpoints-si/` (deleted on completion; git history)
- MP mesh: Monkhorst & Pack, PRB **13**, 5188 (1976); CP2K KPOINTS convention
- GTH: Goedecker–Teter–Hutter, PRB **54**, 1703 (1996); Krack, TCA **114**,
  145 (2005); MOLOPT: VandeVondele & Hutter, JCP **127**, 114105 (2007)
- Twin input: `tools/cp2k_oracle/inputs/gpw-kpoints-si.inp`
