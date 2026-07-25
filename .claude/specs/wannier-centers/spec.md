# Spec: `wannier-centers` — Foster–Boys MLWF foothold

| Field | Value |
| ----- | ----- |
| Slug | `wannier-centers` |
| Branch | `wf/wannier-centers` |
| Queue | `docs/developer/cp2k-method-queue.md` row `wannier-centers` |
| Roadmap | `docs/roadmap/dft.md` §11 Maximally localized Wannier functions (○ → foothold) |
| Dep | `bands-dos-si` ✓ (k-space/response product path open; this foothold is molecular) |
| Dual gate | Correctness below; speed **`perf_exempt`** (no CP2K binary on pin host) |

## Problem

Roadmap §11 lists **maximally localized Wannier functions (MLWF)** as ○.
CP2K can oracle offline when available; the product needs a **shipped foothold**
that localizes occupied mean-field MOs and returns orbital centers [Å] and
spreads [Å²].

Periodic multi-k / GPW-grid MLWF (Marzari–Vanderbilt) is **out of scope** for
this foothold. Molecular Γ-point Foster–Boys on a closed-shell Gaussian RHF
wavefunction reuses existing multipole one-electron integrals and is the
fastest textbook path (H₂ midpoint; H₂O “rabbit ears”).

## Public API

Concrete class (no `*Config` bag, no factory):

```cpp
// include/atomiverse/cpu/scf/wannier.h
namespace ATV::cpu::scf {

struct WannierResult {
    std::vector<std::array<double, 3>> centers; ///< orbital centroids [Å]
    std::vector<double> spreads;                ///< Ω_i [Å²]
    double total_spread = 0;                    ///< Σ Ω_i after localization [Å²]
    double total_spread_initial = 0;            ///< Σ Ω_i of input occupied MOs [Å²]
    int sweeps = 0;                             ///< Jacobi sweeps executed
    bool converged = false;
};

/// Foster–Boys localization of occupied closed-shell MOs (Jacobi).
class WannierCenters {
  public:
    int max_sweeps = 200;   ///< max full pair-sweeps (dimensionless)
    double tol = 1e-10;     ///< stop when max |γ| over pairs < tol [rad]

    /// @param outs  MeanFieldOutputs after a restricted method `run()`
    /// @param frame Geometry [Å] (origin of multipole integrals = (0,0,0))
    WannierResult compute(const MeanFieldOutputs& outs,
                          const cpu::Frame& frame) const;
};

/// Pure math (unit-testable): ¼ atan2(B, A) for one Boys pair [rad].
double boys_pair_angle(double A, double B);

} // namespace
```

Call path:

```cpp
cpu::scf::RHF rhf;
rhf.basis = BasisSet::STO3G;
rhf.run(frame);
cpu::scf::WannierCenters wb;
auto r = wb.compute(rhf.result_handle(), frame);
// r.centers[i], r.spreads[i]
```

## Algorithm

**Foster–Boys** (Foster & Boys, *Rev. Mod. Phys.* **32**, 300 (1960)):
maximize \(\sum_i |\langle\mathbf r\rangle_i|^2\) among occupied spatial MOs
by successive 2×2 Jacobi rotations (equivalently minimize total second-moment
spread \(\Omega=\sum_i(\langle r^2\rangle_i-|\langle\mathbf r\rangle_i|^2)\)).

AO multipole matrices from `scf::MultipoleIntegrals` about origin \(\mathbf C=0\)
(Bohr). Transform to occupied MO basis
\(R_{ij}^a = \mathbf c_i^\top \langle r_a\rangle \mathbf c_j\).
Pair \((i,j)\):

\[
A_{ij}^a = R_{ii}^a - R_{jj}^a,\quad
B_{ij}^a = 2 R_{ij}^a,
\]
\[
A=\sum_a\bigl((A^a)^2-(B^a)^2\bigr),\quad
B=\sum_a 2 A^a B^a,\quad
\gamma=\tfrac14\operatorname{atan2}(B,A).
\]

Rotate columns of \(C_\mathrm{occ}\) by \(\gamma\). Sweep until
\(\max|\gamma|<\texttt{tol}\) or `max_sweeps`.

**Scope (foothold):** closed-shell restricted (RHF/RKS) only; integer
occupations; throw on open-shell / generalized / not-ready. Periodic
Marzari–Vanderbilt / GPW-grid MLWF deferred.

## Tests

| Tier | Path | Pins |
| ---- | ---- | ---- |
| Unit | `tests/cpu/scf/test_wannier.cpp` | `boys_pair_angle` quadrants; 1-orbital identity; synthetic 2-orbital rotation lowers spread; no SCF e2e |
| E2E | `regressions/scf/wannier-centers/` | H₂O / STO-3G RHF: 5 centers; \(\Omega_\mathrm{tot}\) decreases; classic Boys geometry (core + 2 LP + 2 OH); H₂ midpoint algebra pin in unit or same binary |

## Docs

- Tutorial: `docs/tutorial/methods/wannier-centers.md`
- Roadmap §11 row MLWF: foothold note (CPU ✓ foothold, not full periodic)
- Queue row → done on merge

## Non-goals

- Pipek–Mezey, Edmiston–Ruedenberg
- Periodic / multi-k Marzari–Vanderbilt
- GPW/PW frozen-MO path (same math later; multipole AO path is Gaussian first)
- Writing centers onto Frame meta (caller owns the `WannierResult`)
- CP2K binary wall (perf_exempt)
