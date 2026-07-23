---
slug: bands-dos-si
status: open
created: 2026-07-23
track: cp2k-method-loop
domain: pw
grilled: false
---

# bands-dos-si

| Field | Value |
| ----- | ----- |
| Slug | `bands-dos-si` |
| Branch | `wf/bands-dos-si` |
| Domain | plane-wave DFT post-processing (band structure + DOS) |
| Observable | \(\varepsilon_n(\mathbf{k})\) [Ry], indirect gap [eV], DOS [1/Ry] |
| Geometry | Si diamond, 2-atom fcc primitive cell, \(a = 5.43\) Å |
| Basis / PP | NC-UPF ONCVPSP LDA `tests/data/pw/Si_lda_nc.upf` |
| XC | LDA (`XCType::LSDA`) |
| Settings (twin) | SCF: \(E_\mathrm{cut}=15\) Ry, \(k_\mathrm{mesh}=2\); bands: diamond path Γ–X–W–L–Γ–K with 4 pts/segment, 8 bands; DOS: \(k_\mathrm{mesh}=4\), \(\sigma=0.05\) Ry |
| Public API | `cpu::pw::Bands` / `cpu::pw::Dos` (members + `compute(PW, Frame)`) — no public `*Config` |
| Regression | `regressions/pw/bands-dos-si/` |
| Queue row | `docs/developer/cp2k-method-queue.md` `## in_progress` `[~] bands-dos-si` |

## Goal

Ship a **band structure + density-of-states post-processing foothold** on
CPU plane-wave DFT: SCF → freeze \(V_\mathrm{eff}\) → non-self-consistent
diagonalization on a high-symmetry **k-path** (bands) and on a denser
Monkhorst–Pack mesh (DOS). Product API is OOP members on `Bands` / `Dos`;
Frame is geometry I/O only. Enables later Wannier / Berry / GW response
slugs that need \(\varepsilon_n(\mathbf{k})\).

## Domain basis

- KS eigenvalues from the frozen self-consistent effective potential
  (standard “bands” / NSCF step in QE, ABINIT, CP2K `BAND_STRUCTURE`).
- Silicon is an **indirect** semiconductor: VBM at \(\Gamma\), CBM near \(X\)
  (Yin & Cohen, *Phys. Rev. B* **26**, 5668 (1982); experimental gap
  \(\approx 1.12\) eV; LDA underestimates \(\sim 0.45\) eV at converged
  settings).
- DOS: Gaussian-broadened sum over k-weights and bands
  \[
  D(E) = \sum_{n\mathbf{k}} w_{\mathbf{k}}\, g_\sigma(E-\varepsilon_{n\mathbf{k}})
  \]
  with spinless closed-shell convention \(f\in\{0,2\}\) already absorbed into
  band counting (two electrons per band when integrating occupied DOS).
- Units: eigenvalues reported in **Ry** (QE-comparable) on the result
  structs; gaps also available in eV via `qe_units::rydberg_to_ev`.

## Non-goals

- GPW / Gaussian k-path (slug `gpw-kpoints-si`).
- Spin-polarized bands / spin-resolved DOS.
- Berry phase, Wannier centers, DFPT phonons (later queue rows).
- Public `*Config` bags, `make_*` / `build_*` factories.
- Live CP2K / QE / ABINIT at ctest time.
- MPI multi-rank scaling.
- Production-converged ecut (~60 Ry) if it breaks the ≤ 60 s budget.

## Design

### Public API

```cpp
namespace ATV::cpu::pw {

struct BandsResult {
    std::vector<std::array<double, 3>> kfrac;           // fractional k
    std::vector<std::vector<double>> eigenvalues_ry;    // [ik][ib]
    double e_vbm_ry = 0, e_cbm_ry = 0, gap_ry = 0;
    int vbm_k = 0, cbm_k = 0, vbm_band = 0, cbm_band = 0;
    bool scf_converged = false;
    double energy_ry = 0;  // SCF total energy [Ry]
};

struct DosResult {
    std::vector<double> energy_ry;  // bin centres [Ry]
    std::vector<double> dos;        // states / Ry / cell
    double n_elec_to_vbm = 0;       // ∫_{-∞}^{VBM} D(E) dE (×2 spinless)
    double e_vbm_ry = 0, e_cbm_ry = 0, gap_ry = 0;
    bool scf_converged = false;
};

/// High-symmetry path band structure (NSCF on frozen V_eff).
class Bands {
  public:
    int n_bands = 8;  // valence + conduction (Si: 4 occupied valence)
    std::vector<std::array<double, 3>> kfrac; // empty → set_diamond_path
    void set_diamond_path(int n_seg = 4);
    BandsResult compute(const PW& method, const cpu::Frame& frame) const;
};

/// Gaussian-broadened DOS from denser MP mesh eigenvalues.
class Dos {
  public:
    int k_mesh = 4;
    int n_bands = 8;
    double sigma_ry = 0.05;
    int n_bins = 200;
    double e_min_ry = 0; // if e_max <= e_min → auto window from evals
    double e_max_ry = 0;
    DosResult compute(const PW& method, const cpu::Frame& frame) const;
};

} // namespace
```

Parameters are **members** of `Bands` / `Dos`. No `BandsConfig`. Construction
is default; caller sets members then calls `compute(method, frame)`.

### Internal

- Reuse `setup_pw_problem` + restricted fixed-occupation SCF loop (same
  shape as forces’ private `run_scf`).
- `rebuild_k_tables(PWProblem&, …, kfrac)` replaces only k-spheres and
  KB/atom-wfc tables; density and \(V_\mathrm{eff}\) stay frozen.
- Davidson with `n_bands` on each path / DOS-mesh k-point.
- Gap: \(n_\mathrm{occ} = n_e/2\); VBM = max \(\varepsilon_{n_\mathrm{occ}-1}(k)\);
  CBM = min \(\varepsilon_{n_\mathrm{occ}}(k)\).

## Tasks

1. Write this spec + acceptance on `wf/bands-dos-si` before production code.
2. Offline pin in `regressions/pw/bands-dos-si/REFERENCE.md` (literature LDA
   Si gap + twin settings + DOS integral rule; CP2K BAND_STRUCTURE note).
3. Implement `rebuild_k_tables`, `Bands`, `Dos` (public headers + source).
4. Unit tests under `tests/cpu/planewave/` (path geometry, DOS Gaussian
   integral, synthetic gap finder — **no** full e2e SCF in unit tests).
5. E2E regression `regressions/pw/bands-dos-si/`.
6. Wire CMake (library source, unit test, regression under `ATV_BUILD_PW`).
7. Build + run until PASS; perf gate or `perf_exempt`.
8. Commit on `wf/bands-dos-si` only.

## Acceptance criteria

See `acceptance.md` (binding).
