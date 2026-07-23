# Acceptance — `pw-lda-si`

Every item must be verified before `acceptance_met=true` / `status=green`.

## Correctness (dual gate part 1)

1. **Offline pin exists** at `regressions/pw/pw-lda-si/REFERENCE.md` with:
   - ABINIT-class same-pseudo pin: PseudoDojo `ONCVPSP-LDA-SR-noCCv0.4`
     `Si/Si.djrepo` deltafactor \(E_\mathrm{min} = -230.158202\) eV/cell
     \(= -16.916316\) Ry at \(a_0 \approx 5.468\) Å (num_sites=2 diamond),
     retrieved from
     <https://github.com/PseudoDojo/ONCVPSP-LDA-SR-noCCv0.4> (`Si/Si.djrepo`,
     ecut-hint key `8.0`). Cite Hamann PRB 88, 085117 (2013); van Setten
     et al. CPC 226, 39 (2018).
   - In-tree UPF: `tests/data/pw/Si_lda_nc.upf` (same ONCVPSP LDA noCC family).
   - Twin settings for the ≤ 60 s regression: \(a=5.43\) Å, \(E_\mathrm{cut}=15\) Ry,
     \(k_\mathrm{mesh}=2\), SCF tol \(10^{-6}\).
   - Stated absolute tolerance and incomplete-basis rationale (not a golden of
     Atomiverse-only output).
   - Perf note: `perf_exempt` or measured \(T_\mathrm{atv}/T_\mathrm{cp2k}\).

2. **Geometry** (Å), 2-atom fcc primitive diamond Si:

   | Atom | x | y | z |
   | ---- | --: | --: | --: |
   | Si | 0 | 0 | 0 |
   | Si | \(a/4\) | \(a/4\) | \(a/4\) |

   Cell vectors: \((a/2)(0,1,1)\), \((a/2)(1,0,1)\), \((a/2)(1,1,0)\) with
   \(a = 5.43\) Å.

3. **Public API regression** `reg_pw-lda-si`:
   - Includes only public `atomiverse/…` headers; links `atomiverse`.
   - Constructs `cpu::pw::PW` with **members** (`upf_path`, `ecut_ry`, `k_mesh`,
     `tol`, `max_iter`) — no public `*Config` bag, no factories.
   - Runs via `cpu::pw::PWEvaluator` + `cpu::single_point(frame, pwe)`.
   - Asserts `last_result().converged == true`.
   - Asserts total energy [Ry] lies in the incomplete-basis band above the
     ABINIT-class \(E_\mathrm{min}\):
     \[
     E_\mathrm{min}^\mathrm{ABINIT} \le E_\mathrm{ATV}
       \le E_\mathrm{min}^\mathrm{ABINIT} + \Delta
     \]
     with \(E_\mathrm{min}^\mathrm{ABINIT} = -16.916316\) Ry and
     \(\Delta = 1.50\) Ry (tutorial ecut/k incomplete vs PseudoDojo-converged
     limit; catches unit disasters and unphysical energies while allowing
     known incomplete-basis elevation).
   - Asserts Frame `meta.get(kEnergy)` [eV] matches
     `last_result().energy * 13.605693` within \(10^{-4}\) eV.
   - Exit non-zero with `got / want / tol` (or band bounds) on failure;
     runs ≤ 60 s single-rank CPU.

4. **CMake**: case registered under `ATV_BUILD_PW` as
   `atv_add_regression(pw pw-lda-si)` with UPF path compile definition pointing
   at `${CMAKE_SOURCE_DIR}/tests/data/pw/Si_lda_nc.upf`.

5. **Build + run PASS** with:
   ```bash
   cmake -B build -DCMAKE_BUILD_TYPE=Release \
     -DATV_BUILD_TESTS=ON -DATV_BUILD_HF=ON -DATV_BUILD_PW=ON \
     -DATV_BUILD_GPW=ON -DATV_BUILD_REGRESSIONS=ON \
     -DATV_BUILD_CUDA=OFF -DBLA_VENDOR=OpenBLAS
   cmake --build build --target reg_pw-lda-si -j$(nproc)
   ./build/regressions/reg_pw-lda-si
   ```

## Performance (dual gate part 2)

6. If a CP2K binary is available offline: record
   \(r = T_\mathrm{atv}/T_\mathrm{cp2k}\) in REFERENCE.md; require \(r \le 10\).
   Note fair-bench caveat: CP2K Quickstep uses GTH-PADE (different PP); wall
   comparison is same-geometry wall only, not absolute-energy parity.
7. If no CP2K binary: set `perf_exempt=true` in REFERENCE.md with reason.
   Green is allowed under exempt.

## Process

8. Spec written **before** production regression on `wf/pw-lda-si`.
9. Commit(s) on `wf/pw-lda-si` only; no push; never commit on `dev`.
10. No public `*Config` bag and no `make_*`/`build_*` factory introduced.
