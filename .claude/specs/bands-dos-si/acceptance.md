# Acceptance — `bands-dos-si`

Every item must be verified before `acceptance_met=true` / `status=green`.

## Correctness (dual gate part 1)

1. **Offline pin exists** at `regressions/pw/bands-dos-si/REFERENCE.md` with:
   - In-tree UPF: `tests/data/pw/Si_lda_nc.upf` (ONCVPSP LDA, same family as
     `pw-lda-si`).
   - Literature Si LDA **indirect gap** citation: Yin & Cohen, *Phys. Rev. B*
     **26**, 5668 (1982) (~0.45 eV at converged LDA); experimental ~1.12 eV.
   - Twin geometry \(a=5.43\) Å diamond Si (same as `pw-lda-si`).
   - Twin settings: SCF \(E_\mathrm{cut}=15\) Ry, \(k_\mathrm{mesh}=2\),
     XC=`LSDA`; bands: diamond path with `n_seg=4`, `n_bands=8`; DOS:
     `k_mesh=4`, \(\sigma=0.05\) Ry, `n_bands=8`.
   - Stated acceptance windows (incomplete basis) for gap and DOS integral —
     **not** a golden of Atomiverse-only free-hand numbers without external
     science context.
   - CP2K note: Quickstep `BAND_STRUCTURE` / DOS use GTH+Gaussian; absolute
     \(\varepsilon_n\) not comparable to NC-UPF PW; character + literature
     window is the pin. Optional oracle input under `tools/cp2k_oracle/inputs/`.
   - Perf note: `perf_exempt` or measured \(T_\mathrm{atv}/T_\mathrm{cp2k}\).

2. **Public API**:
   - Headers: `atomiverse/cpu/planewave/bands_dos.h` (or equivalent public
     path under `cpu/planewave/`).
   - Types `cpu::pw::Bands`, `cpu::pw::Dos`, `BandsResult`, `DosResult`.
   - Parameters are **members** — no public `*Config`, no `make_*`/`build_*`.
   - `Bands::set_diamond_path(n_seg)` fills fractional k-list Γ–X–W–L–Γ–K.
   - `compute(const PW&, const Frame&)` runs SCF (method members) then NSCF.

3. **Unit tests** (`tests/cpu/planewave/`, no full SCF e2e):
   - Diamond path has expected endpoint count / includes Γ and X.
   - Gaussian DOS of a single eigenvalue integrates to weight within \(10^{-3}\).
   - Synthetic spectrum: VBM/CBM/gap finder returns known gap.

4. **E2E regression** `regressions/pw/bands-dos-si/reg_bands-dos-si.cpp`:
   - Public headers only; links `atomiverse`.
   - Si diamond Frame; `PW` members (`upf_path`, `ecut_ry=15`, `k_mesh=2`,
     `xc=LSDA`, `tol`, `max_iter`) — no Config bag.
   - `Bands` with `n_bands=8`, `set_diamond_path(4)`; assert
     `scf_converged`, sorted \(\varepsilon_n(\mathbf{k})\),
     \(n_\mathrm{bands}=8\) columns; every **local** same-\(k\) gap \(> 0\).
   - Path samples \(\Gamma\) and \(X\); **indirect character**
     \(\varepsilon_c(X) < \varepsilon_c(\Gamma)\).
   - **Γ→X gap window** (incomplete ecut/k) — textbook Si transition
     \(E_{\Gamma\to X}=\varepsilon_c(X)-\varepsilon_v(\Gamma)\), not sparse-path
     \(\max\varepsilon_v-\min\varepsilon_c\):
     \[
     0.05~\mathrm{eV}
       \le E_{\Gamma\to X}^\mathrm{ATV}\times 13.605693~\mathrm{(if\ Ry)}
       \le 2.50~\mathrm{eV}.
     \]
     Literature: Yin & Cohen LDA ~0.45 eV. Measured twin ~0.48 eV.
   - Direct \(\Gamma\) gap \(\ge\) Γ→X − 0.05 eV.
   - `Dos` with `k_mesh=4`, `n_bands=8`, `sigma_ry=0.05`:
     occupied integral \(n_\mathrm{elec\_to\_vbm}\) within **1.0** of
     \(n_e=8\) (tutorial mesh/broadening; rejects unit disasters).
   - Exit non-zero with got/want on failure; ≤ 60 s single-rank CPU.

5. **CMake**: library source listed; unit test registered; regression
   `atv_add_regression(pw bands-dos-si)` under `ATV_BUILD_PW`.

6. **Build + run PASS**:
   ```bash
   cmake -B build -DCMAKE_BUILD_TYPE=Release \
     -DATV_BUILD_TESTS=ON -DATV_BUILD_HF=ON -DATV_BUILD_PW=ON \
     -DATV_BUILD_GPW=ON -DATV_BUILD_REGRESSIONS=ON \
     -DATV_BUILD_CUDA=OFF -DBLA_VENDOR=OpenBLAS
   cmake --build build --target reg_bands-dos-si test_cpu_pw_bands_dos -j$(nproc)
   ./build/tests/test_cpu_pw_bands_dos
   ./build/regressions/reg_bands-dos-si
   ```

## Performance (dual gate part 2)

7. If CP2K binary available offline: fair wall same geometry; record
   \(r=T_\mathrm{atv}/T_\mathrm{cp2k}\) in REFERENCE; require \(r\le 10\).
   Caveat: CP2K Quickstep GTH vs NC-UPF PW is wall-only.
8. If no CP2K binary: `perf_exempt=true` with reason in REFERENCE.md.

## Process

9. Spec written **before** production implementation on `wf/bands-dos-si`.
10. Commit(s) on `wf/bands-dos-si` only; no push; never on `dev`.
11. No public `*Config` bag and no `make_*`/`build_*` factory introduced.
