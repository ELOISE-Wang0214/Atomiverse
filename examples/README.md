# Examples

A flat, numbered tutorial sequence. Each file teaches exactly one concept and
is meant to be read top-to-bottom; numbers define the recommended reading
order.

| # | File | Concept | Requires |
|---|------|---------|----------|
| 01 | `01_scf.cpp` | The bare SCF loop on raw xtensor arrays (no atomiverse link) | — |
| 02 | `02_rhf.cpp` / `02_rhf.cu` | First real chemistry: RHF single-point (CPU via `scf::rhf_solve`; CUDA twin via `single_point + cuda::scf::SCFEvaluator`) | `ATV_BUILD_RHF` (+`ATV_BUILD_CUDA` for the twin) |
| 03 | `03_uhf.cpp` | Open-shell UHF (`scf::uhf_solve`) + spin contamination | `ATV_BUILD_RHF` |
| 04 | `04_uks.cpp` | DFT (UKS-LSDA) with the native XC stack (`scf::uks_solve`) | `ATV_BUILD_RHF` |
| 05 | `05_geom_opt.cpp` | Geometry optimization: `Minimize + LBFGS + scf::SCFEvaluator` (DFT/PBE forces) | `ATV_BUILD_RHF` |
| 06 | `06_aimd.cpp` | Ab-initio MD: `MD + VelocityVerlet + NoseHoover + scf::SCFEvaluator` | `ATV_BUILD_RHF` |
| 07 | `07_md_lj.cpp` / `07_md_lj.cu` | Classical MD (NVE LJ argon): `md::MD + NeighborList + PairLJCut` (CPU canonical; CUDA twin = `07_md_lj_cuda`) | — (`.cu` needs `ATV_BUILD_CUDA`) |
| 08 | `08_pw.cpp` | Plane-wave DFT for a periodic crystal: `Evaluate + pw::PWEvaluator` (silicon) | `ATV_BUILD_PW` |
| 09 | `09_gpw.cpp` | Gaussian-and-plane-waves DFT single point: `gpw::GPWSolver` | `ATV_BUILD_GPW` |
| 10 | `10_gpw_znoedl.cpp` | GPW capstone: 1214-atom ZnO electric-double-layer system | `ATV_BUILD_GPW` |
| 11 | `11_koopmans.cpp` | Koopmans' theorem: HOMO energy → vertical IP via `scf::OrbitalInfoObserver` (5×2 basis sweep) | `ATV_BUILD_RHF` |
| 12 | `12_delta_scf.cpp` | ΔSCF ionization energies beyond Koopmans: neutral RHF vs `SCFInitialGuess`-seeded cation UHF | `ATV_BUILD_RHF` |
| 13 | `13_mp2.cpp` | Closed-shell MP2 correlation energy (`scf::mp2_energy`, AO→MO ovov transform) | `ATV_BUILD_RHF` |
| 14 | `14_fci.cpp` | Exact correlation via tiny full CI (`scf::fci_energy`) on H2/STO-3G | `ATV_BUILD_RHF` |

## Rules for adding an example

1. **One concept per file.** If a file needs two headlines, split it.
2. **Numbered, flat, ordered.** Pick the next free number in the pedagogical
   order; no subdirectories.
3. **CPU `.cpp` is the canonical lesson.** A CUDA `.cu` twin reuses the same
   number and gets a `_cuda`-suffixed target (`02_rhf.cu` → `02_rhf_cuda`,
   `07_md_lj.cu` → `07_md_lj_cuda`).
4. **Examples are not benchmarks.** Anything that measures speed or compares
   against LAMMPS / PySCF / QE belongs in `benchmarks/<domain>/<case>/`.
5. Register the target in `CMakeLists.txt` behind the narrowest option gate
   (`ATV_BUILD_RHF`, `ATV_BUILD_CUDA`, …).

Build: configure with `-DATV_BUILD_EXAMPLES=ON`, then
`cmake --build build --target 02_rhf` (any number works the same way).
