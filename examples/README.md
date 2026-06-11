# Examples

A flat, numbered tutorial sequence. Each file teaches exactly one concept and
is meant to be read top-to-bottom; numbers define the recommended reading
order.

| # | File | Concept | Requires |
|---|------|---------|----------|
| 01 | `01_scf.cpp` | The bare SCF loop on raw xtensor arrays (no atomiverse link) | — |
| 02 | `02_rhf.cpp` / `02_rhf.cu` | First real chemistry: RHF single-point via `Evaluate + scf::SCFEvaluator` (CUDA twin shows the same lesson on device) | `ATV_BUILD_RHF` (+`ATV_BUILD_CUDA` for the twin) |
| 03 | `03_uhf.cpp` | Open-shell UHF | `ATV_BUILD_RHF` |
| 04 | `04_uks.cpp` | DFT (UKS) with the native XC stack | `ATV_BUILD_RHF` |
| 05 | `05_geom_opt.cpp` | Geometry optimization: `Minimize + LBFGS + SCFEvaluator` | `ATV_BUILD_RHF` |
| 06 | `06_aimd.cpp` | Ab-initio MD: `MD + VelocityVerlet + SCFEvaluator` | `ATV_BUILD_RHF` |
| 07 | `07_md_lj.cu` | Classical CUDA MD: `cuda::md::MD + NeighborList + PairLJCut` | `ATV_BUILD_CUDA` |
| 08 | `08_rhf_multi_gpu.cu` | Multi-GPU RHF | `ATV_BUILD_RHF` + `ATV_BUILD_CUDA` |

## Rules for adding an example

1. **One concept per file.** If a file needs two headlines, split it.
2. **Numbered, flat, ordered.** Pick the next free number in the pedagogical
   order; no subdirectories.
3. **CPU `.cpp` is the canonical lesson.** A CUDA `.cu` twin reuses the same
   number and gets a `_cuda`-suffixed target (`02_rhf.cu` → `02_rhf_cuda`)
   unless the concept itself is CUDA-only (07, 08).
4. **Examples are not benchmarks.** Anything that measures speed or compares
   against LAMMPS / PySCF / QE belongs in `benchmarks/<domain>/<case>/`.
5. Register the target in `CMakeLists.txt` behind the narrowest option gate
   (`ATV_BUILD_RHF`, `ATV_BUILD_CUDA`, …).

Build: configure with `-DATV_BUILD_EXAMPLES=ON`, then
`cmake --build build --target 02_rhf` (any number works the same way).
