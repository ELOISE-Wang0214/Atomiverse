# Acceptance — slab-2d-ewald-foothold

`acceptance_met=true` only when **every** item below is verified.

| ID | Criterion | How verified |
|---|---|---|
| A1 | Oracle is **LAMMPS** (`kspace_style ewald` + `kspace_modify slab`) | `REFERENCE.md` has `oracle: lammps` and `oracle_version:` |
| A2 | Offline pin exists (no live LAMMPS in ctest) | Pin path + hard-coded numbers in REFERENCE / reg body |
| A3 | MD-only build flags | configure with HF/PW/GPW/CUDA/Torch **OFF**, MD+KSPACE+TESTS+REGRESSIONS **ON** |
| A4 | Observable: total electrostatic energy \(E\) (eV) of parallel-plate slab toy | regression prints \(E\) and asserts vs pin |
| A5 | Energy tol: \(\lvert E_{\mathrm{atv}} - E_{\mathrm{pin}}\rvert / \lvert E_{\mathrm{pin}}\rvert \le 1\times 10^{-4}\) (relative), or absolute \(5\times 10^{-3}\) eV | reg exit 0 |
| A6 | \(F_z\) on electrode-like planes: sign + magnitude vs pin within \(1\times 10^{-3}\) eV/Å | reg asserts mean \(F_z\) per plane |
| A7 | Dual-gate: correctness PASS + (within_10x OR perf_exempt with reason) | REFERENCE perf section |
| A8 | No public `*Config` bag; no `make_*` / `build_*` factory for this feature | API is `KSpaceEwald` ctor + `slab_factor` member |
| A9 | Architecture: Component on existing MD/Evaluate stack; no new Driver category; no charge solve in KSpace | code review of diff |

## Pin geometry (binding)

- Box: \(12\times 12\times 20\) Å, orthorhombic  
- 4× (+1e) at \(z=5\) Å on a 6 Å xy grid; 4× (−1e) at \(z=15\) Å  
- `coul/long` cutoff 5.9 Å; `gewald` = α = 0.5 Å⁻¹; accuracy 1e-6  
- `kspace_modify slab 3.0`  
- Units: metal (eV, Å, e)  

Pinned \(E\) and \(F_z\) recorded in `regressions/md/slab-2d-ewald-foothold/REFERENCE.md`.
