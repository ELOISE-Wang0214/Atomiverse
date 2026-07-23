# Acceptance — `virtual-site-multi-topology`

Binding checklist. `acceptance_met=true` only when every row is verified.

| ID | Criterion | How verified | Pass |
| -- | --------- | ------------ | ---- |
| A1 | Oracle LAMMPS tip4p / virtual-site path | REFERENCE.md `oracle: lammps` | [x] |
| A2 | Offline pin under `tools/metalwalls_oracle/inputs/lammps/` | files present | [x] |
| A3 | `oracle:` + `oracle_version:` in REFERENCE | grep | [x] |
| A4 | TIP4P dimer E vs pin ≤ 1e-8 eV | `reg_virtual-site-multi-topology` PASS | [x] |
| A5 | OM = qdist ≤ 1e-10 Å | regression dual-gate PASS | [x] |
| A6 | O/H forces ≤ 1e-6 eV/Å; M force ~0 after redistribute | regression PASS | [x] |
| A7 | TIP5P r_OL / θ_LOL geometry dual-gate | regression + unit PASS | [x] |
| A8 | Force partition sum identity ≤ 1e-12 | unit test PASS | [x] |
| A9 | MD-only build flags | HF/PW/GPW/CUDA/Torch OFF | [x] |
| A10 | tests/ unit only; one e2e regression | layout | [x] |
| A11 | Perf ≤10× or perf_exempt documented | ratio 0.0837 | [x] |
