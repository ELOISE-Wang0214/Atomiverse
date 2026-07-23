# Acceptance — charging-work-energy

## Dual gate

| Gate | Criterion | Pass when |
| ---- | --------- | --------- |
| **Correctness** | Offline pin match | \(\lvert E_{\mathrm{atv}} - E_{\mathrm{pin}}\rvert \le 10^{-12}\) eV on parallel-plate pin |
| **Perf** | Fair wall vs oracle binary | **perf_exempt** — no MetalWalls/LAMMPS twin binary for this term; pure \(O(N_{\mathrm{elec}})\) sum |

## Checklist

1. **Oracle named in REFERENCE.md**
   - `oracle: literature`
   - `oracle_version:` MetalWalls JOSS 5(53) 2373 (2020) eq.6; Scalfi et al. PCCP 22 10480 (2020) / arXiv:1911.09351 eq.12
2. **Offline pin path**
   - `regressions/md/charging-work-energy/REFERENCE.md` records \(E_{\mathrm{charge}} = -0.4\) eV with geometry table
3. **MD-only build**
   - Configure with `ATV_BUILD_MD=ON`, `ATV_BUILD_KSPACE=ON`, `ATV_BUILD_TESTS=ON`, `ATV_BUILD_REGRESSIONS=ON`, and HF/PW/GPW/CUDA/Torch **OFF**
4. **Observable**
   - `frame.meta.get(kEnergyChargingWork)` and contribution to `kEnergy` [eV]
5. **Tolerance**
   - \(1\times 10^{-12}\) eV absolute
6. **Public API**
   - `cpu::electrodes::ChargingWork` is a `Potential`; no public `*Config`; no `make_*`/`build_*`
7. **Architecture**
   - Component on MD/Evaluate only; no new Driver category; no charge math in `KSpace*`
8. **Unit tests** under `tests/cpu/electrodes/` (not e2e)
9. **One e2e** under `regressions/md/charging-work-energy/`
10. **Regression PASS** vs pin
11. **perf_exempt=true** with reason recorded in REFERENCE.md

## Pin (binding)

| Quantity | Value | Unit |
| -------- | ----: | ---- |
| \(E_{\mathrm{charge}}\) | \(-0.4\) | eV |
| \(\Delta\psi = \Psi_L - \Psi_R\) | \(1.0\) | V |
| \(Q_L = -Q_R\) | \(0.4\) | e |

## Definition of done

`acceptance_met=true` only when every checklist item is verified on
`wf/charging-work-energy` with MD-only build and regression exit 0.
