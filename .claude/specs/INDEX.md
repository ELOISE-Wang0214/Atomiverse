# Spec Index

- [python-bindings](python-bindings.md) — nanobind + DLPack Python bindings; mp.Frame as sole interface, CUDA as compute backend [in-progress: MD/SCF/optimizer/PW bindings landed; RDF observer binding still TODO]
- [eri-hpc-opt](eri-hpc-opt.md) — Shell-block ERI optimization: -march=native, xsimd, shell-type dispatch to match/exceed PySCF direct SCF [in-progress: 5/10 tasks (phases 1-3); phase 4-5 pending]
- [rust-frame-lowering](rust-frame-lowering.md) — Rust-first canonical frame, typed requirement system, Spec/Plan/Driver split, and staged migration roadmap [draft]
- [cpu-qc-product-roadmap](cpu-qc-product-roadmap.md) — **Roadmap index** linking the QC phase sub-specs [in-progress: phase 1 active, phase 2 satisfied via cpu-driver-taxonomy, phases 3-4 draft, phase 5 unscoped]
- [phase1-productization](phase1-productization.md) — Multi-basis + d-shell Rys + 12-molecule benchmark + recorder/observer components + locked defaults [in-progress: 35/38 tasks done; remaining = f-shell Rys + KS gradient d-shell + denser grid]
- [phase3-hybrid-dft](phase3-hybrid-dft.md) — B3LYP/PBE0 via `alpha_hf` parameter extension — zero new drivers, zero wrappers [draft]
- [phase4-vibrations](phase4-vibrations.md) — New `vib::Phase` enum + numerical Hessian driver + thermochemistry [draft]
- [msu-pe-validation](msu-pe-validation.md) — MSU amorphous-PE tension validation vs LAMMPS; correctness + multi-core scaling final acceptance [in-progress: Phase 0/A/B + C landed; A.2 LAMMPS energy cross-check verified 2026-05-20 via c20h42-energy-xcheck; Phase D sliced into the msu-pe-phase-d-* chain below; Phase E (scaling) remains]
- [msu-pe-phase-d-03-equilibrium](msu-pe-phase-d-03-equilibrium.md) — Phase D.3: equilibrium-state observable comparison (density, ⟨E_pot⟩, RDF peak, ⟨R²⟩) over stage-4 100K NPT [in-progress: implemented; needs re-verify after msu-pe-align-02 realigns the benchmark]
- [msu-pe-phase-d-04-stress-strain](msu-pe-phase-d-04-stress-strain.md) — Phase D.4: stress-strain curve comparison (modulus, yield, hardening ±5%) + combined Phase D exit gate [in-progress: implemented; needs re-verify after msu-pe-align-02 realigns the benchmark]
- [msu-pe-align-03-statistical-thermo](msu-pe-align-03-statistical-thermo.md) — replace D.2 per-step thermo parity with statistical block-average comparison (the Langevin stage 1 is stochastic) [approved]
- [cuda-driver-taxonomy](cuda-driver-taxonomy.md) — Mirror the CPU taxonomy on CUDA (cuda::md::MD → DynamicsDriver, cuda::minimizer::SD/CG → Optimizer components, cuda::scf::SCF → SCFSolver inner). Fixes 2 latent CUDA bugs noted during CPU work. [planned]
