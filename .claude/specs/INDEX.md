# Spec Index

- [python-bindings](python-bindings.md) — nanobind + DLPack Python bindings; mp.Frame as sole interface, CUDA as compute backend [in-progress: MD/SCF/optimizer/PW bindings landed; RDF observer binding still TODO]
- [libxc-xc-fix](libxc-xc-fix.md) — Integrate libxc as ground truth, fix native XC functional + grid to match within 1e-12 [draft: 9 tasks, 8 acceptance criteria]
- [eri-hpc-opt](eri-hpc-opt.md) — Shell-block ERI optimization: -march=native, xsimd, shell-type dispatch to match/exceed PySCF direct SCF [in-progress: 5/10 tasks (phases 1-3); phase 4-5 pending]
- [rust-frame-lowering](rust-frame-lowering.md) — Rust-first canonical frame, typed requirement system, Spec/Plan/Driver split, and staged migration roadmap [draft]
- [cpu-qc-product-roadmap](cpu-qc-product-roadmap.md) — **Roadmap index** linking the QC phase sub-specs [in-progress: phase 1 active, phase 2 satisfied via cpu-driver-taxonomy, phases 3-4 draft, phase 5 unscoped]
- [phase1-productization](phase1-productization.md) — Multi-basis + d-shell Rys + 12-molecule benchmark + recorder/observer components + locked defaults [in-progress: 35/38 tasks done; remaining = f-shell Rys + KS gradient d-shell + denser grid]
- [phase3-hybrid-dft](phase3-hybrid-dft.md) — B3LYP/PBE0 via `alpha_hf` parameter extension — zero new drivers, zero wrappers [draft]
- [phase4-vibrations](phase4-vibrations.md) — New `vib::Phase` enum + numerical Hessian driver + thermochemistry [draft]
- [parallel-md-engine](parallel-md-engine.md) — MPI + OpenMP + xsimd multi-node CPU MD; hot-path vs. analysis separation via `std::thread` pool [in-progress: phases 1–4 landed; **phase 5 partial 2026-05-06** — 15 of 27 tasks done (incl. B.2 LJ xsimd 1.21× NEON, B.3 bonded stack-pack), 7 partial (reference-machine + B.2-bench-guard infeasible), 5 pending (B.4 cell-list cadence + B.9 LAMMPS perf measurement)]
- [msu-pe-validation](msu-pe-validation.md) — MSU amorphous-PE tension validation vs LAMMPS; correctness + multi-core scaling final acceptance [in-progress: 0/21 tasks done]
- [cuda-driver-taxonomy](cuda-driver-taxonomy.md) — Mirror the CPU taxonomy on CUDA (cuda::md::MD → DynamicsDriver, cuda::minimizer::SD/CG → Optimizer components, cuda::scf::SCF → SCFSolver inner). Fixes 2 latent CUDA bugs noted during CPU work. [planned]
