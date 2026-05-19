# Spec Index

- [python-bindings](python-bindings.md) — nanobind + DLPack Python bindings; mp.Frame as sole interface, CUDA as compute backend [in-progress]
- [libxc-xc-fix](libxc-xc-fix.md) — Integrate libxc as ground truth, fix native XC functional + grid to match within 1e-12 [draft]
- [eri-hpc-opt](eri-hpc-opt.md) — Shell-block ERI optimization: -march=native, xsimd, shell-type dispatch to match/exceed PySCF direct SCF [draft]
- [rust-frame-lowering](rust-frame-lowering.md) — Rust-first canonical frame, typed requirement system, Spec/Plan/Driver split, and staged migration roadmap [draft]
- [cpu-qc-product-roadmap](cpu-qc-product-roadmap.md) — **Roadmap index** linking the QC phase sub-specs [draft]
- [phase1-productization](phase1-productization.md) — Multi-basis + d-shell Rys + 12-molecule benchmark + recorder/observer components + locked defaults [draft]
- [phase3-hybrid-dft](phase3-hybrid-dft.md) — B3LYP/PBE0 via `alpha_hf` parameter extension — zero new drivers, zero wrappers [draft]
- [phase4-vibrations](phase4-vibrations.md) — New `vib::Phase` enum + numerical Hessian driver + thermochemistry [draft]
- [parallel-md-engine](parallel-md-engine.md) — MPI + OpenMP + xsimd multi-node CPU MD; hot-path vs. analysis separation via `std::thread` pool [in-progress: phases 1–4 landed; **phase 5 partial 2026-05-06** — 15 of 27 tasks done (incl. B.2 LJ xsimd 1.21× NEON, B.3 bonded stack-pack), 7 partial (reference-machine + B.2-bench-guard infeasible), 5 pending (B.4 cell-list cadence + B.9 LAMMPS perf measurement)]
- [msu-pe-validation](msu-pe-validation.md) — MSU amorphous-PE tension validation vs LAMMPS; correctness + multi-core scaling final acceptance [draft]
- [cuda-driver-taxonomy](cuda-driver-taxonomy.md) — Mirror the CPU taxonomy on CUDA (cuda::md::MD → DynamicsDriver, cuda::minimizer::SD/CG → Optimizer components, cuda::scf::SCF → SCFSolver inner). Fixes 2 latent CUDA bugs noted during CPU work. [planned]
