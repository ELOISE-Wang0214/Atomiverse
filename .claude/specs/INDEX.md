# Spec Index

- [aimd](aimd.md) — first-class Born-Oppenheimer AIMD by reusing cpu::md::MD + SCFEvaluator: set_qm categorical setter + faithful Kolafa-ASPC density warm-start + NVE energy-conservation validation [approved]

- [rust-frame-lowering-01-molrs-bridge](rust-frame-lowering-01-molrs-bridge.md) — molrs CXX bridge: add Frame zero-copy ingress + create/write, expose the Frame handle from molrs-python [approved]
- [rust-frame-lowering-02-cpu-materialize](rust-frame-lowering-02-cpu-materialize.md) — cpu::Frame materializer for molrs.Frame ↔ cpu::Frame, add charge field, extract shared SystemScalars [approved]
- [rust-frame-lowering-03-python-path](rust-frame-lowering-03-python-path.md) — run() accepts molrs.Frame via a shared handle; delete _convert.py and the CpuFrame binding [approved]
- [rust-frame-lowering-04-cuda-torch](rust-frame-lowering-04-cuda-torch.md) — remove SCFWorkspace/PairList from cuda::Frame, align the three Frames, add cuda/torch materializers [approved]
- [msu-pe-validation](msu-pe-validation.md) — MSU amorphous-PE tension validation vs LAMMPS; correctness + multi-core scaling final acceptance [in-progress: Phase 0/A/B + C landed; A.2 LAMMPS energy cross-check verified 2026-05-20 via c20h42-energy-xcheck; Phase D sliced into the msu-pe-phase-d-* chain below; Phase E (scaling) remains]
- [msu-pe-phase-d-03-equilibrium](msu-pe-phase-d-03-equilibrium.md) — Phase D.3: equilibrium-state observable comparison (density, ⟨E_pot⟩, RDF peak, ⟨R²⟩) over stage-4 100K NPT [in-progress: implemented; needs re-verify after msu-pe-align-02 realigns the benchmark]
- [msu-pe-phase-d-04-stress-strain](msu-pe-phase-d-04-stress-strain.md) — Phase D.4: stress-strain curve comparison (modulus, yield, hardening ±5%) + combined Phase D exit gate [in-progress: implemented; needs re-verify after msu-pe-align-02 realigns the benchmark]
- [msu-pe-align-03-statistical-thermo](msu-pe-align-03-statistical-thermo.md) — replace D.2 per-step thermo parity with statistical block-average comparison (the Langevin stage 1 is stochastic) [code-complete: 7/8 verified; ac-008 needs the operator reference run]
- [cuda-driver-taxonomy-02-md](cuda-driver-taxonomy-02-md.md) — reparent cuda::md::MD onto cuda::DynamicsDriver, make Zero* phase-agnostic, tag MD components role()/capabilities() [approved]
- [cuda-driver-taxonomy-03-minimizer](cuda-driver-taxonomy-03-minimizer.md) — delete the SD/CG driver classes, promote SteepestDescent to an Optimizer component, add cuda::minimizer::Minimize [approved]
- [cuda-driver-taxonomy-04-scf](cuda-driver-taxonomy-04-scf.md) — demote cuda::scf::SCF to the SCFSolver inner solver (internals → cuda::scf::internal::), add cuda::scf::SCFEvaluator + cuda::Evaluate + cuda::single_point [approved]
- [cuda-driver-taxonomy-05-scrub-and-docs](cuda-driver-taxonomy-05-scrub-and-docs.md) — scrub dangling CUDA scf::Phase refs (enum kept — cpu::pw plane-wave SCF owns it), update CLAUDE.md + notes [approved]
