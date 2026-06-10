# pw — plane-wave DFT benchmarks (bulk Si vs CP2K)

Benchmarks the plane-wave SCF (`cpu::pw::pw_scf` / `cuda::pw::pw_scf`)
against CP2K Quickstep on bulk silicon (2-atom fcc primitive cell,
LDA). Layout follows the standard scheme:

```
validation/   science: Si equation of state (a0, B0) vs CP2K
performance/  speed: SCF time-to-solution CPU vs CUDA vs CP2K
```

Build with `-DATV_BUILD_PW=ON -DATV_BUILD_RHF=ON -DATV_BUILD_BENCHMARKS=ON`
(add `-DATV_BUILD_CUDA=ON` on the GPU node for the cuda backend rows).

## validation — equation of state (scientific part)

`bench_pw_si_cp2k_eos` scans a = 5.25 … 5.55 Å, fits E(V) by least
squares (quadratic in V), and compares the equilibrium lattice
constant a0 and bulk modulus B0 against the CP2K reference fitted with
the same code. PASS tolerances: |Δa0| ≤ 0.08 Å, |ΔB0| ≤ 25 GPa —
generous enough to absorb the legitimate cross-code differences
(pseudopotential construction, LDA parametrization, basis), tight
enough to catch real physics bugs (wrong Ewald, broken k-sampling,
XC/Hartree double counting all shift a0/B0 far more).

```bash
./bench_pw_si_cp2k_eos                  # CPU, ecut 30 Ry, 2x2x2 MP (smoke)
./bench_pw_si_cp2k_eos --backend cuda --ecut 60 --kmesh 4   # production
```

Production settings are `--ecut 60 --kmesh 4`: the ONCV pseudo needs
~60 Ry for sub-mRy total-energy convergence, and the 2×2×2 mesh leaves
several-percent a0 errors. On the GH200 the production scan takes
~15 min; on a single CPU core it is hours — run the science part on
the cuda backend (the parity test guarantees both backends produce
identical energies).

The CP2K side is generated once per machine with
`validation/reference/cp2k/run_eos.sh /path/to/cp2k.ssmp` (any recent
CP2K; the deck inlines basis + pseudopotential, so no data directory
is required). Until that CSV exists the bench compares against the
same-pseudopotential ABINIT deltafactor reference (a0 = 5.469 Å,
B0 = 90 GPa — see `validation/reference/SOURCES.md`) and prints a
notice.

Why these observables: total energies are not comparable across codes
(different pseudopotentials/discretizations); a0 and B0 are the
standard pseudopotential-insensitive EOS observables. LDA reference:
a0 ≈ 5.40 Å, B0 ≈ 96 GPa (experiment: 5.431 Å, 98.8 GPa).

## performance — SCF time-to-solution

`bench_pw_scf_perf` times complete SCF runs (tol 1e-7, same physical
system) and appends CSV rows; `performance/cp2k/run_perf.sh` appends
backend=cp2k rows in the same schema:

```
backend,n_atoms,a_angstrom,ecut_ry,k_mesh,threads,iterations,converged,energy_ry,wall_s,s_per_iter
```

```bash
./bench_pw_scf_perf --backend cpu  --ecut-list 20,30,40 --out pw_scf_perf.csv
./bench_pw_scf_perf --backend cuda --ecut-list 20,30,40 --out pw_scf_perf.csv
performance/cp2k/run_perf.sh ~/bin/cp2k.ssmp pw_scf_perf.csv
```

Protocol & caveats (read before quoting numbers):

- The compared quantity is **wall time to a converged total energy of
  the same physical system**, each code at its own production-typical
  discretization — not FLOP-for-FLOP parity. CP2K's `CUTOFF` is the
  *density*-grid cutoff (Ry) of GPW, not the plane-wave wavefunction
  cutoff in the atomiverse rows; do not compare rows at equal
  `ecut_ry` across backends.
- Atomiverse pw is single-threaded per process; CP2K ssmp uses OpenMP
  (`threads` column records it). Per-core comparisons should divide
  accordingly.
- Run CPU rows for both codes on the same node. The CUDA rows
  (GH200) quantify the device speedup of the Davidson hot path over
  the single-threaded CPU backend.
