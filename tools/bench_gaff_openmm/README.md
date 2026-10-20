# GAFF PEO: Atomiverse CUDA versus OpenMM CUDA

`bench_gaff` and `compare.py` consume the same LAMMPS `units real` data and
coefficients. The harness checks the initial potential energy and every force
component before accepting timing. OpenMM evaluates its own `System`; the
reference numbers are never generated from Atomiverse's formulas.

Build `bench_gaff` with CUDA, MD, potentials, neighbor lists, KSPACE and
benchmarks enabled. No molrs installation or Python binding is required for
the C++ consumer. Run all experiments in an Arrhenius GPU batch allocation:

```sh
export DATA=/path/to/peo/system.data FF=/path/to/peo/system.ff
export RUN_SCRIPT=tools/bench_gaff_openmm/run.sh
sbatch -A <gpu-account> tools/bench_gaff_gpu_job.sh
```

The job uses an isolated `.openmm-env` and records the installed package
versions. `OPENMM_ENV` can select an existing environment. Compilation and
package installation happen before timing. Each engine process has a
10-second deadline; three timed samples use their median. `STEPS` defaults to
100 and `REPLICAS` to `1 3` (replication
along each axis). `PRECISION=mixed BUILD=build-gaff-mixed` builds a separate
mixed-precision binary. Keep both engines on the same GPU.

Matched defaults: unconstrained NVE, 0.25 fs, zero initial velocities, 10 warmup
steps, unshifted Lorentz–Berthelot LJ without tail correction, cutoff 9 Å,
Ewald alpha 0.35 Å⁻¹, order-5 PME, mesh dimensions rounded from 1 Å spacing to
2/3/5/7-smooth sizes. Amber graph-distance exceptions use LJ `(0,0,1/2)` and
Coulomb `(0,0,5/6)`, minimum-image geometry, and no exception cutoff.
Bond, angle and torsion forces use periodic geometry in both engines.
The benchmark reserves a neighbor skin of 8% of the cutoff (0.72 Å), matching
OpenMM CUDA's default padding. This does not change the force cutoff.

Replication follows minimum-image bonded paths and maps boundary-crossing
endpoints to adjacent replicas, preserving bond/angle/torsion geometry.

Supported input: contiguous atom/type IDs, orthorhombic PEO data, harmonic
bonds/angles, Fourier torsions with multiplicities 1–4 and phases 0 or pi.
Those torsions map exactly to the existing `DihedralMultiHarmonic` polynomial.
Unsupported sections, impropers, coefficients or phases raise errors. This
harness does not establish support for every GAFF molecule or NPT workflow.

OpenMM's leapfrog velocity starts one half-kick behind the zero full-step
velocity so the positions follow the same velocity-Verlet trajectory as ATV.
Both timed calls synchronize and retrieve final state. Initial forces/energy
and final periodic positions are checked (final position limits: 1e-4 Å double,
1e-2 Å mixed over the short trajectory). The bounded fixtures explicitly
reserve the full possible cluster CSR in build; only actual neighbors are
evaluated, with no allocation in timed steps. The separate profile
run enables GPU event timing; its wall time is not used in the comparison.
CUDA PME uses real FFTs with a half spectrum. Its fixed-cell influence
coefficients are prepared in build, with conjugate averaging at Nyquist modes
to preserve the full-grid real potential and energy. The separate
`cuda-pme-half-spectrum` regression checks even/odd and skew meshes against
an independent full Fourier sum and analytic force gradient.
Spline/grid/FFT precision follows the selected work precision;
energy and self-energy reductions remain double in both tested modes.
Screening and exception corrections retain the library special functions.
System-level acceptance thresholds are identical across optimization experiments.
Outputs under `logs/gaff-openmm/` include raw forces, energies, actual OpenMM
mesh, input/binary SHA-256 hashes, process durations and ms/step ratios.
`benchmark_seconds` records the whole comparison case, including both engines.
A timing is accepted only when all four initial/final energy, initial force
and final position checks pass (`accepted: true` in the JSON).
Double-mode acceptance uses maximum force error divided by maximum reference
force below `3e-5`, and energy error below `max(1e-3 eV, 3e-5 |E|)`.
Mixed-mode inputs (coordinates and cell) are rounded once to float-representable values for
both engines, so initial parity compares the same positions. Mixed-mode force
tolerance is `1e-4`: ATV stores coordinates in float whereas
OpenMM mixed retains position corrections, so the modes have different
coordinate precision despite the matching names.

Architecture reference: OpenMM source commit
`3a2326d0dd091229eaa4c0d7b5ad1001d5fe376e`, CUDA cluster AABB hierarchy,
nonbonded parameter fusion and sparse exception handling. The implementation
uses ATV's own CSR, masks, Frame and Component lifecycle. The fused force path
packs nonempty J columns across cluster boundaries only when neighbor masks
are rebuilt; timed force calls reuse the build-reserved packed buffers.

Validated Arrhenius GH200 measurements (OpenMM 8.6.1, CUDA 13, three samples
of 100 steps; times in ms/step):

| Atoms | Precision | Atomiverse | OpenMM | ATV / OpenMM | Whole comparison |
| ---: | --- | ---: | ---: | ---: | ---: |
| 1,840 | double | 0.12623 | 0.10162 | 1.242 | 2.01 s |
| 49,680 | double | 0.55344 | 0.62508 | 0.885 | 6.41 s |
| 1,840 | mixed | 0.11730 | 0.07393 | 1.586 | 2.36 s |
| 49,680 | mixed | 0.33558 | 0.27986 | 1.199 | 6.83 s |

All four cases pass the initial/final acceptance checks. Large-system double
is faster than this OpenMM reference; mixed step time remains about 20%
higher. Small systems retain a larger latency gap. These
results cover the PEO force-field subset described above, on one GPU.
