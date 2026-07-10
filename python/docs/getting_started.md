# Getting Started

## Installation

**From source (development):**

```bash
git clone https://github.com/molcrafts/Atomiverse
cd Atomiverse
pip install --no-build-isolation -ve python/
```

Requirements: Python ≥ 3.9, NumPy, `molrs` (the frame/data layer — installed
separately), C++17 compiler, CMake ≥ 3.18, a Rust toolchain (molrs is built
via corrosion).

Optional backends at build time:

- `CMAKE_ARGS='-DATV_BUILD_PW=ON'` — plane-wave DFT (`pw_scf`; needs FFTW3).
- `CMAKE_ARGS='-DATV_BUILD_CUDA=ON'` — CUDA backend (GPU host only).

## What this package exposes (and what it doesn't)

Atomiverse binds **drivers and components only** — there is no Atomiverse
Frame type. Simulation state crosses the boundary as a **`molrs.Frame`**:
each `run()` borrows the frame, advances a transient execution frame in C++,
and writes positions / velocities / forces / energy back into the same
`molrs.Frame`.

Current surface:

- **MD** — the `MD` driver and its components under the canonical C++
  (LAMMPS-style) names: `VelocityVerlet`, `InitVelocity`, `NoseHoover`,
  `PairLJCut`, `NeighborList`, `BondHarmonic`, `AngleHarmonic`,
  `DihedralMultiHarmonic`, `Thermo`.
- **Plane-wave DFT** — `pw_scf` / `PWResult` (PW-enabled wheels only; a
  default wheel raises `RuntimeError` with rebuild guidance).

**Roadmap:** Gaussian-SCF and geometry-optimization bindings (the C++
`SCFEvaluator` family and the `Minimize` driver with `LBFGS`/`CG`
components) are future work with their own spec. They are deliberately not
exposed rather than exposed broken.

## Quickstart: Argon MD

```python
import numpy as np
import molrs
import atomiverse as av

# 1. Build an FCC lattice of 500 Ar atoms as a molrs.Frame
n_cells = 5
a0 = 5.26        # lattice constant [Å]
basis = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]
positions = [
    ((ix + bx) * a0, (iy + by) * a0, (iz + bz) * a0)
    for ix in range(n_cells) for iy in range(n_cells) for iz in range(n_cells)
    for bx, by, bz in basis
]
n = len(positions)
xs, ys, zs = (np.array(c, dtype=np.float64) for c in zip(*positions))

frame = molrs.Frame()
atoms = molrs.Block()
atoms.insert("x", xs)
atoms.insert("y", ys)
atoms.insert("z", zs)
atoms.insert("mass", np.full(n, 39.948))   # Ar [amu]
atoms.insert("type", np.zeros(n, dtype=np.int32))
frame["atoms"] = atoms
frame.box = molrs.Box.cube(n_cells * a0)

# 2. Configure the MD driver
epsilon, sigma, cutoff = 0.0104, 3.4, 8.5   # LJ-Ar [eV, Å, Å]

sim = av.MD(dt=0.001)                        # timestep [fs]
sim.add(av.InitVelocity(temperature=10.0, seed=42))
sim.on("force").add(av.NeighborList(cutoff), order=-10) \
               .add(av.PairLJCut(epsilon, sigma, cutoff))
sim.on("end_of_step").add(av.Thermo(100))

# 3. Run (build() is called automatically)
sim.run(frame, 1000)

# 4. Read results back from the SAME molrs.Frame
x = np.asarray(frame["atoms"]["x"])
pe = float(np.asarray(frame["scalars"]["energy"])[0])
ke = float(np.asarray(frame["scalars"]["kinetic"])[0])
print(f"PE = {pe:.4f} eV,  KE = {ke:.4f} eV,  x[0] = {x[0]:.3f} Å")
```

## Key concepts

### molrs.Frame in, molrs.Frame out

The frame must carry an `atoms` block with `x, y, z, type` columns (`mass`
required for MD). After `run()` the frame holds updated `x/y/z`, `vx/vy/vz`,
`fx/fy/fz`, and the system scalars (`scalars` block: `energy`, `kinetic`,
`conserved`). Missing required fields fail loud with a `RuntimeError` naming
the field.

### MD driver

Components register through two paths, mirroring the C++ API:

- `sim.add(comp)` — components with default phase bindings
  (`VelocityVerlet`, `InitVelocity`, `NoseHoover`).
- `sim.on(phase).add(comp, order=0)` — explicit phase binding, where
  *phase* is one of `"initial_integrate"`, `"pre_force"`, `"force"`,
  `"post_force"`, `"final_integrate"`, `"end_of_step"`.

| Component | Phase | Description |
|-----------|-------|-------------|
| `VelocityVerlet(dt)` | initial/final integrate | Symplectic integrator |
| `InitVelocity(T, seed)` | build-time | Maxwell-Boltzmann init |
| `NoseHoover(T_target, tau_T, dt)` | initial/final integrate | NVT thermostat |
| `NeighborList(cutoff, skin=-1)` | force (order −10) | Cell-list pair builder |
| `PairLJCut(ε, σ, cutoff)` | force | Lennard-Jones 12-6 |
| `BondHarmonic(K, r0)` | force | Bond stretching |
| `AngleHarmonic(K, θ0)` | force | Angle bending |
| `DihedralMultiHarmonic(Cn)` | force | Dihedral potential |
| `Thermo(freq)` | end_of_step | Thermodynamic output |

### Units

| Quantity | Unit |
|----------|------|
| Distance | Å (Angstrom) |
| Energy | eV |
| Time | fs |
| Mass | amu |
| Temperature | K |
| Force | eV/Å |

## Plane-wave DFT (PW-enabled wheels)

```python
result = av.pw_scf(frame, "Si.upf", ecut_ry=30.0, k_mesh=4)
print(result.energy)      # [Ry]
print(result.converged, result.iterations)
```

An empty `upf_path` runs a toy local potential (useful as a smoke test;
results are not physical). On a default wheel (`ATV_BUILD_PW=OFF`) the call
raises `RuntimeError` telling you to rebuild with `ATV_BUILD_PW=ON`.

## Timing

```python
sim.enable_timing()
sim.build()
sim.run(frame, 1000)
print(sim.summary())   # per-phase wall-clock breakdown
```

## Next steps

- **Examples** — `python/examples/pw_qe_compare.py` (plane-wave DFT vs
  Quantum ESPRESSO for Si diamond).
- **Benchmarks** — `python/benchmarks/bench_md.py` for system-size scaling.
- **Tests** — `python/tests/test_surface.py` documents the exact exported
  surface and its fail-loud contract.
