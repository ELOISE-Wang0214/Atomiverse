# Getting Started

## Installation

**From PyPI (stable):**

```bash
pip install atomiverse
```

**From source (development):**

```bash
git clone https://github.com/molcrafts/Atomiverse
cd Atomiverse
pip install --no-build-isolation -ve python/
```

Requirements: Python ≥ 3.9, NumPy, C++17 compiler, CMake ≥ 3.18.

Optional: CUDA toolkit for GPU backend (`ATV_BUILD_CUDA=ON pip install python/`).

---

## Quickstart: Argon MD

Run a Lennard-Jones argon simulation in five steps.

```python
import numpy as np
import atomiverse as av

# 1. Build an FCC lattice of 500 Ar atoms
n_cells = 5
a0 = 5.26        # lattice constant [Å]
n = n_cells ** 3 * 4
basis = [(0,0,0),(0.5,0.5,0),(0.5,0,0.5),(0,0.5,0.5)]
positions = [
    ((ix + bx) * a0, (iy + by) * a0, (iz + bz) * a0)
    for ix in range(n_cells) for iy in range(n_cells) for iz in range(n_cells)
    for bx, by, bz in basis
]
xs, ys, zs = (np.array(c, dtype=np.float32) for c in zip(*positions))

frame = av.CpuFrame(n)
frame.set_positions(xs, ys, zs)
frame.set_mass(np.full(n, 39.948, dtype=np.float32))   # Ar [amu]
frame.set_type(np.zeros(n, dtype=np.int32))
frame.set_box(np.eye(3, dtype=np.float32) * (n_cells * a0))

# 2. Configure the MD driver
epsilon, sigma, cutoff = 0.0104, 3.4, 8.5   # LJ-Ar [eV, Å, Å]
dt = 0.001   # [fs]

sim = av.MD()
sim.add_velocity_verlet(av.VelocityVerlet(dt))
sim.add_init_velocity(av.InitVelocity(temperature=10.0, seed=42))
sim.on(av.Phase.force) \
   .add_neighbor_list(av.NeighborList(cutoff), order=-10) \
   .add_lj(av.LJ(epsilon, sigma, cutoff))
sim.on(av.Phase.end_of_step) \
   .add_thermo(av.Thermo(freq=100))

# 3. Build and run
sim.build()
sim.run(frame, n_steps=1000)

# 4. Read results
print(f"PE = {frame.energy:.4f} eV")
print(f"KE = {frame.kinetic:.4f} eV")
print(f"x[0] = {frame.get_x()[0]:.3f} Å")
```

---

## Key Concepts

### CpuFrame

`CpuFrame` holds system state as NumPy arrays. All `set_*` methods copy data
in; all `get_*` methods copy data out (arrays are owned by the frame).

```python
frame = av.CpuFrame(n_atoms)
frame.set_positions(x, y, z)   # float32 arrays, length n_atoms
frame.set_mass(mass)            # float32
frame.set_type(types)           # int32
frame.set_box(box)              # float32, shape (3, 3)
```

Topology (optional):

```python
frame.set_bonds(atomi, atomj, btype)
frame.set_angles(atomi, atomj, atomk, atype)
frame.set_dihedrals(atomi, atomj, atomk, atoml, dtype)
```

### MD Driver

`MD` is the simulation driver. Components are registered with `add_*` (auto
phase) or `on(phase).add_*(component, order=0)` (explicit phase).

| Component | Phase | Description |
|-----------|-------|-------------|
| `VelocityVerlet(dt)` | initial/final integrate | Symplectic integrator |
| `InitVelocity(T, seed)` | build-time | Maxwell-Boltzmann init |
| `NoseHoover(T, tau, dt)` | initial/final integrate | NVT thermostat |
| `NeighborList(cutoff)` | force (-10) | Cell-list pair builder |
| `LJ(ε, σ, cutoff)` | force | Lennard-Jones 12-6 |
| `HarmonicBond(K, r0)` | force | Bond stretching |
| `HarmonicAngle(K, θ0)` | force | Angle bending |
| `MultiHarmonicDihedral(Cn)` | force | Dihedral potential |
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

---

## Thermostats

NVT ensemble with Nosé-Hoover:

```python
sim = av.MD()
sim.add_velocity_verlet(av.VelocityVerlet(dt=0.001))
sim.add_init_velocity(av.InitVelocity(T=300.0, seed=42))
sim.add_nose_hoover(av.NoseHoover(T_target=300.0, tau=100.0, dt=0.001))
sim.on(av.Phase.force) \
   .add_neighbor_list(av.NeighborList(cutoff), -10) \
   .add_lj(av.LJ(epsilon, sigma, cutoff))
sim.build()
sim.run(frame, n_steps=10000)
```

---

## Extracting Trajectories

Results are read back from the frame after `run()`:

```python
sim.run(frame, n_steps)

x  = frame.get_x()   # np.ndarray float32, shape (n_atoms,)
vx = frame.get_vx()
fx = frame.get_fx()
PE = frame.energy     # float
KE = frame.kinetic    # float
```

---

## Timing

```python
sim.enable_timing()
sim.build()
sim.run(frame, 1000)
print(sim.summary())   # per-phase wall-clock breakdown
```

---

## Next Steps

- **Examples** — `python/examples/` for PySCF + Atomiverse comparison experiments
- **Benchmarks** — `python/benchmarks/bench_md.py` for system-size scaling
- **CUDA backend** — set `ATV_BUILD_CUDA=ON` when building on a GPU host
