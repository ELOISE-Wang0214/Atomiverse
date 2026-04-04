---
title: Frame Block-Model Refactor (molcrafts naming alignment)
status: draft
created: 2026-04-04
---

## Summary

Restructure `cpu::Frame` and align `cuda::Frame` naming to match molpy's dict-of-blocks model. Field names follow molcrafts conventions. This is a prerequisite for the Python binding layer (see `python-bindings.md`).

**Goal:** `frame.atoms.x`, `frame.bonds.atomi`, `frame.atoms.type` — same names as molpy blocks, on both CPU and CUDA.

## Current state vs target

### cpu::Frame (needs restructuring)

```cpp
// CURRENT: flat struct, Atomiverse-specific names
struct Frame {
    xt::xtensor<Real, 1> x, y, z, vx, vy, vz, fx, fy, fz, mass;
    xt::xtensor<Int, 1>  type_id;                       // ← diverges from molpy
    xt::xtensor<Int, 1>  bond_i, bond_j, bond_type;     // ← flat, diverges
    xt::xtensor<Int, 1>  angle_i, angle_j, angle_k, angle_type;
    xt::xtensor<Int, 1>  dihedral_i, dihedral_j, dihedral_k, dihedral_l, dihedral_type;
    Real energy, kinetic, stress[6];
    Mat3 box;
};

// TARGET: sub-structs matching molpy blocks, molcrafts names
struct Frame {
    struct Atoms {
        xt::xtensor<Real, 1> x, y, z;
        xt::xtensor<Real, 1> vx, vy, vz;
        xt::xtensor<Real, 1> fx, fy, fz;
        xt::xtensor<Real, 1> mass;
        xt::xtensor<Int, 1>  type;          // was type_id
        Int size() const;                    // was n_atoms()
    } atoms;

    struct Bonds {
        xt::xtensor<Int, 1> atomi, atomj;   // was bond_i, bond_j
        xt::xtensor<Int, 1> type;            // was bond_type
        Int size() const;
    } bonds;

    struct Angles {
        xt::xtensor<Int, 1> atomi, atomj, atomk;  // was angle_i/j/k
        xt::xtensor<Int, 1> type;
        Int size() const;
    } angles;

    struct Dihedrals {
        xt::xtensor<Int, 1> atomi, atomj, atomk, atoml;
        xt::xtensor<Int, 1> type;
        Int size() const;
    } dihedrals;

    Real energy = 0, kinetic = 0;
    Real stress[6] = {};
    Mat3 box{};
    const void* pair_list = nullptr;
};
```

### cuda::Frame (mostly aligned, needs renames)

```cpp
// CURRENT: per-atom flat, topology already sub-structed
struct Frame {
    Real *x, *y, *z, *vx, *vy, *vz, *fx, *fy, *fz, *mass;
    Int* type_id;           // ← rename to type
    // ...
    struct BondList {       // already sub-struct
        Int *i, *j, *type;  // ← rename i→atomi, j→atomj
    } bonds;
    struct AngleList {
        Int *i, *j, *k, *type;  // ← rename
    } angles;
    struct DihedralList {
        Int *i, *j, *k, *l, *type;  // ← rename
    } dihedrals;
};

// TARGET: group per-atom into atoms sub-struct, rename topology fields
struct Frame {
    struct Atoms {
        Real *x, *y, *z, *vx, *vy, *vz, *fx, *fy, *fz, *mass;
        Int *type;
        Int* d_size;        // was d_n_atoms
    } atoms;

    struct Bonds {
        Int *atomi, *atomj, *type;   // was i, j, type
        Int* d_size;
    } bonds;

    struct Angles {
        Int *atomi, *atomj, *atomk, *type;
        Int* d_size;
    } angles;

    struct Dihedrals {
        Int *atomi, *atomj, *atomk, *atoml, *type;
        Int* d_size;
    } dihedrals;

    PairList pairs;          // unchanged
    SCFWorkspace scf;        // unchanged

    Real *energy, *kinetic, *stress, *box;  // device scalars, unchanged

    // Host-cached metadata
    Int n_atoms_, n_bonds_, n_angles_, n_dihedrals_;
    Int n_atoms() const { return n_atoms_; }
};
```

## Rename mapping

All renames are mechanical find-and-replace within the sub-struct context.

| scope | old | new | files |
|-------|-----|-----|-------|
| per-atom | `frame.x` | `frame.atoms.x` | ~42 |
| per-atom | `frame.vx` | `frame.atoms.vx` | ~42 |
| per-atom | `frame.fx` | `frame.atoms.fx` | ~42 |
| per-atom | `frame.mass` | `frame.atoms.mass` | ~20 |
| per-atom | `frame.type_id` | `frame.atoms.type` | 9 |
| per-atom | `frame.n_atoms()` | `frame.atoms.size()` (cpu) / `frame.n_atoms()` (cuda, keep) | all |
| bonds | `frame.bond_i` | `frame.bonds.atomi` | 4 (cpu) |
| bonds | `frame.bonds.i` | `frame.bonds.atomi` | cuda (already sub-struct) |
| angles | `frame.angle_i` | `frame.angles.atomi` | 4 (cpu) |
| angles | `frame.angles.i` | `frame.angles.atomi` | cuda |
| dihedrals | `frame.dihedral_i` | `frame.dihedrals.atomi` | 4 (cpu) |
| dihedrals | `frame.dihedrals.i` | `frame.dihedrals.atomi` | cuda |

## Execution plan

### Phase 1: cpu::Frame header + transfer (low risk, unlocks Python binding)

1. Rewrite `include/atomiverse/cpu/frame.h` with sub-structs
2. Update `src/cuda/transfer.cu` (cpu::to_cuda, cuda::to_cpu)
3. Update `src/torch/transfer.cu`
4. Build — all downstream files will fail with clear errors

### Phase 2: mechanical rename (bulk, parallelizable by domain)

Domain batches (independent, can be done in parallel):

| batch | files | pattern |
|-------|-------|---------|
| MD integrators | velocity_verlet, init, kinetic_energy, nose_hoover, barostat, deform | `f.x` → `f.atoms.x` |
| Potentials | lj, harmonic_bond, harmonic_angle, multi_harmonic_dihedral | `f.x` → `f.atoms.x`, bonds/angles/dihedrals |
| Neighbor list | neighbor_list, sort_and_cluster, grid, pair_search, fused_search, export_pairs | `pos_x` → stays (kernel args, not Frame fields) |
| SCF | hamiltonian, eigen_solver, density_builder, convergence_check, rhf_fock, rhf_density | `f.x` → `f.atoms.x` |
| CPU backend | all `src/cpu/` files | same pattern |
| Observers/recorders | thermo, stress, rdf, density, mulliken, trajectory, xyz | `f.x` → `f.atoms.x` |

### Phase 3: cuda::Frame topology renames

1. Rename `BondList::i` → `BondList::atomi` (etc.) in `cuda/frame.cuh`
2. Update all CUDA potentials that access `bonds.i`, `angles.i`, `dihedrals.i`
3. Update `cuda/transfer.cu`

### Phase 4: verify

1. `cmake --build build -j$(nproc)` — zero errors
2. `ctest --test-dir build --output-on-failure` — all tests pass
3. Examples still run

## Scope estimate

| metric | count |
|--------|-------|
| Header changes | 3 (cpu/frame.h, cuda/frame.cuh, torch/frame.cuh) |
| Source files touched | ~60-70 |
| Nature of changes | Mechanical rename (`.x` → `.atoms.x`, `.bond_i` → `.bonds.atomi`) |
| Risk | Low — compile errors catch all misses; no logic change |
| Test coverage | Existing tests validate correctness post-rename |

## Design decisions

1. **`atoms.size()` vs `n_atoms()`**: cpu::Frame uses `atoms.size()` (block-native). cuda::Frame keeps `n_atoms()` (host-cached, avoids D2H). Both return the same value.

2. **Scalars stay flat**: `energy`, `kinetic`, `stress` remain top-level Frame fields, not inside a block. They're system-wide, not per-block. molpy puts them in `metadata` dict, which is unstructured — no sub-struct needed in C++.

3. **`pair_list` stays top-level**: PairList is runtime state set by NeighborList, not persistent block data. No molpy equivalent.

4. **SCFWorkspace stays top-level**: Electronic structure workspace is driver-specific context, not block data.

5. **No constructor change for cpu::Frame**: Still takes an atom count to allocate arrays. The count is `atoms.size()` after construction, not a top-level `n_atoms_` field. Topology blocks default to empty.

## Non-goals

- Changing storage layout (SoA → AoS). Layout stays as-is; only names change.
- Adding new fields. This is a pure rename/restructure.
- Changing cuda::Frame ownership model. Device memory lifecycle unchanged.
- Python binding work. That's `python-bindings.md`.
