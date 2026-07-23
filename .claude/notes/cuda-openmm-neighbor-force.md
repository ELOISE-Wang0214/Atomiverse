<!-- mol:note:topic:cuda-openmm-neighbor-force -->
## [2026-07-23] CUDA classical MD = OpenMM **engine** (not a style flag)

### Decision (operator, revised same day)

| Layer | Style |
|-------|--------|
| Public / facade API, CPU MD | may stay LAMMPS-**named** |
| **CUDA** neighbor rebuild + nonbonded force | **OpenMM cluster/tile engine end-to-end** |

This is an **engine architecture** rewrite, not `ATV_CUDA_NONBONDED=openmm|lammps`
or “flip default back to cluster CSR”.

**OpenMM end-to-end on CUDA means:**

1. Rebuild: spatial reorder (Morton/Hilbert) → **clusters of 32** → AABB →
   cluster-pair neighbor CSR → **real interaction masks** (never production
   `memset 0xFF`).
2. Force: **tile / shared-memory** kernels over cluster pairs for LJ and
   real-space Coulomb; PME remains separate (`KSpacePME` / cuFFT).
3. multi-GPU consumes the **same** cluster contract (no `vl_*` as force truth).
4. Delete Verlet-direct as CUDA MD default after the rewrite lands.

### Live specs

`.claude/specs/cuda-openmm-md-01-contract.md` … `05-delete.md`  
(supersedes `.claude/specs/cuda-openmm-neighbor-force.md`)

### Forbidden narratives

- “We A/B’d our half OpenMM shell vs optimized Verlet, so OpenMM is slower.”
- Dual-path `NonbondedStyle` as the architecture.
- Treating `cuda-pme-accuracy` CoulLong→`vl_*` as long-term design.

### Performance

Measure after tile force lands (phase 03+): rebuild ms, force ms, step/s on
GH200; record in notes. Micro-opts **inside** OpenMM engine only.
<!-- mol:note:topic:cuda-openmm-neighbor-force -->
