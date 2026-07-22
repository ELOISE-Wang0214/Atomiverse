---
title: CUDA SCF gradients (CUDA SCF is energy/density-only today)
slug: hf-parity-20-cuda-scf-gradients
status: draft
chain: hf-parity
created: 2026-07-22
note: >
  Roadmap draft from hf-deferred-backlog P5. GPU-node work (aarch64 GH200 — see
  CLAUDE.md GPU build recipe). Run `/mol:spec hf-parity-20-cuda-scf-gradients`.
  science.required.
---

# CUDA SCF gradients

## Summary
CUDA SCF is energy/density-oriented; there is no CUDA SCF gradient. Mirror the
CPU analytic gradient (incl. the shipped spherical c2s path) on the CUDA
backend so AIMD/opt with SCF forces runs on GPU. Parent: hf-deferred-backlog P5.

## Domain basis
Same gradient theory as CPU (Pulay + Hellmann-Feynman + 2e derivative); the work
is porting the derivative-integral + contraction kernels to CUDA and keeping the
CUDA SCF solver's narrow `SolverContext` isolation (per CLAUDE.md). Orthogonality
rule: gradients live in the same directory on both backends.

## Design (direction — GPU node only)
- Port `IntegralDerivPack` + `accumulate_two_electron_gradient` analogues to
  `cuda/scf/`; reuse device ERI-derivative kernels; c2s sandwich on device.

## Files (best-effort)
- `cuda/scf/` gradient kernels + `internal/solver_context.cuh`
- `tests/cuda/`; `regressions/scf/` cross-backend parity (CPU vs CUDA gradient)

## Tasks (provisional)
- [ ] RED: CUDA gradient vs CPU gradient parity on a tiny system (GPU test).
- [ ] Port 1e derivative + 2e contraction kernels.
- [ ] c2s sandwich on device for spherical.
- [ ] Cross-backend parity regression.

## Testing
GPU-node build (`GPU/buildenv-gcccuda/2026.03-cu13.0`); CUDA-vs-CPU gradient
parity ≤ ~1e-8. Runs only on the GH200 node.

## Out of scope
CUDA KS-XC gradient depth (fold in if the XC grid is ported).

## Depends on
Shipped CPU spherical gradient; existing CUDA SCF energy/density path.
