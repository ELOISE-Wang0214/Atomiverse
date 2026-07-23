---
spec: cuda-pme-accuracy
created: 2026-07-23
criteria:
  - id: ac-001
    summary: CUDA PME reproduces the NaCl Madelung energy
    type: runtime
    pass_when: |
      On the GH200 node, `ctest -R test_cuda_pme_coul` passes:
      CudaPmeCoul.MadelungNaCl has |E/E_expected - 1| <= 3e-3 with the existing
      parameters (d=2.82, n=4, cutoff=1.9d, alpha=3.5/cutoff, grid_spacing=0.5,
      order=6) and the DEFAULT precision (ATV_CUDA_REAL=float). The tolerance is
      NOT relaxed — CPU reaches 5e-4 on identical parameters.
    status: pending
  - id: ac-002
    summary: the fix is localized to the CUDA reciprocal-space path
    type: code
    pass_when: |
      Only src/cuda/longrange/ (pme.cu and, if needed, its header) changes.
      No CPU code, no change to ATV_CUDA_REAL defaults, no test-tolerance edit.
    status: pending
  - id: ac-003
    summary: the Ewald sibling stays green
    type: runtime
    pass_when: |
      test_cuda_ewald_coul still passes on the GPU node (Ewald and PME share the
      alpha real/reciprocal split convention).
    status: pending
out_of_scope:
  - ATV_CUDA_REAL default change; CPU PME/Ewald changes; relaxing the Madelung
    tolerance instead of fixing the reciprocal-space math.
---

# Acceptance — cuda-pme-accuracy

Done means: the CUDA PME reciprocal-space energy matches the analytic NaCl
Madelung reference to the same gate the CPU implementation already meets, with
the default float precision and without touching the tolerance. Ruled out before
this spec opened: unit renames (cancel in the ratio), float input precision
(fails identically under `-DATV_CUDA_REAL=double`), missing energy terms
(reciprocal + self are present; net-charge correction is zero for neutral NaCl),
and the Frame/scalars bridge (all other CUDA energy tests are green).
