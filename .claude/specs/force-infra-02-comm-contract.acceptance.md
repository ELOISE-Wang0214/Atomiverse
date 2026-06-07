---
slug: force-infra-02-comm-contract
criteria:
  - id: ac-001
    summary: ComponentBase carries comm_forward/comm_reverse + pack/unpack (CPU + CUDA)
    type: code
    pass_when: |
      Both cpu and cuda ComponentBase declare comm_forward/comm_reverse (default 0) and the
      four pack/unpack forward/reverse callbacks (default no-op), mirroring LAMMPS
      pair.h:50-52,207-210. Existing components compile unchanged.
    status: pending
  - id: ac-002
    summary: CPU + CUDA halo iterate participants; position/force/charge are built-ins
    type: code
    pass_when: |
      src/cpu/comm/halo_exchange.cpp and MultiGpu forward/reverse_comm no longer hardcode the
      field set; they iterate registered participants sized by sum of widths. Position-forward,
      force-reverse, and charge-forward are the built-in participants.
    status: pending
  - id: ac-003
    summary: CoulLong uses ghost charge from the contract; Allgatherv removed
    type: code
    pass_when: |
      cpu::CoulLong reads ghost charges supplied by the comm contract; the MPI_Allgatherv
      tag→charge map (tag_charge_/refresh_charge_map/use_charge_map_) is deleted.
    status: pending
  - id: ac-004
    summary: Synthetic extra field round-trips exactly (forward broadcast + reverse accumulate)
    type: code
    pass_when: |
      A test participant registers a 1-wide forward field and a 1-wide reverse field; forward
      sets ghost values to owners', reverse adds ghost contributions back into owners exactly
      once. Verified CPU (n8) + CUDA (multi-GPU).
    status: pending
  - id: ac-005
    summary: CoulLong real-space results unchanged via the new charge channel
    type: scientific
    evaluator_hint: operator-action — CPU mpiexec n1/n8
    pass_when: |
      CoulLong real-space coulomb energy + per-atom forces on n8 match n1 and the pre-change
      (Allgatherv) results within documented tolerance.
    status: pending
  - id: ac-006
    summary: Position/force halo bit-identical; suites green
    type: scientific
    evaluator_hint: operator-action — CPU n1/n8 + GH200 multi-GPU
    pass_when: |
      test_mpi_lj_energy_bruteforce (n1/n8) and test_cuda_multi_gpu (LJ/PE, 1/2/4-GPU) are
      bit-identical to baseline; full CPU + CUDA ctest green.
    status: pending
---

# Acceptance criteria

- **ac-001** — ComponentBase comm contract on CPU & CUDA (code).
- **ac-002** — halo iterates participants; pos/force/charge built-in (code).
- **ac-003** — CoulLong uses contract charge; Allgatherv deleted (code).
- **ac-004** — synthetic field forward/reverse round-trip exact (code).
- **ac-005** — CoulLong real-space unchanged via new charge channel (scientific; operator).
- **ac-006** — pos/force halo bit-identical; suites green (scientific; operator).
