# Architecture Blueprint

> 跑 `/mol:map` 填充本蓝图。
> Run `/mol:map` to populate this blueprint.

Consumed by `librarian` during `/mol:spec` Step 4.5 for reuse candidates and placement recommendations.

## Execution model

The CPU MD timestep (phase skeleton, per-phase responsibilities, component →
(phase, order) binding map, and the naming conventions) is documented in
**[`cpu-md-execution-model.md`](cpu-md-execution-model.md)** — the LAMMPS
*Developer flow* analogue. Start there to understand how one MD step runs.
