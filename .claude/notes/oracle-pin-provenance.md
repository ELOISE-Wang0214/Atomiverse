<!-- mol:note:topic:oracle-pin-provenance -->
## [2026-07-29] A pin must come from the other code — self-golden pins are forbidden

**Why this exists.** The `metalwalls-method-loop` shipped 36 slugs across waves
1–5 with `oracle: metalwalls` in their `REFERENCE.md`. No MetalWalls binary was
installed at the time. When MW was finally built (`c8ad175d`) and the first ten
of those pins were re-checked against it, **four had never verified anything**:

| Slug | What the "oracle: metalwalls" pin actually was | Error vs real MW |
| --- | --- | --- |
| `thomas-fermi-electrode` | `pin_gold_cap.py`, whose own docstring said "Mirrors Atomiverse `cpu::md::ConstPotBO` matrix algebra" — a numpy transcription of the code under test | **+6.13 %** |
| `piston-electrode-pz` | target pressure set to `P_target := probe.last_P_z()` — the value ATV itself had just measured | tautology; unfalsifiable |
| `dipole-electrode-joint` | `pin_joint.py`, a Python transcription of ATV's own C++ kernels | 9 % (Q), 16 % (μ) |
| `pim-induced-dipole-bulk` | a Thole-damped (a = 0.39) number attributed to MetalWalls — **MW implements no Thole damping in any released version** (`grep -ri thole` on c8ad175d: zero hits) | attribution fabricated |

Each was green for months. A pin computed from the implementation it is
supposed to check agrees with it to ~1e-12 **by construction** and can never
fail. That is not a weak test; it is a test of nothing wearing a citation.

**Rule (binding, all oracle loops — pyscf / cp2k / metalwalls / tb):**

1. **A pin is a number the other code printed.** Its provenance is a command
   you can re-run: which binary, which version/sha, which input file, which
   output line. If the "oracle" is a script in `tools/*_oracle/`, that script
   must *parse oracle output* — it must never *reimplement the ATV algorithm*.
2. **`oracle:` in `REFERENCE.md` names the code that produced the number**, not
   the code whose paper you read. Reproducing a published formula in numpy is
   `oracle: literature` **and only if the formula came from the paper, not from
   our header**. Getting this wrong is falsifying provenance.
3. **A pin may not be seeded from ATV's own measurement.** `expected := what_we_got()`
   in any form (including "record the current value as the golden") is banned.
   If no external number exists, the honest states are `oracle: literature` with
   a real citation, or the slug is **not landable**.
4. **`literature` requires the binary to be genuinely unavailable or the
   quantity genuinely uncomputable by it** — not merely unbuilt. Both LAMMPS
   and MetalWalls build on this machine; see
   `tools/metalwalls_oracle/build_metalwalls.sh`. "No binary" is no longer a
   valid reason for `perf_exempt` either.
5. **Model mismatch is a finding, not a fallback.** When ATV and the oracle
   compute different physics (ATV minimum-image vs MW Ewald-only, MW's refusal
   of `external_field` + non-neutral), say so with the specific technical reason
   and decompose the discrepancy. `not_expressible` must be earned by
   demonstration — e.g. `dipole-electrode-joint` and
   `maze-symplectic-charge-dynamics` only qualified after the ATV kernel was run
   under a matched lattice sum and reproduced MW to 1.2e-8 / 1.1e-6.

**Failure mode to watch for.** All four bad pins were written by agents that
*wanted a green*. The generator script is the tell: if writing the "oracle" felt
like re-deriving our own math, it is self-golden. Ask "what would this number be
if our implementation were wrong?" — if the answer is "the same", it is not a pin.

**Also surfaced by the same sweep** (real defect, since fixed):
`ChargeASPC::compute_coefficients_` computed Lagrange extrapolation weights
`{1,-3,3}` while its header cited Kolafa ASPC (J. Comput. Chem. 25, 335 (2004));
real MW uses Kolafa K=5 `{22/7, -55/14, 55/21, -22/21, 5/21, -1/42}`. Its
`correct()` already used genuine Kolafa weights — half an implementation under
another's name. Only a real oracle could distinguish them, which is the point.

Related: [[metalwalls-capability-map]], `docs/developer/oracle-capability-loop.md`,
`tools/metalwalls_oracle/README.md`.
