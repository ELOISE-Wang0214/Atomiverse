# Spec: `piston-electrode-pz`

| Field | Value |
| ----- | ----- |
| slug | `piston-electrode-pz` |
| branch | `wf/piston-electrode-pz` |
| domain | classical electrochemical MD (E6 electrode mechanics) |
| map slice | E6 — piston electrode + constant \(P_z\) (capability map §1 / §5) |
| method | Rigid electrode piston at target normal pressure \(P_z\) |
| geometry | Parallel plate: frozen lower wall + LJ fluid + movable upper wall |
| observable | Instantaneous and mean piston \(P_z\) [eV/Å³]; force balance residual |
| public API | `cpu::md::PistonElectrode` **Component** on `cpu::md::MD` |
| regression | `regressions/md/piston-electrode-pz/` |
| oracle | **metalwalls** (MW-only; rigid electrode piston ≠ generic barostat; no LAMMPS twin for electrode \(P_z\)) |
| queue | `docs/developer/metalwalls-method-queue.md` wave-2 `[~] piston-electrode-pz` |
| deps | E1 freeze + electrode forces (pair forces on wall atoms) |

## Goal

Land a **Component** that treats one electrode (atom type) as a **rigid
piston** free to translate only along \(z\), with an external normal force
\(F_{\mathrm{ext}} = \pm P_{\mathrm{target}} A\) so that equilibrium samples
constant normal pressure \(P_z\).

MetalWalls (JOSS 2020, 10.21105/joss.02373): *“constant pressure \(P_z\)
simulations … treating the electrodes as rigid pistons.”* Map row: *rigid
piston ≠ generic barostat*.

Architecture (binding):

- Component on existing `cpu::md::MD` — **not** a new Driver category.
- Parameters = constructor arguments / public members only — **no** public
  `*Config` bag, **no** `make_*` / `build_*` factory.
- Frame I/O only; selection by atom `type`.
- Does **not** solve electrode charges (const-pot is orthogonal; pair with
  `ConstPotBO` + `rebuild_matrix_every_step` when charges track movable walls).
- Does **not** live inside `KSpace*`.

## Algorithm (v1)

For electrode atoms of type \(t\):

1. **Area** \(A = L_x L_y\) from the current box (orthorhombic).
2. **Atom force** \(F_z^{\mathrm{atoms}} = \sum_{i\in t} F_{z,i}\) after the
   physical force stage (pair / bonded / long-range).
3. **External pressure force** (fluid between walls):
   - `Side::Upper` (wall above fluid): \(F_{\mathrm{ext}} = -P_{\mathrm{target}} A\)
   - `Side::Lower` (wall below fluid): \(F_{\mathrm{ext}} = +P_{\mathrm{target}} A\)
4. **Net** \(F_{\mathrm{net}} = F_z^{\mathrm{atoms}} + F_{\mathrm{ext}}\).
5. **Mass-weighted redistribution** so the electrode is a rigid \(z\)-body:
   \(F_{z,i} \leftarrow F_{\mathrm{net}}\, m_i / M\), \(F_{x,i}=F_{y,i}=0\),
   with \(M=\sum m_i\). Equalize \(v_z\) to COM velocity; zero \(v_x,v_y\).
6. **Instantaneous piston pressure** (diagnostic, positive compressive when
   fluid pushes the wall outward against \(P_{\mathrm{target}}\)):
   - Upper: \(P_z = F_z^{\mathrm{atoms}} / A\)
   - Lower: \(P_z = -F_z^{\mathrm{atoms}} / A\)

Bindings (default):

- `post_force @ +900` — after physical forces (+ biases), before `FixFreeze`
  (+1000) so freeze remains available for immobile walls; this piston
  **owns** its type and must not share that type with `FixFreeze`.
- `post_integrate @ 0` — restore rigid structure (equal \(v_z\), freeze \(xy\)
  relative offsets).
- `final_integrate @ +200` — re-equalize velocities after final half-kick.

Optional public member `tau_damp` [fs]: if `>0`, add COM friction
\(F_{\mathrm{damp}} = -M\, v_{\mathrm{com},z}/\tau\) into \(F_{\mathrm{net}}\)
for fast textbook equilibration (not part of the MW force identity).

## Non-goals

- Full electrode NPT with both walls + charge solve CG (map CG is
  `const-pot-cg-large-electrode`).
- Box \(L_z\) rescaling as the primary pressure control (use piston positions).
- Thomas–Fermi metallicity (E7).
- CUDA, SCF/PW/GPW, AIM/DAIM.
- Live MetalWalls / LAMMPS at ctest time.

## Tasks

1. Write this spec + `acceptance.md` (binding acceptance list).
2. Offline pin notes + `REFERENCE.md` with `oracle: metalwalls` and
   `oracle_version:` (algorithm/literature when no MW binary).
3. Implement `cpu::md::PistonElectrode`:
   - `include/atomiverse/cpu/dynamics/piston_electrode.h`
   - `src/cpu/dynamics/piston_electrode.cpp`
   - CMake `ATV_BUILD_MD` source list
4. Unit tests under `tests/cpu/dynamics/test_piston_electrode.cpp`.
5. E2E `regressions/md/piston-electrode-pz/` vs offline pin.
6. MD-only build + run unit + regression.
7. Perf gate: `perf_exempt` if no MW twin binary; else fair wall ≤ 10×.
8. Commit on `wf/piston-electrode-pz` only.

## Acceptance

See `acceptance.md` (binding). `acceptance_met=true` only when every criterion
is verified.
