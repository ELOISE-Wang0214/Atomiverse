# Constant-charge electrode control — parallel-plate conq twin

*(Technical-report draft; narrative house style TBD — rewrite later.)*

| Field | Value |
| ----- | ----- |
| Slug | `const-charge-electrode` |
| Regression | `regressions/md/const-charge-electrode/` |
| Public API | `cpu::md::ConstPotBO::add_electrode_charge` (Component on `cpu::md::MD`, default `pre_force`) |
| Roadmap / queue | `docs/developer/metalwalls-method-queue.md` → wave 5; map §1 "Constant-charge electrodes (control)" |
| Oracle | **LAMMPS** ELECTRODE package (`fix electrode/conq`; madelung conp log anchor) |

## Summary

**Constant-charge** electrode control is the dual of constant potential:
prescribe each electrode's **total charge** \(Q_e\), solve the per-site
distribution from the same Gaussian elastance matrix, and read the uniform
electrode potential \(\Psi_e\) as the **response** (LAMMPS ELECTRODE
`fix electrode/conq`, whose outputs `f_conq[e]` / `v_dv` map to
`electrode_potential(e)` / `delta_psi()` here). This is *not* "fixed uniform
charge + freeze": the total charge redistributes over the sites
self-consistently in the field of the electrolyte every step.

## Method

With indicator \(d_e\) of const-Q group \(e\) and unknown uniform group
potentials \(\psi_e\), the BO solve becomes the constrained system

\[
Q = A^{-1}\bigl(\Psi_{\text{base}} + \textstyle\sum_e \psi_e d_e - b\bigr),
\qquad d_e^{\mathsf T} Q = Q_e,
\]

reduced to the small SPD system \(M\psi = r\),
\(M_{ef} = d_e^{\mathsf T}A^{-1}d_f\),
\(r_e = Q_e - d_e^{\mathsf T}A^{-1}(\Psi_{\text{base}} - b)\), on both the
dense-Cholesky and CG paths. Registering any const-Q group supersedes the
`electroneutrality` projection (choose targets summing to
\(-q_{\text{elyt}}\) for global neutrality). Mixed const-\(\Psi\) /
const-\(Q\) groups are supported.

## Usage

```cpp
cpu::md::ConstPotBO conq(/*eta [1/Å]=*/2.0);
conq.add_electrode_charge(/*type=*/1, /*Q [e]=*/-0.5);
conq.add_electrode_charge(/*type=*/2, /*Q [e]=*/+0.5);
conq.xy_image_range = 250;
sim.add(std::move(conq));            // pre_force @ -50
// after a step:
// conq.total_charge(e)         — stays exactly the target ("Q fixed")
// conq.electrode_potential(e)  — the V response (LAMMPS f_conq[e])
// conq.delta_psi()             — LAMMPS v_dv
```

## Oracle contract

The LAMMPS anchor is the madelung conp log pin
(\(\Delta\Psi = 1\) V ⇔ \(Q_{\mathrm{top}} = +0.279214485147238\) e) plus the
package's documented inverse contract (`graph-il/in.conq`: *"symm on and off
give different electrode potentials, but identical potential difference"*):
conq fed the conp totals must return the conp \(\Delta\Psi\). The regression
gates that inverse to 1e−8 V, the charge scale to 1% vs the LAMMPS log, the
per-site distribution duality to 1e−10 e, and the empty-capacitor response
against the discrete Gaussian-lattice Madelung plate value
\(\Delta\Psi = 2k_e[2\pi\sigma d + \eta\sqrt{2/\pi} + M_\square + g(\eta)] =
651.8103\) V (1.5%; the bare-sheet \(4\pi k_e\sigma d\) omits the on-site
lattice term). Fair wall vs LAMMPS is **`perf_exempt`** (local `lmp`
lacks `PKG_ELECTRODE`).
