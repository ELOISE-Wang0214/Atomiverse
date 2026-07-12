# Reference — RATTLE rigid water: rigidity and equipartition

## What this case pins

Two things, and the second is the one that matters.

### 1. The constraint invariant

A constraint is not judged by a number it reproduces; it is judged by an **invariant it
never breaks**. Every constrained distance — the two O–H bonds and the H–H distance
that encodes the H–O–H angle — must stay at its target for the whole trajectory.

Measured: worst relative deviation **5 × 10⁻¹¹** over 5 ps at `dt = 2 fs`. (A *flexible*
water could not be integrated at 2 fs at all: the O–H stretch has a period of ~9 fs.
Removing that vibration is the entire practical point of rigid water.)

Geometry: TIP3P (Jorgensen et al., *J. Chem. Phys.* **79**, 926 (1983)) —
r(OH) = 0.9572 Å, θ(HOH) = 104.52°, hence d(HH) = 2 r sin(θ/2) = 1.5139 Å.

### 2. Equipartition — which is really a test of the degree-of-freedom count

The **equipartition theorem**: at equilibrium each quadratic degree of freedom carries
`½ k_B T`, so

$$
T = \frac{2\,\text{KE}}{N_\text{dof}\,k_B}.
$$

(Allen & Tildesley, *Computer Simulation of Liquids*, 2nd ed., §2.2; Frenkel & Smit,
§4.3.)

A **rigid** water molecule has **6** degrees of freedom — 3 translational, 3 rotational
— not 9. The three constraints remove the three vibrational modes. For N molecules with
centre-of-mass motion removed:

```
Ndof = 6 N_mol - 3        NOT   3 N_atoms - 3 = 9 N_mol - 3
```

So the assertion is: thermostat the box at 300 K, and the temperature computed with the
**correct** Ndof must come back at 300 K.

Measured, 125 rigid waters, CSVR thermostat at 300 K, averaged over 6 ps:

```
Ndof (rigid, 6/mol) = 747      [3N-3 would be 1122]
<T>                 = 295.97 K   vs  300 K setpoint   (1.3 %, tol 3 %)
<T> with naive Ndof = 197.05 K   <- 33 % cold
```

**That last line is the bug this case exists to prevent.** Before this work, five
separate sites in the engine hardcoded `Ndof = 3N − 3`. With rigid water that ratio is
`(6N−3)/(9N−3) ≈ 2/3`, so the engine would have reported the water as ~33 % colder than
it was — and worse, since a thermostat's target kinetic energy is `Ndof k_B T / 2`, it
would have *driven the system to the wrong temperature* and reported success. Nothing
crashes. The trajectory is simply wrong.

The 3 % tolerance is the statistical spread of a thermostatted mean over 6 ps of a
125-molecule box, not a fudge factor: a broken DOF count misses by 33 %, an order of
magnitude outside it.

## What this case is NOT

**It is not a TIP3P water simulation, and it cannot be one today.** Two engine
limitations force that, and they are worth stating plainly rather than papering over:

1. **`PairLJCut` carries a single ε/σ for all atoms** — there is no per-type mixing
   table. In TIP3P the Lennard-Jones site sits on oxygen *only*; here the hydrogens are
   forced to wear the oxygen's σ = 3.15 Å too. At liquid density (O–O spacing 3.1 Å)
   neighbouring molecules then overlap catastrophically, and no constraint solver can
   hold a molecule that is being blown apart. The box here is therefore **dilute**
   (4.6 Å spacing), which equipartition does not care about — it holds at any density.

2. **No charges.** TIP3P is mostly electrostatics; the engine's PME exists but is not
   wired into this case. Without charges there are no hydrogen bonds and this is not
   liquid water — it is a fluid of rigid rotors.

Neither limitation weakens what the case actually pins. Rigidity is a property of the
constraint solver, and equipartition is a property of the DOF bookkeeping; both are
independent of how good the force field is. A real water benchmark belongs in
`benchmarks/`, against an external engine, once per-type LJ and charges are wired.

## Cross-references

The constraint algebra itself (RATTLE: position stage, velocity stage, iteration to
tolerance) is unit-tested in `tests/cpu/md/test_fix_rattle.cpp`, including the second
invariant this case does not check — that the relative velocity along every constrained
bond vanishes, `d|r|²/dt = 2 r·(v_i − v_j) = 0`. Skip the velocity stage and the bond
*lengths* stay right while the molecule quietly absorbs energy along its own bonds.

- H. C. Andersen, "RATTLE: A 'velocity' version of the SHAKE algorithm for molecular
  dynamics calculations", *J. Comput. Phys.* **52**, 24 (1983).
- J.-P. Ryckaert, G. Ciccotti, H. J. C. Berendsen, "Numerical integration of the
  Cartesian equations of motion of a system with constraints", *J. Comput. Phys.* **23**,
  327 (1977) — the original SHAKE.
