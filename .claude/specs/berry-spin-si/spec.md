# berry-spin-si — collinear-spin King–Smith–Vanderbilt polarization (PW)

**Loop:** cp2k-method-loop, wave 10. **Domain:** response. **Backend:** CPU PW.
**Depends on:** `berry-polarization` (closed-shell KSV strings) — done.
**Sibling residual (out of scope):** metallic Berry (smearing), GPW Berry, CUDA.

## What

Lift the `BerryPolarization::compute` throw on `method.spin_polarized` and
evaluate **collinear-spin** Berry-phase formal polarization on the CPU PW
method. Each spin channel \(s\in\{\uparrow,\downarrow\}\) has its own occupied
manifold of size \(N_s\); string phases \(\varphi_{i,s}\) are averaged over the
transverse mesh independently, and the electronic fractional polarization is

\[
p^\mathrm{el}_i = \sum_{s=\uparrow,\downarrow}
  \frac{\varphi_{i,s}}{2\pi}
  \qquad\text{(spin factor \(f=1\) per channel)}.
\]

Ionic bookkeeping is spin-independent (\(p^\mathrm{ion}_i =
\sum_a Z^\mathrm{val}_a \tau_{a,i}\)). Total formal polarization wraps with
quantum \(1\) (one electron per lattice vector; nspin=2 convention — cf.
closed-shell quantum \(2\) with \(f=2\)).

When \(N_\uparrow = N_\downarrow\) (multiplicity 1, even \(N_e\)) the physics
must reproduce the closed-shell pins of `berry-polarization` on diamond Si
(modulo branch choice of the smaller quantum).

Literature: King-Smith & Vanderbilt, *Phys. Rev. B* **47**, 1651 (1993);
Resta, *Rev. Mod. Phys.* **66**, 899 (1994); Spaldin, *J. Solid State Chem.*
**195**, 2 (2012); Pick–Cohen–Martin ASR, *Phys. Rev. B* **1**, 910 (1970).

CP2K twin: same modern-theory-of-polarization observable as
`si_berry_polarization_note.inp` (`&MOMENTS PERIODIC T`); absolute phases
are basis-dependent (GTH/GPW vs ONCV PW); **symmetry pins are the product**.

## Public surface (architecture-conform)

Extend existing `cpu::pw::BerryPolarization` — **no** new class, **no**
`*Config`, **no** factories:

```cpp
cpu::pw::PW method;
method.spin_polarized = true;
method.multiplicity = 1;   // N↑−N↓ = 0 for Si; or 2S+1 for open shell
cpu::pw::BerryPolarization berry;
auto r = berry.compute(method, frame);
// r.phase_el  — mean of spin phases (so p_el = 2*phase_el/2π still holds
//               when N↑=N↓; when N↑≠N↓, p_el = (φ↑+φ↓)/2π exactly)
// r.p_el_frac / r.p_ion_frac / r.p_frac / r.dipole_e_ang
```

Optional diagnostic members on the result (spin-resolved phases) are fine
if needed for tests; keep them public fields on `BerryPolarizationResult`,
not a bag type.

## Internal work

1. Shared **collinear fixed-occupation SCF** for post-processors:
   `internal::run_collinear_scf` next to `run_restricted_scf` (second
   non-DFPT consumer of the collinear loop shape; DFPT keeps its specialized
   response path). Leaves frozen total \(\rho\), spin XC, and per-spin
   \(V_\mathrm{eff}(r)\) on the problem for NSCF.
2. `BerryPolarization::compute`:
   - Reject metallic smearing (unchanged).
   - If `!spin_polarized`: existing restricted path (unchanged).
   - If `spin_polarized`: collinear SCF → for each requested direction and
     each spin with \(N_s > 0\), NSCF Davidson on strings with frozen
     \(V_{\mathrm{eff},s}\), occupied overlap determinants, string phase
     average; combine as above. Empty spin channel contributes phase 0.
3. Multiplicity / parity checks mirror `PW::run` / DFPT spin (throw on
   incompatible \(N_e\) and \(2S+1\)).

## Case

Bulk diamond Si (same geometry/settings as `berry-polarization`):

| Setting | Value |
| ------- | ----- |
| Structure | diamond Si, 2-atom fcc, \(a=5.43\) Å |
| PP | `tests/data/pw/Si_lda_nc.upf` |
| ecut / k_mesh / XC | 15 Ry / 2 / LSDA |
| Spin | `spin_polarized=true`, `multiplicity=1` (\(N_\uparrow=N_\downarrow=4\)) |
| Strings | `k_perp=2`, `k_string=4` |

Pins (literature, basis-independent):

1. Centrosymmetric formal polarization — every total \(p_i\) within 0.05 of
   the nearest multiple of 0.5 (Spaldin 2012).
2. Born effective charge \(Z^*=0\) for elemental crystal under sublattice
   displacement \(\delta=0.02\) along \(\mathbf a_3\) (\(|Z^*|<0.5\);
   Pick–Cohen–Martin 1970).
3. Ionic bookkeeping exact (same as closed-shell).
4. SCF converges on both geometries.

## Out of scope

- Metallic (smeared) Berry strings
- Noncollinear spin
- GPW / Gaussian Berry
- CUDA
- Wannier / MLWF spin (`periodic-wannier-uks` sibling)
- Magnetic open-shell crystal product pin (multiplicity>1) — optional
  diagnostic only if free; Si multiplicity-1 is the acceptance gate
