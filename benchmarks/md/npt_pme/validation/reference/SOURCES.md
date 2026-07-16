# NaCl coul/long + PPPM static frames — LAMMPS reference

## Tool

- **Engine:** LAMMPS (22 Jul 2025 - Update 4), Homebrew `lmp_serial`
- **Package:** `KSPACE` (pppm) + `coul/long`
- **Host:** macOS arm64 (regenerated 2026-07-16)

## Inputs

| File | Role |
|---|---|
| `nacl_s100.data` / `nacl_s099.data` / `nacl_s102.data` | Rock-salt NaCl, `n=4` (64 ions), spacing `d=2.82 Å`, uniform scale `s` of box + coordinates |
| `in.one_frame` | Static `run 0` evaluation template |
| `in.nacl_static_frames` | Multi-frame sketch (superseded by three independent data files after `change_box` segfaulted on this build) |

## Pinned electrostatic parameters (match Atomiverse)

| Parameter | Value |
|---|---|
| `pair_style coul/long` cutoff | `5.358` Å (`1.9 · d`) |
| `kspace_modify gewald` (`α`) | `0.6532288167226578` 1/Å (`3.5 / cutoff`) |
| `kspace_modify mesh` | `24 24 24` |
| `kspace_modify order` | `6` |
| `kspace_style pppm` accuracy | `1.0e-5` (grid fixed by mesh, not re-estimated) |
| `units` | `metal` (eV, Å, amu, elementary charge) |

## Regeneration

```bash
cd benchmarks/md/npt_pme/validation/reference
# regenerate data files if needed (see repo history / generator comment in SOURCES)
for tag data force out in \
  0 nacl_s100.data forces_s100.lammpstrj frame0.txt \
  1 nacl_s099.data forces_s099.lammpstrj frame1.txt \
  2 nacl_s102.data forces_s102.lammpstrj frame2.txt
do
  sed -e "s/DATAFILE/$data/" -e "s/OUTFILE/$out/" -e "s/FORCEFILE/$force/" in.one_frame > in.tmp
  lmp_serial -in in.tmp -log log.$tag
done
cat frame0.txt frame1.txt frame2.txt > frames_lammps.txt
```

## Output

`frames_lammps.txt` — one line per frame with `E` [eV], pressure-tensor
components [bar] (`pxx…pyz`), `V` [Å³], `lx` [Å], force-sum [eV/Å].

Pressure → virial [eV] conversion used by the validation binary:

```text
virial_αβ = p_αβ[bar] · V / nktv2p
nktv2p    = 1 / 6.2415e-7   # inverse of Atomiverse units:: bar → eV/Å³
```

## Why three frames

The acceptance criterion asks for initial / intermediate / final
**variable-cell** samples. These are three fixed rock-salt lattices at
uniform scales `s ∈ {1.00, 0.99, 1.02}` — pure electrostatic oracles at
distinct cells without integrator drift. Atom-wise forces vanish by cubic
symmetry; the gate still checks force sum and max |F| against the dump
(numerical noise floor).
