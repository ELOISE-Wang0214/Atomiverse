# Offline pin — GPW DFT-D4 + ATM (`gpw-d4-atm`)

## Geometry

Same periodic water dimer as `gpw-d4-water` (`POSCAR` here is a copy).

## Hard pins (dftd4 3.7.0)

```bash
export PATH=/software/sse2/el9_epyc9005/easybuild/prefix/software/dftd4/3.7.0-gfbf-2025b/bin:$PATH
export LD_LIBRARY_PATH=/software/sse2/el9_epyc9005/easybuild/prefix/software/dftd4/3.7.0-gfbf-2025b/lib:$LD_LIBRARY_PATH
# + FlexiBLAS / GCCcore from module load buildenv-gcc/2026.03-mpich

dftd4 -f pbe --json pin_periodic_atm.json POSCAR
# s9=1 (default pbe-D4-ATM)
# Dispersion energy: -9.9554475834240E-04 Eh

dftd4 -f pbe --mbdscale 0.0 --json pin_periodic_s9_0.json POSCAR
# two-body only: -9.9458842745847E-04 Eh
```

Arrhenius **CP2K 2026.1** lacks DFTD4 — no QS ENERGY twin on that host.

## Files

| File | Content |
| ---- | ------- |
| `POSCAR` | VASP-format water dimer (12 Å cubic) |
| `pin_periodic_atm.json` | dftd4 3.7.0 energy with s9=1 |
| `pin_periodic_s9_0.json` | two-body context |
| `pin_isolated_atm.json` | isolated s9=1 energy |
