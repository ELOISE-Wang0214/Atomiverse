# Specs index

| slug | status | summary |
|---|---|---|
| [gpw-gaussian-plane-wave](gpw-gaussian-plane-wave.md) | draft | GPW DFT (Gaussian basis + uniform-grid Hartree/XC + analytic GTH PP) targeting CP2K Quickstep parity; goal = getting-started tutorial passes, then migrate tutorials/benchmarks into `benchmarks/gpw/`. Depends on native-xc-libxc-parity for XC. |
| [native-xc-libxc-parity](native-xc-libxc-parity.md) | draft | Native `cpu/xc/` functional library reproducing libxc bit-for-bit (PBE etc., validated vs libxc.so oracle) and CP2K-native PADE (validated vs CP2K). Powers both `cpu::scf` and `cpu::gpw`. |
| [gpw-edl-aimd](gpw-edl-aimd.md) | approved | GPW production bridge to run CP2K's ZnO electric-double-layer AIMD (1214-atom periodic NVT, Zn/O/H/Na/Cl): ExtXYZ-via-molrs load, Zn/Na/Cl GTH+MOLOPT validation, analytic GPW forces, CSVR thermostat, large-system scale/perf, DFT-D3. Builds on gpw-gaussian-plane-wave. |
