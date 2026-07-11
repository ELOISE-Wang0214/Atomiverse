# Specs index

| slug | status | summary |
|---|---|---|
| [gpw-edl-aimd](gpw-edl-aimd.md) | approved | GPW production bridge → 1214-atom ZnO electric-double-layer AIMD. **Refactored 2026-07-11**: the blocker was structural, not physical — every GPW test that touched a real molecule read its PP/basis data from a CP2K install at test time (`ATV_CP2K_DATA_DIR`, 14 tests), so with the gate off GPW had 12 pure-math tests and zero force validation. New M0 makes GTH-PBE + DZVP-MOLOPT-SR-GTH in-tree constexpr tables (published in Krack 2005 / VandeVondele-Hutter 2007, mirroring `cpu/scf/basis_data/`); the four CP2K-parity criteria are replaced by reference-free gates (FD force oracle, energy conservation, CSVR H̃ drift) plus pinned anchor constants. Remaining physics: M5 `gth_matrix_element` atom screening (O(N³) → feasible at N=1214), then the capstone. |
