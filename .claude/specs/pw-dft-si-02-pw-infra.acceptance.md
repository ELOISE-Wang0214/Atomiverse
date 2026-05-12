---
slug: pw-dft-si-02-pw-infra
criteria:
  - id: ac-001
    summary: G-vector count matches reference for Si at 30 Ry
    type: code
    pass_when: |
      Si FCC cell a=10.26 Bohr, E_cut=15 Hartree → count ∈ [250, 300]; all G satisfy |G|²/2 ≤ E_cut.
    status: pending
  - id: ac-002
    summary: FFTW3 round-trip error < 1e-12
    type: code
    pass_when: |
      12×12×12 grid, random fill, R2C→C2R, max|diff| < 1e-12.
    status: pending
  - id: ac-003
    summary: UPF parser loads SG15 Si correctly
    type: code
    pass_when: |
      z_val()==4, l_max()>=1, mesh arrays non-empty, r monotonic.
    status: pending
  - id: ac-004
    summary: Ewald energy converges on cutoff doubling
    type: scientific
    pass_when: |
      |E_total(E_cut) - E_total(2*E_cut)| < 1e-6 Hartree for Si +4 charges.
    status: pending
  - id: ac-005
    summary: k-point weights sum to 1
    type: code
    pass_when: |
      4×4×4 MP → 64 points, Σweights=1.0 within 1e-14.
    status: pending
  - id: ac-006
    summary: Si crystal has correct atoms and positions
    type: code
    pass_when: |
      n_atoms()==2, type_id==14, positions at (0,0,0) and (a/4,a/4,a/4).
    status: pending
---
