---
slug: pw-dft-si-01-core-types
criteria:
  - id: ac-001
    summary: Mat3 inverse × Mat3 == identity
    type: code
    pass_when: |
      For 5 random non-singular 3x3 matrices, frobenius_norm(I - M * M.inverse()) < 1e-8.
    status: verified
    last_checked: 2026-05-12
  - id: ac-002
    summary: Mat3 determinant matches cofactor expansion
    type: code
    pass_when: |
      For 5 random matrices, |M.det() - cofactor_ref(M)| < 1e-8.
    status: verified
    last_checked: 2026-05-12
  - id: ac-003
    summary: Mat3 × Vec3 matches reference
    type: code
    pass_when: |
      For 5 random (M, v) pairs, each component of M*v matches Σ_j M(i,j)*v(j) within 1e-8.
    status: verified
    last_checked: 2026-05-12
  - id: ac-004
    summary: Vec3 dot/cross/norm match analytic
    type: code
    pass_when: |
      For 5 random (u,v): u.dot(v)==Σu_i*v_i, u.cross(v) matches formula, u.norm()==√u.dot(u).
    status: verified
    last_checked: 2026-05-12
  - id: ac-005
    summary: Vec3 normalized() == 1
    type: code
    pass_when: |
      For 5 random non-zero vectors, |v.normalized().norm() - 1.0| < 1e-12.
    status: verified
    last_checked: 2026-05-12
  - id: ac-006
    summary: Box to_frac/to_cart roundtrip
    type: code
    pass_when: |
      For 5 random (cell, point) triples, to_cart(to_frac(point)) == point within 1e-10.
    status: verified
    last_checked: 2026-05-12
  - id: ac-007
    summary: Box wrap folds into [0, L)
    type: code
    pass_when: |
      For 100 random points in a cubic cell, each fractional component of wrapped point is in [0, 1).
    status: verified
    last_checked: 2026-05-12
  - id: ac-008
    summary: Box shortest_vector within half-box
    type: code
    pass_when: |
      For 100 random (p1, p2, cell) triples, each fractional component of shortest_vector is in [-0.5, 0.5].
    status: verified
    last_checked: 2026-05-12
  - id: ac-009
    summary: Box recip_cell × cell == 2π·I
    type: code
    pass_when: |
      For 3 random cells, frobenius_norm(recip_cell() * cell - 2π*identity()) < 1e-10.
    status: verified
    last_checked: 2026-05-12
  - id: ac-010
    summary: kAngstrom2Bohr converts correctly
    type: code
    pass_when: |
      5.43 * kAngstrom2Bohr ≈ 10.26 within 1e-6 relative.
    status: verified
    last_checked: 2026-05-12
  - id: ac-011
    summary: kRy2Hartree == 0.5, kRy2eV == 13.605693
    type: code
    pass_when: |
      Exact compile-time constexpr match.
    status: verified
    last_checked: 2026-05-12
  - id: ac-012
    summary: All existing Mat3 call sites compile with rewritten interface
    type: code
    pass_when: |
      Full project builds successfully with new Mat3/Vec3/Box. No old-API usage remains.
    status: verified
    last_checked: 2026-05-12
---
