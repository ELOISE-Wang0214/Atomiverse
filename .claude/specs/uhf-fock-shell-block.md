# UHF Fock shell-block port (`build_jk_uhf_direct` → block)

**Status:** draft (review-derived implementation plan, 2026-07-02).
**Origin:** multi-axis CPU review — perf axis 🔴 `operators.cpp:185-282`: the UHF
Fock builder never got the shell-block treatment the RHF twin documents as
~150× (operators.cpp:336-338). Every UHF / UKS-hybrid-exchange SCF iteration
pays one `native::eri()` call per AO quartet + xtensor `operator()` accumulation.
**Gate:** implement on the cluster where the SCF benchmarks run; the 1e-9
elementwise AO-vs-block parity test is the correctness gate (solver-level
tolerances are too loose to catch degeneracy-factor mistakes).

## Context (verified 2026-07-02)

- Slow path: `build_jk_uhf_direct`, `src/cpu/scf/operators.cpp:192-289` (one
  `native::eri()` per quartet at :225-226; thread-local Jv/Kav/Kbv at :202-204).
- Fast twin: `build_g_rhf_direct_block` (operators.cpp:500-670) with
  `raw_sym_add` (:341-346). A K-only block twin already exists:
  `build_k_direct_block` (:348-498) — **no 0.5 factor, exactly the UHF K
  normalization** (comment at :407).
- ERI engine: `native::eri_block()` (native_integrals.h:200-201) — (ss|ss) →
  xsimd kernel, else Rys. Block descriptors `native::ShellBlock`
  (native_integrals.h:47-55), `ShellBlockPairList` (shell_pair.h:62-94).
  Buffer bound `kMaxNCart = 10` → `blk[10^4]` = 80 KB stack (proven at :535).
- Production call sites (only two): `internal/hamiltonian_builder.cpp:73`
  (UHF `F_σ = Hcore + J − K_σ`) and `:158` (unrestricted-KS hybrid
  `F_σ = Hcore + J − α_hf·K_σ + Vxc_σ`). `prepare()` already builds
  `sp_block_` unconditionally (:20) — no plumbing needed.
- Adjacent, same disease: `build_j_direct` (:291-331) is still AO-level and is
  the production Coulomb builder for **all pure-KS** runs
  (hamiltonian_builder.cpp:127, 169). `build_j_direct_block` is a trivial
  subset of this port — do it in the same spec as a second commit.

## Plan

New `build_jk_uhf_direct_block(const ShellBlockPairList& sp, P_total, P_alpha,
P_beta, J, K_alpha, K_beta)` in operators.cpp (decl in operators.h next to
:70-73). **One fused ERI pass** — do NOT compose from
`build_g_rhf_direct_block` + 2× `build_k_direct_block` (that is 3-4 ERI passes).

Construct mapping (RHF-block → UHF-block):

| RHF-block construct (operators.cpp) | UHF-block equivalent |
|---|---|
| `Graw(nT, vector<double>(n*n))` thread arena (:509) | Three arenas: `Jraw`, `Karaw`, `Kbraw` |
| `const double* Pp = P.data()` (:510) | `Ptp` / `Pap` / `Pbp` (column-major, symmetric) |
| `omp parallel for schedule(dynamic)` over canonical block pairs, `right <= left` (:513-522) | identical |
| Block Schwarz `bp_l.schwarz * bp_r.schwarz < thr` (:526) | identical (+ keep per-element `abs(val) < thr` skip :557/:621) |
| `double blk[kMaxNCart^4] = {}; eri_block(...)` (:535-536) | identical |
| `left != right`: `diag_bra`/`diag_ket` AO pruning (:538-540) | identical |
| `left == right`: `ao_pairs[36]` sub-canonicalization (:592-616) | identical, verbatim |
| Coulomb `raw_sym_add(Gp, a,b, fcd·P(c,d)·val)` (:563-566, :626-628) | same into `Jraw`, reading `Ptp` |
| Exchange `hval = 0.5*val; Gp[..] -= hval·P(..)` (:569-587) | `Kraw[..] += val·Pσ(..)` — **sign +, no 0.5** — duplicated for α (`Pap`) and β (`Pbp`); line-for-line the scatter of `build_k_direct_block` (:407-425, :460-485) |
| Diagonal-quartet exchange: 3 `raw_sym_add` with `c!=d && a!=b` guard and the omitted `K(a,d)` duplicate (:651-655; rationale :147-182) | copy verbatim per spin — do not re-derive |
| Final reduce: flat loop per thread arena (:663-669) | same ×3 |

Then swap the two hamiltonian_builder call sites (:73, :158) — identical
signature, one token each. Keep the AO-level `build_jk_uhf_direct` as the
parity-test oracle (mirroring `build_g_rhf_direct` / `compare_g_ao_vs_block`).

## Accumulation formulas (8-fold symmetry; canonical a≥b, c≥d, (ab)≥(cd); val = (ab|cd); P symmetric)

Off-diagonal quartet ((ab) ≠ (cd)), `f_ab = (a==b?1:2)`, `f_cd = (c==d?1:2)`:

```
J(a,b) += f_cd·P_tot(c,d)·val ;  J(b,a) same if a≠b      [sym_add]
J(c,d) += f_ab·P_tot(a,b)·val ;  J(d,c) same if c≠d      [sym_add]
for σ ∈ {α, β}:
  Kσ(a,c) += Pσ(b,d)·val ;  Kσ(c,a) += Pσ(b,d)·val        [explicit double-add, ALSO when a==c]
  if a≠b:            Kσ(b,c) += Pσ(a,d)·val ;  Kσ(c,b) += Pσ(a,d)·val
  if c≠d:            Kσ(a,d) += Pσ(b,c)·val ;  Kσ(d,a) += Pσ(b,c)·val
  if c≠d and a≠b:    Kσ(b,d) += Pσ(a,c)·val ;  Kσ(d,b) += Pσ(a,c)·val
```

Diagonal quartet ((ab) == (cd), i.e. a==c, b==d):

```
J(a,b) += f_cd·P_tot(c,d)·val (+ mirror)                  [single J term]
for σ ∈ {α, β}:
  Kσ(a,c) += Pσ(b,d)·val               (+ mirror if a≠c)  [sym_add]
  if a≠b:          Kσ(b,c) += Pσ(a,d)·val (+ mirror)
  if c≠d and a≠b:  Kσ(b,d) += Pσ(a,c)·val (+ mirror)
  — Kσ(a,d) intentionally OMITTED (duplicates the Kσ(b,c) mirror; see
    operators.cpp:147-155, :171-182)
```

Critical subtlety (operators.cpp:147-155): the K double-add is explicit, not
`sym_add` — when the target is diagonal (e.g. a==c) it must still be added
twice (the mirror permutation is a distinct ERI permutation). Already correct
in the AO UHF code (:241-261) and `build_k_direct_block` — copy, don't reinvent.

## Correctness pins

- `tests/cpu/scf/test_uhf_solve.cpp`: H −0.4666±5e-3 (:86); He UHF==RHF at
  1e-3 (:23, :108-116); Li ∈ (−8,−7) (:146-153); H2 PES UHF ≤ RHF+1e-5 (:196-207);
  spin-contamination gates 1e-6/1e-4/1e-3.
- `tests/cpu/scf/test_uks_solve.cpp` pins PySCF LSDA at 1e-6..1e-4 Ha but
  routes through `build_j_direct` (pure KS) — only gates this port if the
  `build_j_direct_block` follow-up is folded in.
- **Coverage hole:** the unrestricted-hybrid site (hamiltonian_builder.cpp:153-165)
  has NO test (no open-shell B3LYP/PBE0 anywhere). Add an open-shell hybrid
  smoke gate or stage the `:158` swap as a separate commit.
- **Required new test:** `compare_jk_ao_vs_block` in
  `tests/cpu/scf/test_integrals.cpp` (model: `compare_g_ao_vs_block` :375-407,
  elementwise 1e-9), H2O/CH4/H2 with a **dense random symmetric
  P_alpha ≠ P_beta** — the existing tridiagonal test density is too weak for
  spin discrimination.

## Size / risks

~330 insertions / ~6 modified lines (operators.cpp +~210, operators.h +~45,
hamiltonian_builder.cpp ~4, test +~70). Additive, revertible per call site.

1. Degeneracy/branch structure (highest): a mistake shifts E by ~1e-3 Ha and
   can pass the loose solver tests — the 1e-9 parity test is the real gate.
2. Spin-asymmetry copy-paste (`Pap` in a β line) is invisible on closed-shell
   fixtures — parity test must use P_α ≠ P_β; doublet tests (H, Li) stay green.
3. Precondition: both Pσ symmetric (true by construction; document in Doxygen).
4. Thread-local memory ×3: def2-TZVP n≈500, 64 threads → ~380 MB. Acceptable;
   note for a future blocked reduction.
5. Expected win: same mechanism as the ~150× RHF speedup (Rys-table
   amortization + raw pointers); UHF scatters slightly more per quartet.
