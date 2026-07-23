# Spec — electrode-charge-predictor (E2)

## Goal

Ship an **ASPC-like electrode charge predictor** for constant-potential MD so
iterative (CG) electrode charge solves need **fewer iterations** than a cold
start (zero / no history). MetalWalls map **E2**. Component surface on
`cpu::md::MD` — **not** a new Driver category; **not** inside `KSpace*`.

## Rationale

MetalWalls (and Kolafa ASPC for smooth MD fields) extrapolates electrode
charges \(Q(t)\) from prior steps to warm-start the BO/CG charge solve when the
electrolyte configuration moves. Dense Cholesky always costs \(O(N^3)\) once;
the predictor pays off on the **iterative** path that large electrodes use.

## Public API (minimal)

| Symbol | Role |
| --- | --- |
| `cpu::electrodes::ChargeASPC` | Kolafa ASPC for charge **vectors** (order \(K\), record / predict / correct / reset). No `*Config`, no factories. |
| `cpu::md::ConstPotBO` | Existing BO Component gains optional **`use_cg`**, **`predictor_order`**, **`cg_tol`**, **`max_cg_iters`**, **`last_iters()`**. When `use_cg` and `predictor_order > 0`, ASPC warm-starts CG from prior \(Q\). |

Parameters are ctor args / public members only.

## Architecture constraints

- Frame I/O only (`kAtom_charge` write on electrode sites).
- No public `*Config` bags; no `make_*` / `build_*` factories.
- Charge solve stays on the Component — never welded into `KSpace*`.
- Units: eV, Å, e, V on the public surface.
- MD-only build (HF/PW/GPW/CUDA/Torch OFF).

## Oracle

| Field | Value |
| --- | --- |
| oracle | **metalwalls** (MW-only predictor path; LAMMPS ELECTRODE has no ASPC twin) |
| offline pin | `regressions/md/electrode-charge-predictor/REFERENCE.md` |
| live third-party in ctest | **forbidden** |

## Observable

On a short parallel-plate const-pot toy with moving electrolyte:

1. CG from **cold start** (zero \(Q\), no ASPC) iteration count \(N_{\mathrm{cold}}\).
2. CG with **ASPC order 2** after history filled: \(N_{\mathrm{warm}}\).
3. Acceptance: \(N_{\mathrm{warm}} < N_{\mathrm{cold}}\) (mean over post-fill steps).
4. CG \(Q\) matches dense Cholesky \(Q\) within stated tol on the same RHS.

## Out of scope

Large-electrode CG preconditioning (sibling `const-pot-cg-large-electrode`),
piston electrodes, Thomas–Fermi, SCF/HF ASPC density, CUDA.
