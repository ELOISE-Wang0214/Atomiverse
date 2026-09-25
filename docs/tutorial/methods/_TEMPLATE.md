# <Method name> — <system / basis>

*(Technical-report draft; narrative house style TBD — rewrite later.)*

| Field | Value |
| ----- | ----- |
| Slug | `<slug>` |
| Public API | `…` |
| Roadmap | link / row name |

## Summary

What the method computes, typical use, and where it sits relative to RHF/RKS/MP2/….
Three to six sentences. No pedagogy.

## Method

Equations, spin treatment, kernel/response level, important approximations.
Use display math where useful. State what is **not** included.

## Implementation

Public surface only (`cpu::Frame`, concrete method members; no `*Config` bags).

```cpp
// minimal call sequence — match the capability case
```

```bash
cmake --build build --target reg_<slug>
reg_<slug>
```

Name the public types readers call; do not dump repository paths.

## Validation (vs PySCF)

Source: case `<slug>` pin (offline pin; never live in ctest).

| Quantity | Atomiverse / pin | PySCF | Tol | Notes |
| -------- | ---------------: | ----: | --- | ----- |
| … | … | … | … | … |

Geometry, basis, `cart`/spherical, `conv_tol`, PySCF version — all from the pin notes.

## Performance

```bash
# fair wall vs the oracle twin for this case
reg_<slug>   # or the matching fair-bench invocation for the domain
```

| Twin | \(T_{\mathrm{PySCF}}\) [s] | \(T_{\mathrm{ATV}}\) [s] | \(r\) | ≤10× | Caveat |
| ---- | -------------------------: | -----------------------: | ----: | ---- | ------ |
| … | … | … | … | ✓/✗ | … |

Gate: \(r = T_{\mathrm{ATV}}/T_{\mathrm{PySCF}} \le 10\) unless exempted.

## Limitations

Scope, known failures, speed debt, follow-ons.

## References

- Offline pin: case `<slug>`
- Literature …
