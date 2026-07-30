# Frame column access: `span` vs `get` (and Python)

**Status:** binding (2026-07-29). Router: `CLAUDE.md` iron law
`frame-column-access`.

## Decision

Operator call (discussion → “你说了算”): keep **C++** Frame column access as
**`span`**. Do not rebrand the engine API as `get_column` for “pythonic”
feel. Friendlier table-style names, if any, land on the **Python binding**.

## Why `span` on C++

1. **Return type is a view.** `Span<T>` is non-owning continuous storage over
   a column — same family as C++20 `std::span`. The name documents “range,
   not a copied value.”
2. **Pairs with `meta.get`.** Meta is “one typed result”; columns are “index
   over all atoms.” Two verbs prevent confusing SoA columns with scalar meta.
3. **Hot-path contract.** Components resolve once at launch entry and iterate
   pointers; resolution-budget / generation checks attach to that resolve.
   “Get column” invites map-value mental models and re-resolve-in-loop smells.
4. **Engine is C++20.** Pythonic naming optimizes readers who think in pandas;
   production consumers of `include/atomiverse/` are C++ / HPC.

## Why not rename to `get_column` on C++

- Cosmetics over a stable pillar name already used tree-wide.
- Blurs scalar `get` vs column views.
- Real Python users never see C++ headers — they see the wheel. Renaming C++
  does not make the library more pythonic for them.

## Python (and other foreign surfaces)

Allowed and preferred when wrapping Frame:

```text
frame.get_column("x")   # or frame["x"] / frame.x
# implementation: buffer protocol / nanobind over the same SoA storage
```

Binding names are **not** a license to alias or replace `span` in public C++
headers, tutorials that teach C++, or `regressions/` (C++ out-of-library
consumers).

## Related

- `.claude/notes/frame-only-token.md` — Frame is the only public I/O token
- `.claude/notes/frame-region-fields.md` — membership lives on Frame columns
- Architecture: Frame = `blocks` + `meta` + `box`
