# molrs product gate (derived dependency)

**Status:** binding (2026-07-29). Router: `cmake/AtvProductModules.cmake`.

## Rule

**Products control molrs — molrs does not control products.**

| User / product switch | Effect |
|----------------------|--------|
| `ATV_BUILD_ANALYSIS` | RDF/MSD/VACF/… sources **and** derives `ATV_BUILD_MOLRS=ON` |
| `ATV_BUILD_IO` | extxyz / frame_writer / XYZ recorders **and** derives MOLRS |
| `ATV_BUILD_PYTHON` | forces IO (exchange) → derives MOLRS |
| `ATV_BUILD_GPW` / `HF` / `SEMIEMPIRICAL` | Element / AM1-BCC / ECP → derives MOLRS |

`ATV_BUILD_MOLRS` is **computed** (CACHE FORCE), not a peer product option.
Never write `FATAL_ERROR "ANALYSIS requires MOLRS"` — enable MOLRS from the
product instead.

## Architecture note

molrs is three concerns under one Rust crate today:

1. **Analysis** (RDF/MSD/…) — pure “engine measures, molrs analyzes”
2. **I/O + Frame exchange** (`from_molrs` / `to_molrs`)
3. **`molrs::Element`** for table lookups (GPW/EAM/SW/ECP)

(1)+(2) match “only when analysis/recorder modules are on.”  
(3) still pulls molrs for method modules that use Element; long-term may move
Element into ATV if pure-engine builds must avoid Rust.

## Wrong

```cmake
if(ATV_BUILD_ANALYSIS AND NOT ATV_BUILD_MOLRS)
  message(FATAL_ERROR "ANALYSIS requires MOLRS")  # inverted
endif()
option(ATV_BUILD_MOLRS "…" ON)  # peer product with its own default
```

## Right

```cmake
if(ATV_BUILD_ANALYSIS OR ATV_BUILD_IO OR …)
  set(ATV_BUILD_MOLRS ON CACHE BOOL "derived" FORCE)
endif()
# then if(ATV_BUILD_ANALYSIS) add rdf sources
#      if(ATV_BUILD_MOLRS) corrosion + exchange
```
