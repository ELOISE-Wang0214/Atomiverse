# Atomiverse

Atomiverse is a C++/CUDA molecular simulation codebase with a shared Driver API
for molecular dynamics and self-consistent electronic-structure workflows.

The repository documentation is now split into three layers:

- user guide
- developer manual
- API reference

## Documentation

Start from these source pages:

- `docs/user-guide/index.md`
- `docs/developer/index.md`
- `docs/api/index.md`

Generate the docs with:

```bash
cmake -B build/docs -DATV_DOCS_ONLY=ON
cmake --build build/docs --target docs
cmake --build build/docs --target docs-open
```

The generated HTML lives in `docs/html/`.

All examples and build recipes in this repository assume **C++17 or later**.

## Deep-Dive Tutorial Order

If you are new to computational chemistry, read the tutorials in this order:

1. [`docs/tutorial/01-scf.md`](docs/tutorial/01-scf.md)
2. [`docs/tutorial/02-rhf.md`](docs/tutorial/02-rhf.md)
3. [`docs/tutorial/03-uhf.md`](docs/tutorial/03-uhf.md)
4. [`docs/tutorial/04-uks.md`](docs/tutorial/04-uks.md)

The first tutorial explains the SCF iteration pattern with a toy matrix model.
The later tutorials show how the same feedback loop becomes RHF, UHF, and UKS.

## Example Programs

These examples are the main entry points for the tutorial material:

| File | Purpose |
|------|---------|
| `examples/toy-scf.cpp` | Explicit CPU toy SCF example |
| `examples/scf.cu` | Driver-based GPU toy SCF example |
| `examples/h2-rhf.cpp` | CPU RHF example using Atomiverse integral and solver APIs |
| `examples/h2_rhf.cu` | GPU RHF example with CPU reference comparison |

## Notes on Scope

The tutorials are intentionally written to match the code that exists in this
repository. In particular:

- `h2-rhf.cpp` is not a handwritten standalone chemistry solver; it already uses
  Atomiverse library APIs
- the GPU RHF path currently reuses CPU-built integrals and focuses the GPU work
  on the SCF loop itself
- the SCF driver structure is shared between the toy SCF and RHF examples
