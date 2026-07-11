# Getting Started

This page takes you from a fresh checkout to a real, converged calculation. It
covers only the **C++** workflow. By the end you will have built and run an
electronic-structure calculation and know how to read its output.

Atomiverse offers several methods through one common API: classical molecular
dynamics, geometry optimization, Gaussian-basis SCF (Hartree-Fock and DFT for
molecules), plane-wave DFT (for crystals), and ab-initio MD. You do **not** need
to learn all of them at once — build one example, run it, read the matching
[tutorial](../tutorial/index.md), then move on.

Everything assumes **C++17 or later**.

## 1. What you need

- A C++17 compiler and **CMake ≥ 3.18**.
- **BLAS/LAPACK** (OpenBLAS) — the SCF solvers use dense linear algebra.
- A **Rust toolchain** — the molecular-I/O dependency (molrs) is built from
  Rust source automatically by the build.
- *Optional:* **FFTW3** for plane-wave DFT (`ATV_BUILD_PW`); the **CUDA
  toolkit** + a GPU for the CUDA backend (`ATV_BUILD_CUDA`).

The main CMake options you will use:

| Option | Turns on |
|---|---|
| `ATV_BUILD_RHF` | Gaussian SCF: RHF, UHF, UKS (and AIMD / geometry opt with SCF forces) |
| `ATV_BUILD_PW` | Plane-wave DFT (needs FFTW3) |
| `ATV_BUILD_CUDA` | the CUDA backend (needs a GPU) |
| `ATV_BUILD_EXAMPLES` | the example programs |
| `ATV_BUILD_TESTS` | the test suite |

## 2. Set up the build environment

**On a generic Linux machine** with `gcc`, OpenBLAS and CMake already on your
`PATH`, you can skip straight to step 3.

**On an HPC cluster that uses Lmod modules** (this is the project's setup on the
Arrhenius login node), two things commonly trip people up:

1. The toolchain is not loaded by default — load it in a **login shell**:
   ```bash
   module load buildenv-gcc/2026.03-mpich      # GCC 14.3 + OpenBLAS + FFTW
   ```
2. The module sets the *link-time* library path but not the *run-time* one.
   Mirror it before building or running anything, or binaries fail with
   `libopenblas.so.0: cannot open shared object file` /
   `GLIBCXX_3.4.32 not found`:
   ```bash
   export OPENBLAS_ROOT="$OPENBLASDIR"
   export LD_LIBRARY_PATH="$LIBRARY_PATH:$LD_LIBRARY_PATH"
   export RUSTUP_TOOLCHAIN=stable-x86_64-unknown-linux-gnu   # x86_64 login node only
   ```

> The `RUSTUP_TOOLCHAIN` line is only needed on the x86_64 login node, whose
> default Rust toolchain is (incorrectly) the aarch64 one. Omit it elsewhere.

## 3. Configure and build

The smallest useful starting point is the Gaussian-SCF build with examples:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release \
      -DATV_BUILD_CUDA=OFF -DATV_BUILD_RHF=ON -DATV_BUILD_EXAMPLES=ON \
      -DBLA_VENDOR=OpenBLAS
cmake --build build --target 02_rhf -j
```

Add `-DATV_BUILD_PW=ON` if you also want the plane-wave example, and
`-DATV_BUILD_TESTS=ON` to build the test suite.

> **molrs path (worktrees / non-standard checkouts):** molrs is found at
> `../molrs` by default. If your checkout is not beside a `molrs/` directory
> (e.g. you are in a git worktree), pass an absolute path:
> `-DATV_MOLRS_PATH=/abs/path/to/molrs`.

## 4. Run your first calculation

```bash
./build/examples/02_rhf
```

This computes the restricted Hartree-Fock energy of the hydrogen molecule in
the STO-3G basis. The tail of the output is:

```
=== RHF Result ===
Converged:        yes (3 iterations)
E_electronic:     -1.8318636465 Hartree
V_nuclear:         0.7151043391 Hartree
E_total:          -1.1167593074 Hartree

--- Textbook comparison (Szabo & Ostlund) ---
  Reference: -1.1168 Hartree
  This:      -1.1168 Hartree
```

The total energy is the sum of the electronic energy and the nuclear-repulsion
energy, and it matches the published value. That is a full electronic-structure
calculation: atomic coordinates in, total energy out. The matching theory and a
line-by-line walk-through are in
[Chapter 2: Restricted Hartree-Fock](../tutorial/ch02-rhf.md).

## 5. Which example next?

Build any target the same way (`cmake --build build --target <name>`). Pair each
with its tutorial chapter.

| Program | You want to learn | Build flag |
|---|---|---|
| `01_scf` | the SCF feedback loop, no chemistry | — |
| `02_rhf` | your first real chemistry result (RHF) | `ATV_BUILD_RHF` |
| `03_uhf` | open-shell / unrestricted spin | `ATV_BUILD_RHF` |
| `04_uks` | Kohn-Sham DFT on a grid | `ATV_BUILD_RHF` |
| `05_geom_opt` | relaxing a molecule's geometry | `ATV_BUILD_RHF` |
| `06_aimd` | ab-initio molecular dynamics | `ATV_BUILD_RHF` |
| `07_md_lj` | classical MD (LJ argon, NVE) — the canonical CPU lesson | always built |
| `07_md_lj_cuda` | the CUDA twin of `07_md_lj` | `ATV_BUILD_CUDA` |
| `08_pw` | plane-wave DFT for a crystal | `ATV_BUILD_PW` |
| `09_gpw` | GPW DFT (Gaussian + plane waves), H₂O | `ATV_BUILD_GPW` |

## 6. The one data structure to remember

The user-facing description of your system lives in
[`cpu::Frame`](https://github.com/MolCrafts/Atomiverse/blob/dev/include/atomiverse/cpu/frame/frame.h). You fill it with
atom positions, masses, atom types, the simulation box, and optionally topology
(bonds/angles). Positions are in Ångström; energies come back on
`frame.scalars.energy` in eV (the SCF/PW evaluators also report Hartree/Rydberg
on their result objects). CPU workflows operate on the frame directly; CUDA
workflows copy it to a device frame internally and copy results back, so you do
not manage GPU memory by hand for a standard run.

## 7. Recommended reading order

1. This page.
2. [Chapter 0: A Tour of the CPU Backend](../tutorial/ch00-architecture.md) —
   the Driver / Component / Frame vocabulary.
3. [Chapter 2 → 5](../tutorial/ch02-rhf.md): SCF → RHF → UHF → UKS (the quantum
   chemistry arc).
4. [Chapter 8: Plane-Wave DFT](../tutorial/ch08-plane-wave.md) for periodic systems.
5. [Chapter 10: Classical MD](../tutorial/ch10-classical-md.md),
   [Chapter 15: AIMD](../tutorial/ch15-cross-method.md), and
   [Chapter 12: Geometry Minimization](../tutorial/ch12-geometry-optimization.md).

If classical force fields are what brought you here, start at Chapter 10
instead; if electronic structure is, start at Chapter 2.
