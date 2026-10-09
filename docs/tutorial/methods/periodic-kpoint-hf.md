# Periodic k-point Hartree–Fock

Restricted and unrestricted Hartree–Fock with a Monkhorst–Pack mesh on
an orthorhombic cell. The Coulomb and exchange matrices are the PySCF
FFT density-fitting kernels (discrete Fourier transform of the AO
products, Coulomb kernel `4π/G²`, Ewald probe charge on the exchange
`G = 0` term). The nuclear repulsion is the Ewald sum.

```cpp
#include "atomiverse/cpu/scf/krhf.h"
#include "atomiverse/cpu/scf/kuhf.h"
#include "atomiverse/cuda/scf/krhf.cuh"

cpu::Frame frame = /* orthorhombic cell, positions in Å */;
cpu::scf::KRHF krhf{cpu::scf::BasisSet::STO3G, /*k_mesh=*/{2, 1, 1}};
krhf.run(frame);                 // kEnergy and kEnergyElectronic [eV]
double e_ha = krhf.energy();     // Hartree

// Mesh totals, not a per-k copy. Triplet H₂ with cell.spin = 2 on a
// 2-point mesh is Nα = 3, Nβ = 1 (PySCF KUHF.nelec).
cpu::scf::KUHF kuhf{cpu::scf::BasisSet::STO3G, {2, 1, 1},
                    /*n_alpha=*/3, /*n_beta=*/1};
kuhf.run(frame);

// Restricted SCF on the device. The step (Coulomb, exchange, Fock,
// diagonalize, density) stays on device; l > 0 throws (`krhf_ao_lgt0`).
// There is no device KUHF.
cuda::scf::KRHF gpu{cpu::scf::BasisSet::STO3G, {2, 1, 1}};
gpu.run(frame);
```

`k_mesh = {1,1,1}` is the single Γ point of this kernel (not the
truncated real-space lattice behind `RHF::periodic`). Shells with
angular momentum above s throw (`krhf_ao_lgt0`). `forces` throws
(`krhf_nuclear_forces`) — there is no analytic nuclear gradient on
this path. `provides_force()` is false.

H₂ / STO-3G in a 4 Bohr cube (shorter than the STO-3G pair reach) has
PySCF total energies −1.4471748215213558 Ha for `KRHF` with
`k_mesh = {2,1,1}`, and −0.7658470471886634 Ha for triplet `KUHF` on that
mesh (`cell.spin = 2`, Nα = 3, Nβ = 1).
