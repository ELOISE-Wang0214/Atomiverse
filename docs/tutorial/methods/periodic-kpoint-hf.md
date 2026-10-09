# Periodic k-point Hartree–Fock

Restricted Hartree–Fock with a Monkhorst–Pack mesh on an orthorhombic
cell. The Coulomb and exchange matrices are the PySCF FFT density-fitting
kernels (discrete Fourier transform of the AO products, Coulomb kernel
`4π/G²`, Ewald probe charge on the exchange `G = 0` term). The nuclear
repulsion is the Ewald sum.

```cpp
#include "atomiverse/cpu/scf/krhf.h"
#include "atomiverse/cuda/scf/krhf.cuh"

cpu::Frame frame = /* orthorhombic cell, positions in Å */;
cpu::scf::KRHF krhf{cpu::scf::BasisSet::STO3G, /*k_mesh=*/{1, 1, 2}};
krhf.run(frame);                 // frame.meta kEnergy [eV]
double e_ha = krhf.energy();     // Hartree

// Same arguments on the device. The SCF step (J, K, diagonalize,
// density) stays on device; l > 0 throws.
cuda::scf::KRHF gpu{cpu::scf::BasisSet::STO3G, {1, 1, 2}};
gpu.run(frame);
```

`k_mesh = {1,1,1}` is the single Γ point of this kernel (not the
truncated real-space lattice behind `RHF::periodic`). Shells with
angular momentum above s throw. `forces` throws — there is no analytic
nuclear gradient on this path.

The reference energy for H₂ / STO-3G in a 4 Bohr cube with
`k_mesh = {1,1,2}` is the PySCF `KRHF` total energy
−1.354143879833 Ha.
