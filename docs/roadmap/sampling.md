# Sampling and Free Energy

Legend and caveats: [Roadmap index](index.md).

Monte Carlo exists as a first-class Driver category (`cpu::mc::MC` on
`SamplingDriver`) with the inline Metropolis accept/reject test and a
transaction/rollback mechanism for trial moves — but it ships exactly one
ensemble (NVT) and one mover (single-particle translation), and it is CPU-only.
It reuses the MD force pipeline for trial energies, with the neighbor list in
always-rebuild mode. Validated against NIST LJ reference data
(`benchmarks/mc/lj/validation/`).

Enhanced sampling and free-energy machinery do not exist at all. Neither does
any biasing/collective-variable infrastructure to build them on — there is no
CV layer, no bias component, no PLUMED interface.

## 1. Monte Carlo Methods

| Method                                      | CPU | CUDA |
| ------------------------------------------- | :-: | :--: |
| Metropolis Monte Carlo                      |  ✓  |   ○  |
| Canonical Monte Carlo, NVT                  |  ✓  |   ○  |
| Isothermal–isobaric Monte Carlo, NPT        |  ○  |   ○  |
| Grand-canonical Monte Carlo                 |  ○  |   ○  |
| Gibbs-ensemble Monte Carlo                  |  ○  |   ○  |
| Semi-grand-canonical Monte Carlo            |  ○  |   ○  |
| Configurational-bias Monte Carlo            |  ○  |   ○  |
| Continuous-fractional-component Monte Carlo |  ○  |   ○  |
| Hybrid Monte Carlo                          |  ○  |   ○  |
| Hamiltonian Monte Carlo                     |  ○  |   ○  |
| Replica-exchange Monte Carlo                |  ○  |   ○  |
| Simulated tempering                         |  ○  |   ○  |
| Wang–Landau sampling                        |  ○  |   ○  |
| Multicanonical sampling                     |  ○  |   ○  |
| Nested sampling                             |  ○  |   ○  |
| Kinetic Monte Carlo                         |  ○  |   ○  |
| Rejection-free kinetic Monte Carlo          |  ○  |   ○  |
| Reaction-ensemble Monte Carlo               |  ○  |   ○  |

The missing ensembles are mover-shaped, not architecture-shaped: NPT-MC needs a
volume mover, GCMC an insertion/deletion mover. The `Transaction` proposal
channel (`ctx.mc.proposal`) was designed for exactly that, but only
`TranslateMove` is written.

## 2. Enhanced-Sampling Methods

| Method                                  | CPU | CUDA |
| --------------------------------------- | :-: | :--: |
| Umbrella sampling                       |  ○  |   ○  |
| Adaptive umbrella sampling              |  ○  |   ○  |
| Multidimensional umbrella sampling      |  ○  |   ○  |
| Well-tempered metadynamics              |  ○  |   ○  |
| Multiple-walker metadynamics            |  ○  |   ○  |
| Bias-exchange metadynamics              |  ○  |   ○  |
| Parallel-bias metadynamics              |  ○  |   ○  |
| Funnel metadynamics                     |  ○  |   ○  |
| Adaptive biasing force                  |  ○  |   ○  |
| Extended-system adaptive biasing force  |  ○  |   ○  |
| Adaptive biasing potential              |  ○  |   ○  |
| Accelerated weight histogram            |  ○  |   ○  |
| Variationally enhanced sampling         |  ○  |   ○  |
| OPES                                    |  ○  |   ○  |
| OPES Metad                              |  ○  |   ○  |
| OPES Expanded                           |  ○  |   ○  |
| Gaussian accelerated molecular dynamics |  ○  |   ○  |
| Temperature-accelerated MD              |  ○  |   ○  |
| Hyperdynamics                           |  ○  |   ○  |
| Integrated tempering sampling           |  ○  |   ○  |
| Solute tempering                        |  ○  |   ○  |
| REST2                                   |  ○  |   ○  |
| Replica-exchange molecular dynamics     |  ○  |   ○  |
| Hamiltonian replica exchange            |  ○  |   ○  |
| Multidimensional replica exchange       |  ○  |   ○  |
| Expanded-ensemble sampling              |  ○  |   ○  |
| Weighted-ensemble sampling              |  ○  |   ○  |
| Transition-path sampling                |  ○  |   ○  |
| Transition-interface sampling           |  ○  |   ○  |
| Forward-flux sampling                   |  ○  |   ○  |
| Milestoning                             |  ○  |   ○  |

## 3. Free-Energy Methods

| Method                                  | CPU | CUDA |
| --------------------------------------- | :-: | :--: |
| Free-energy perturbation                |  ○  |   ○  |
| Thermodynamic integration               |  ○  |   ○  |
| Multistate Bennett acceptance ratio     |  ○  |   ○  |
| Umbrella integration                    |  ○  |   ○  |
| Potential-of-mean-force calculation     |  ○  |   ○  |
| Adaptive integration of mean force      |  ○  |   ○  |
| Jarzynski equality                      |  ○  |   ○  |
| Crooks fluctuation theorem              |  ○  |   ○  |
| Nonequilibrium switching                |  ○  |   ○  |
| Bidirectional nonequilibrium estimation |  ○  |   ○  |
| Widom insertion                         |  ○  |   ○  |
| Absolute solvation free energy          |  ○  |   ○  |
| Relative binding free energy            |  ○  |   ○  |
| Absolute binding free energy            |  ○  |   ○  |
| Double-decoupling method                |  ○  |   ○  |
| Alchemical lambda dynamics              |  ○  |   ○  |
| Constant-pH free-energy calculation     |  ○  |   ○  |

Alchemical free energy additionally needs soft-core potentials, which are also
absent (see [Force Fields](force-fields.md#2-classical-nonbonded-potentials)).
