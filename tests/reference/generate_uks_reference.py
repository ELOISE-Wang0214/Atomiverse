#!/usr/bin/env python3
"""Generate UKS reference values using PySCF.

Run: python3 generate_uks_reference.py

Requires: pip install pyscf
"""
from pyscf import gto, dft

systems = [
    ('H',  'H 0 0 0',               1,  {'spin': 1}),
    ('He', 'He 0 0 0',              2,  {}),
    ('H2', 'H 0 0 0; H 0 0 0.74',  2,  {}),
    ('Li', 'Li 0 0 0',              3,  {'spin': 1}),
]

functionals = ['lda,vwn5', 'pbe,pbe']

print("UKS Reference Values (STO-3G)")
print("=" * 60)

for name, atom, nelec, kwargs in systems:
    mol = gto.M(atom=atom, basis='sto-3g', **kwargs)
    print(f"\n{name} (nelec={nelec}, {kwargs})")
    print("-" * 40)

    for xc in functionals:
        mf = dft.UKS(mol)
        mf.xc = xc
        mf.grids.level = 4  # standard grid
        mf.verbose = 0
        e = mf.kernel()
        print(f"  {xc:12s}: E = {e:.10f} Ha")

    # Also compute UHF for comparison
    from pyscf import scf
    mf_uhf = scf.UHF(mol)
    mf_uhf.verbose = 0
    e_uhf = mf_uhf.kernel()
    print(f"  {'UHF':12s}: E = {e_uhf:.10f} Ha")

print("\n\nH2 Dissociation Curve (UKS-LSDA vs UHF)")
print("=" * 60)
for r in [0.5, 0.6, 0.74, 1.0, 1.5, 2.0, 2.5, 3.0]:
    mol = gto.M(atom=f'H 0 0 0; H 0 0 {r}', basis='sto-3g')

    mf_uhf = scf.UHF(mol)
    mf_uhf.verbose = 0
    e_uhf = mf_uhf.kernel()

    mf_uks = dft.UKS(mol)
    mf_uks.xc = 'lda,vwn5'
    mf_uks.verbose = 0
    e_uks = mf_uks.kernel()

    print(f"  R={r:5.2f} A: UHF={e_uhf:12.8f}  UKS={e_uks:12.8f}  diff={e_uks-e_uhf:12.8f}")
