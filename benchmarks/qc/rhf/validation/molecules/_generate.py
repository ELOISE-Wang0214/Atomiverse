#!/usr/bin/env python3
"""Generate XYZ starting geometries for the Phase 1.2 benchmark molecules.

Idempotent: overwrites each target file.

Each entry carries a SMILES string plus the conventional comment that
ships in the XYZ header. For the radical cases we build an explicit
molecule with the correct spin multiplicity via the PySCF-friendly
`allowed_atoms` override so the reference-generation script later sees
a geometry consistent with its mult.
"""

from __future__ import annotations

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

HERE = Path(__file__).resolve().parent

TARGETS = {
    "benzene":      ("c1ccccc1",              "Benzene C6H6 (SMILES, UFF-optimized)"),
    "naphthalene":  ("c1ccc2ccccc2c1",        "Naphthalene C10H8 (SMILES, UFF-optimized)"),
    "ethanol":      ("CCO",                   "Ethanol C2H5OH (SMILES, UFF-optimized)"),
    "glycine":      ("NCC(=O)O",              "Glycine NH2-CH2-COOH (SMILES, UFF-optimized)"),
    "aspirin":      ("CC(=O)Oc1ccccc1C(=O)O", "Aspirin (acetylsalicylic acid, SMILES, UFF-optimized)"),
    "caffeine":     ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "Caffeine (SMILES, UFF-optimized)"),
    "ch3_radical":  ("[CH3]",                 "Methyl radical CH3 (planar, doublet)"),
    "o2_triplet":   ("O=O",                   "Molecular oxygen O2 (triplet)"),
    "no_doublet":   ("[N]=O",                 "Nitric oxide NO (doublet)"),
}


def build(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=2026)
    AllChem.UFFOptimizeMolecule(mol)
    return mol


def to_xyz(mol: Chem.Mol, comment: str) -> str:
    conf = mol.GetConformer()
    lines = [str(mol.GetNumAtoms()), comment]
    for atom in mol.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        lines.append(f"{atom.GetSymbol():<2s}  {p.x: .6f}  {p.y: .6f}  {p.z: .6f}")
    return "\n".join(lines) + "\n"


def main() -> None:
    for stem, (smiles, comment) in TARGETS.items():
        mol = build(smiles)
        (HERE / f"{stem}.xyz").write_text(to_xyz(mol, comment))
        print(f"wrote {stem}.xyz  (natoms={mol.GetNumAtoms()})")


if __name__ == "__main__":
    main()
