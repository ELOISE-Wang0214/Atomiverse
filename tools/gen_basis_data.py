#!/usr/bin/env python3
"""Emit C++ headers under include/atomiverse/cpu/scf/basis_data/.

Each header declares one inline function returning the
`std::vector<AtomBasis>` table for a named basis set. The table is
indexed by atomic number Z (0 unused).

Basis data is sourced from PySCF's `mol._basis` dict so the numbers
match the SCF reference bit-for-bit, including the "unnormalized"
contraction-coefficient convention that downstream
`build_shell_list` (CPU) and `CINTgto_norm` (libcint bridge) expect.

Usage:
    python tools/gen_basis_data.py 6-31g  pople_631g  H C N O F
    python tools/gen_basis_data.py 6-31g* pople_631gs H C N O F
    python tools/gen_basis_data.py def2-svp  def2_svp  H B C N O F

The second positional argument is the C++ identifier stem — it becomes
the header filename and the enclosing function name.
"""

from __future__ import annotations

import argparse
import sys
import textwrap
from pathlib import Path

try:
    from pyscf import gto
except ImportError as exc:  # pragma: no cover
    sys.stderr.write("pyscf is required to regenerate basis headers\n")
    raise

REPO_ROOT = Path(__file__).resolve().parent.parent
HEADER_DIR = REPO_ROOT / "include" / "atomiverse" / "cpu" / "scf" / "basis_data"

SYMBOL_TO_Z = {
    "H": 1,  "He": 2, "Li": 3, "Be": 4, "B": 5,  "C": 6,
    "N": 7,  "O": 8,  "F": 9,  "Ne": 10, "Na": 11, "Mg": 12,
    "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18,
}

L_NAME = {0: "s", 1: "p", 2: "d", 3: "f"}


def fetch_atom_basis(basis: str, symbol: str):
    """Return the PySCF `_basis[symbol]` entry for a neutral atom."""
    # PySCF needs a valid electronic structure. Use a closed-shell trick:
    # for odd-electron atoms, attach a bonding partner at a dummy distance.
    atom_spec = f"{symbol} 0 0 0"
    z = SYMBOL_TO_Z[symbol]
    spin = z % 2          # 0 for even, 1 for odd
    try:
        mol = gto.M(atom=atom_spec, basis=basis, spin=spin, verbose=0)
    except RuntimeError:
        # Fall back: anchor a neutral H at a large distance.
        mol = gto.M(atom=f"{symbol} 0 0 0; H 100 0 0",
                    basis=basis, spin=(spin + 1) % 2, verbose=0)
    return mol._basis[symbol]


def emit_shell(l: int, prims: list[list[float]], indent: str) -> str:
    """Format one ContractedShell initializer.

    PySCF stores a shell as `[l, [alpha_0, c_0], [alpha_1, c_1], ...]`,
    and for segmented-contraction bases every primitive has a single
    coefficient column. We reject the (rare) "general contraction"
    multi-column case — those bases would need a different schema.
    """
    alphas = []
    coeffs = []
    for prim in prims:
        alpha = prim[0]
        c = prim[1:]
        if len(c) != 1:
            raise ValueError(
                f"general-contraction shell (l={l}, ncols={len(c)}) not supported"
            )
        alphas.append(alpha)
        coeffs.append(c[0])

    alpha_str = ", ".join(f"{a!r}" for a in alphas)
    coeff_str = ", ".join(f"{c!r}" for c in coeffs)
    return (
        f"{indent}ContractedShell{{ /*l=*/{l} /* {L_NAME.get(l, '?')} */, "
        f"{{{alpha_str}}}, {{{coeff_str}}} }},"
    )


def emit_header(basis_name: str, stem: str, symbols: list[str]) -> str:
    max_z = max(SYMBOL_TO_Z[s] for s in symbols)
    table_slots = max_z + 1   # position 0 unused

    # Build Z -> (symbol, shells) map so the emitted table has one entry per Z
    # up to max_z, with empty AtomBasis{} for unsupported positions.
    by_z: dict[int, tuple[str, list]] = {}
    for sym in symbols:
        by_z[SYMBOL_TO_Z[sym]] = (sym, fetch_atom_basis(basis_name, sym))

    lines = []
    lines.append("    return std::vector<AtomBasis>{")
    for z in range(table_slots):
        if z == 0:
            lines.append("        AtomBasis{},  // Z=0 unused")
            continue
        if z not in by_z:
            lines.append(f"        AtomBasis{{}},  // Z={z} unsupported")
            continue
        sym, shells = by_z[z]
        lines.append(f"        // Z={z} ({sym}): {len(shells)} shells")
        lines.append("        AtomBasis{{")
        for sh in shells:
            l = sh[0]
            prims = sh[1:]
            lines.append(emit_shell(l, prims, indent=" " * 12))
        lines.append("        }},")
    lines.append("    };")

    header = f"""\
/// @file cpu/scf/basis_data/{stem}.h
/// @brief {basis_name} basis-set table (auto-generated; do not edit).
///
/// Emitted by `tools/gen_basis_data.py` from PySCF's `mol._basis`
/// dictionary, so numerical values match the PySCF ground-truth path
/// bit-for-bit. Contraction coefficients use the unnormalized
/// convention (physical GTO normalization is applied by
/// `build_shell_list`).
///
/// Supported elements: {', '.join(symbols)}.
///
/// Regenerate:
///   python tools/gen_basis_data.py {basis_name} {stem} {' '.join(symbols)}
#pragma once

#include <vector>

#include "atomiverse/cpu/scf/basis_data/basis_table.h"

namespace ATV {{
namespace cpu {{
namespace scf {{
namespace basis_data {{

/// @brief Get the {basis_name} table, indexed by atomic number Z.
/// @return Reference to a static table with AtomBasis entries for
///         Z=1..{max_z}; position 0 and any unsupported Z are empty.
inline const std::vector<AtomBasis>& {stem}() {{
    static const std::vector<AtomBasis> data = []() {{
{chr(10).join(lines)}
    }}();
    return data;
}}

}} // namespace basis_data
}} // namespace scf
}} // namespace cpu
}} // namespace ATV
"""
    return header


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("basis_name",
                        help="PySCF basis name, e.g. '6-31g' or 'def2-svp'")
    parser.add_argument("stem",
                        help="C++ identifier stem, e.g. 'pople_631g'")
    parser.add_argument("symbols", nargs="+",
                        help="Element symbols to include (e.g. H C N O F)")
    args = parser.parse_args()

    unknown = [s for s in args.symbols if s not in SYMBOL_TO_Z]
    if unknown:
        parser.error(f"unknown element symbol(s): {unknown}")

    HEADER_DIR.mkdir(parents=True, exist_ok=True)
    out_path = HEADER_DIR / f"{args.stem}.h"
    header = emit_header(args.basis_name, args.stem, args.symbols)
    out_path.write_text(header)
    sys.stdout.write(f"wrote {out_path.relative_to(REPO_ROOT)}\n")


if __name__ == "__main__":
    main()
