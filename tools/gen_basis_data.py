#!/usr/bin/env python
"""Emit an Atomiverse basis_data/*.h static table from PySCF's basis library.

Usage (run manually inside the oracle venv — NEVER wired into CMake/ctest):

    python tools/gen_basis_data.py <pyscf-basis-name> <header-stem> <El> [El ...]
    python tools/gen_basis_data.py --check <header-path> <pyscf-basis-name> \
        <header-stem> <El> [El ...]

Examples:

    python tools/gen_basis_data.py sto-3g sto3g H He Li Be B C N O > \
        include/atomiverse/cpu/scf/basis_data/sto3g.h
    python tools/gen_basis_data.py --check \
        include/atomiverse/cpu/scf/basis_data/sto3g.h sto-3g sto3g H He Li Be B C N O

The emitted numbers come from `pyscf.gto.M(...)._basis` — the raw,
UNNORMALIZED contraction convention (physical GTO normalization is applied
by the C++ ShellList constructor). Output is byte-deterministic: identical
inputs produce identical bytes, and `--check` exits non-zero with a unified
diff when the target header does not match the regenerated text.

Unsupported elements (not present in the requested list, or absent from the
PySCF basis library for the requested basis) raise a hard error — no
fallback, per the repository's no-fallback rule.
"""

import sys

L_LETTER = "spdfgh"

# Periodic table symbols indexed by Z (enough for the basis_data tables).
SYMBOLS = [
    "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn",
]
Z_OF = {sym: z for z, sym in enumerate(SYMBOLS)}

WIDTH = 80


def fetch_basis(basis_name, symbols):
    """Return {symbol: [[l, [e, c], ...], ...]} from PySCF's basis library."""
    from pyscf import gto

    out = {}
    for sym in symbols:
        raw = gto.basis.load(basis_name, sym)
        if not raw:
            raise SystemExit(
                f"error: element {sym} not defined in PySCF basis "
                f"'{basis_name}' — no fallback"
            )
        out[sym] = raw
    return out


def num(x):
    """Shortest exact decimal for a float (Python repr convention)."""
    return repr(float(x))


def wrap_list(values, first_indent, cont_indent):
    """Render {a, b, c} wrapped at WIDTH columns.

    `first_indent` is the column where the opening brace sits; continuation
    lines are indented with `cont_indent` spaces (one past the brace).
    """
    parts = [num(v) for v in values]
    lines = []
    cur = "{"
    col = first_indent + 1
    for i, p in enumerate(parts):
        piece = p + (", " if i < len(parts) - 1 else "")
        if col + len(piece.rstrip()) > WIDTH and cur != "{":
            lines.append(cur.rstrip())
            cur = " " * cont_indent
            col = cont_indent
        cur += piece
        col += len(piece)
    cur += "}"
    lines.append(cur)
    return lines


def render_shell(l, exps, coeffs):
    """Correct multi-line renderer (used instead of shell_lines tail hack)."""
    letter = L_LETTER[l]
    head = f"                ContractedShell{{/*l=*/{l} /* {letter} */,"
    one = (
        head
        + " {"
        + ", ".join(num(e) for e in exps)
        + "}, {"
        + ", ".join(num(c) for c in coeffs)
        + "}},"
    )
    if len(one) <= WIDTH:
        return [one]
    ind = 32
    out = [head]
    e_lines = wrap_list(exps, ind, ind + 1)
    e_lines[0] = " " * ind + e_lines[0]
    e_lines[-1] += ","
    out.extend(e_lines)
    c_lines = wrap_list(coeffs, ind, ind + 1)
    c_lines[0] = " " * ind + c_lines[0]
    c_lines[-1] += "},"
    out.extend(c_lines)
    return out


def render(basis_name, stem, symbols):
    zs = []
    for sym in symbols:
        if sym not in Z_OF:
            raise SystemExit(f"error: unknown element symbol '{sym}'")
        zs.append(Z_OF[sym])
    order = sorted(zip(zs, symbols))
    max_z = order[-1][0]
    by_z = dict(order)

    data = fetch_basis(basis_name, [s for _, s in order])

    syms_comma = ", ".join(s for _, s in order)
    syms_space = " ".join(s for _, s in order)

    lines = []
    a = lines.append
    a(f"/// @file cpu/scf/basis_data/{stem}.h")
    a(f"/// @brief {basis_name} basis-set table (auto-generated; do not edit).")
    a("///")
    a("/// Emitted by `tools/gen_basis_data.py` from PySCF's `mol._basis`")
    a("/// dictionary, so numerical values match the PySCF ground-truth path")
    a("/// bit-for-bit. Contraction coefficients use the unnormalized")
    a("/// convention (physical GTO normalization is applied by")
    a("/// `ShellList` constructor).")
    a("///")
    a(f"/// Supported elements: {syms_comma}.")
    a("///")
    a("/// Regenerate:")
    a(f"///   python tools/gen_basis_data.py {basis_name} {stem} {syms_space}")
    a("#pragma once")
    a("")
    a("#include <vector>")
    a("")
    a('#include "atomiverse/cpu/scf/basis_data/basis_table.h"')
    a("")
    a("namespace ATV {")
    a("namespace cpu {")
    a("namespace scf {")
    a("namespace basis_data {")
    a("")
    a(f"/// @brief Get the {basis_name} table, indexed by atomic number Z.")
    a("/// @return Reference to a static table with AtomBasis entries for")
    a(f"///         Z=1..{max_z}; position 0 and any unsupported Z are empty.")
    a(f"inline const std::vector<AtomBasis>& {stem}() {{")
    a("    static const std::vector<AtomBasis> data = []() {")
    a("        return std::vector<AtomBasis>{")
    a("            AtomBasis{}, // Z=0 unused")
    for z in range(1, max_z + 1):
        if z not in by_z:
            a(f"            AtomBasis{{}}, // Z={z} unsupported")
            continue
        sym = by_z[z]
        shells = data[sym]
        a(f"            // Z={z} ({sym}): {len(shells)} shells")
        a("            AtomBasis{{")
        for shell in shells:
            l = shell[0]
            prims = shell[1:]
            exps = [p[0] for p in prims]
            coeffs = [p[1] for p in prims]
            lines.extend(render_shell(l, exps, coeffs))
        a("            }},")
    a("        };")
    a("    }();")
    a("    return data;")
    a("}")
    a("")
    a("} // namespace basis_data")
    a("} // namespace scf")
    a("} // namespace cpu")
    a("} // namespace ATV")
    return "\n".join(lines) + "\n"


def main(argv):
    if len(argv) >= 2 and argv[1] == "--check":
        if len(argv) < 6:
            raise SystemExit(__doc__)
        path, basis_name, stem = argv[2], argv[3], argv[4]
        symbols = argv[5:]
        text = render(basis_name, stem, symbols)
        with open(path, encoding="utf-8") as fh:
            on_disk = fh.read()
        if text == on_disk:
            print(f"OK: {path} is byte-identical to regenerated output")
            return 0
        import difflib

        diff = difflib.unified_diff(
            on_disk.splitlines(keepends=True),
            text.splitlines(keepends=True),
            fromfile=path,
            tofile="regenerated",
        )
        sys.stdout.writelines(diff)
        print(f"MISMATCH: {path} differs from regenerated output")
        return 1

    if len(argv) < 4:
        raise SystemExit(__doc__)
    basis_name, stem = argv[1], argv[2]
    symbols = argv[3:]
    sys.stdout.write(render(basis_name, stem, symbols))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
