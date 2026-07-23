#!/usr/bin/env python3
"""Fair wall-time: offline PySCF vs Atomiverse binary (never a ctest dep).

Gate for the method loop (docs/developer/pyscf-method-loop.md):
  within_10x  ⇔  T_atv_median / T_pyscf_median ≤ max_ratio  (default 10)

Examples:
  # PySCF-only timing (kernel wall after import)
  python tools/pyscf_oracle/fair_bench.py --method rhf --geom water --basis 6-31g*

  # Fair A/B
  python tools/pyscf_oracle/fair_bench.py --method ump2 --geom o2 --basis sto-3g \\
      --spin 2 --repeat 3 -- -- ./build/regressions/reg_ump2-o2-sto3g

Protocol:
  - Same geometry / basis / method class on both sides when possible.
  - Report median of --repeat runs (default 3).
  - Tiny molecules may be startup-dominated; still record ratio for the gate.
  - Extend --method as new families land (mp2, tda, ump2, …).
"""
from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import sys
import time


def _median(xs: list[float]) -> float:
    return float(statistics.median(xs)) if xs else float("nan")


def _mol_code(geom: str, basis: str, charge: int, spin: int, cart: bool | None) -> str:
    presets = {
        "water": """
O  0.000000  0.000000  0.117300
H  0.000000  0.757200 -0.469200
H  0.000000 -0.757200 -0.469200
""",
        "h2": """
H  0.000000  0.000000  0.000000
H  0.000000  0.000000  0.740000
""",
        "o2": """
O  0.000000  0.000000  0.000000
O  0.000000  0.000000  1.207500
""",
        "h2co": """
C  0.000000  0.000000  0.000000
O  0.000000  0.000000  1.210000
H  0.000000  0.943000 -0.587000
H  0.000000 -0.943000 -0.587000
""",
        "h3": """
H  0.000000  0.000000  0.000000
H  0.000000  0.000000  0.900000
H  0.000000  0.850000  0.450000
""",
    }
    atom = presets.get(geom)
    if atom is None:
        raise SystemExit(f"unknown geom {geom!r}; choose {sorted(presets)}")
    if cart is None:
        cart = not basis.lower().startswith("def2")
    return f"""
from pyscf import gto
mol = gto.M(atom='''{atom}''', basis={basis!r}, charge={charge}, spin={spin},
            unit='Angstrom', cart={cart}, verbose=0)
"""


def pyscf_kernel_seconds(method: str, geom: str, basis: str, charge: int, spin: int,
                         cart: bool | None, xc: str, nstates: int) -> float:
    setup = _mol_code(geom, basis, charge, spin, cart)
    body = {
        "rhf": "from pyscf import scf\nmf=scf.RHF(mol); mf.conv_tol=1e-10; mf.kernel()",
        "uhf": "from pyscf import scf\nmf=scf.UHF(mol); mf.conv_tol=1e-10; mf.kernel()",
        "rohf": "from pyscf import scf\nmf=scf.ROHF(mol); mf.conv_tol=1e-10; mf.kernel()",
        "rks": f"from pyscf import dft\nmf=dft.RKS(mol); mf.xc={xc!r}; mf.conv_tol=1e-10; mf.kernel()",
        "uks": f"from pyscf import dft\nmf=dft.UKS(mol); mf.xc={xc!r}; mf.conv_tol=1e-10; mf.kernel()",
        "mp2": "from pyscf import scf, mp\nmf=scf.RHF(mol); mf.conv_tol=1e-10; mf.kernel(); mp.MP2(mf).run()",
        "ump2": "from pyscf import scf, mp\nmf=scf.UHF(mol); mf.conv_tol=1e-10; mf.kernel(); mp.UMP2(mf).run()",
        "tda": (
            "from pyscf import scf, tdscf\nmf=scf.RHF(mol); mf.conv_tol=1e-10; mf.kernel(); "
            f"td=tdscf.TDA(mf); td.nstates={nstates}; td.kernel()"
        ),
        "tdhf": (
            "from pyscf import scf, tdscf\nmf=scf.RHF(mol); mf.conv_tol=1e-10; mf.kernel(); "
            f"td=tdscf.TDHF(mf); td.nstates={nstates}; td.kernel()"
        ),
        "rks_tda": (
            f"from pyscf import dft, tdscf\nmf=dft.RKS(mol); mf.xc={xc!r}; mf.conv_tol=1e-10; mf.kernel(); "
            f"td=tdscf.TDA(mf); td.nstates={nstates}; td.kernel()"
        ),
        "roks": (
            f"from pyscf import dft\nmf=dft.ROKS(mol); mf.xc={xc!r}; mf.conv_tol=1e-10; mf.kernel()"
        ),
        "roks_grad": (
            f"from pyscf import dft\nmf=dft.ROKS(mol); mf.xc={xc!r}; mf.conv_tol=1e-10; mf.kernel(); "
            "mf.nuc_grad_method().kernel()"
        ),
    }.get(method)
    if body is None:
        raise SystemExit(f"unsupported --method {method!r}")
    code = setup + "\nimport time\nt0=time.perf_counter()\n" + body + "\nprint(time.perf_counter()-t0)\n"
    r = subprocess.run([sys.executable, "-c", code], check=True, capture_output=True, text=True)
    return float(r.stdout.strip().splitlines()[-1])


def bench_cmd(cmd: list[str]) -> float:
    t0 = time.perf_counter()
    # Capture child I/O so PASS/FAIL banners do not pollute JSON stdout.
    subprocess.run(cmd, check=True, capture_output=True)
    return time.perf_counter() - t0


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--method", default="rhf",
                   choices=["rhf", "uhf", "rohf", "rks", "uks", "roks", "roks_grad", "mp2", "ump2", "tda", "tdhf", "rks_tda"])
    p.add_argument("--geom", default="water")
    p.add_argument("--basis", default="sto-3g")
    p.add_argument("--charge", type=int, default=0)
    p.add_argument("--spin", type=int, default=0)
    p.add_argument("--xc", default="b3lyp")
    p.add_argument("--nstates", type=int, default=3)
    p.add_argument("--cart", action=argparse.BooleanOptionalAction, default=None)
    p.add_argument("--repeat", type=int, default=3)
    p.add_argument("--max-ratio", type=float, default=10.0,
                   help="within_10x if T_atv/T_pyscf <= this (default 10)")
    p.add_argument("atv_cmd", nargs=argparse.REMAINDER,
                   help="optional: -- ./path/to/reg_binary [args…]")
    args = p.parse_args()

    out: dict = {
        "ok": True,
        "method": args.method,
        "geom": args.geom,
        "basis": args.basis,
        "max_ratio": args.max_ratio,
        "repeat": args.repeat,
    }

    try:
        pt = [pyscf_kernel_seconds(args.method, args.geom, args.basis, args.charge, args.spin,
                                   args.cart, args.xc, args.nstates)
              for _ in range(args.repeat)]
        out["pyscf_wall_s"] = pt
        out["pyscf_median_s"] = _median(pt)
    except Exception as e:
        out["ok"] = False
        out["pyscf_error"] = str(e)
        print(json.dumps(out, indent=2, sort_keys=True))
        return 2

    cmd = list(args.atv_cmd or [])
    if cmd and cmd[0] == "--":
        cmd = cmd[1:]
    if cmd:
        try:
            at = [bench_cmd(cmd) for _ in range(args.repeat)]
            out["atv_wall_s"] = at
            out["atv_median_s"] = _median(at)
            if out["pyscf_median_s"] > 0:
                ratio = out["atv_median_s"] / out["pyscf_median_s"]
                out["ratio_atv_over_pyscf"] = ratio
                out["within_10x"] = bool(ratio <= args.max_ratio)
        except Exception as e:
            out["ok"] = False
            out["atv_error"] = str(e)

    print(json.dumps(out, indent=2, sort_keys=True))
    if not out.get("ok"):
        return 1
    if cmd and out.get("within_10x") is False:
        return 3  # correct path may still exist; gate failed on speed
    return 0


if __name__ == "__main__":
    sys.exit(main())
