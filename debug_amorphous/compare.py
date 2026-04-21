#!/usr/bin/env python3
"""compare.py — diff LAMMPS vs Atomiverse thermo logs for amorphous PE.

Reads:
  log.PE_nc10_cl1000.txt   (LAMMPS, "step temp press" thermo style)
  log.atomiverse.txt       (Atomiverse, our printer's "step T P PE KE Lx Ly Lz")

Output:
  - prints a side-by-side text table (sampled at LAMMPS-print intervals)
  - writes compare.png with three subplots: T(step), P(step), Lx(step)
    overlaid for each engine and each LAMMPS stage.

Atomiverse pressure is sign-flipped from the physical (LAMMPS) convention
internally; the runner prints the LAMMPS-equivalent (see run_atomiverse.cpp).
Pressure unit on both sides: atm. Temperature: K.

Usage: python compare.py [--lammps log.PE_nc10_cl1000.txt] \
                          [--atomiverse log.atomiverse.txt] \
                          [--out compare.png]
"""
from __future__ import annotations

import argparse
import re
from dataclasses import dataclass, field
from pathlib import Path

# ────────────────────────────────────────────────────────────────────


@dataclass
class Stage:
    label: str
    step: list[int] = field(default_factory=list)
    temp: list[float] = field(default_factory=list)
    press: list[float] = field(default_factory=list)
    pot_e: list[float] = field(default_factory=list)
    kin_e: list[float] = field(default_factory=list)
    lx: list[float] = field(default_factory=list)


# ── LAMMPS log parser ──────────────────────────────────────────────


_LAMMPS_HEADER = re.compile(
    r"^\s*Step\s+Temp\s+Press\b", re.IGNORECASE)
_LAMMPS_RUN_LINE = re.compile(
    r"^\s*\d+\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s*$")


def parse_lammps(path: Path) -> list[Stage]:
    stages: list[Stage] = []
    cur: Stage | None = None
    in_thermo = False
    stage_idx = 0
    with path.open() as f:
        for line in f:
            if _LAMMPS_HEADER.search(line):
                stage_idx += 1
                cur = Stage(label=f"LAMMPS Stage {stage_idx}")
                stages.append(cur)
                in_thermo = True
                continue
            if in_thermo and _LAMMPS_RUN_LINE.match(line):
                tok = line.split()
                cur.step.append(int(tok[0]))
                cur.temp.append(float(tok[1]))
                cur.press.append(float(tok[2]))
                continue
            # Closing a thermo block: any non-data line ends the run.
            if in_thermo and line.strip() and not _LAMMPS_RUN_LINE.match(line):
                in_thermo = False
    return stages


# ── Atomiverse log parser ─────────────────────────────────────────


def parse_atomiverse(path: Path) -> list[Stage]:
    stages: list[Stage] = []
    cur: Stage | None = None
    with path.open() as f:
        for line in f:
            ls = line.strip()
            if ls.startswith("# Stage"):
                cur = Stage(label=ls.lstrip("# ").strip())
                stages.append(cur)
                continue
            if ls.startswith("#") or not ls:
                continue
            tok = ls.split()
            if len(tok) < 8 or cur is None:
                continue
            try:
                cur.step.append(int(tok[0]))
                cur.temp.append(float(tok[1]))
                cur.press.append(float(tok[2]))
                cur.pot_e.append(float(tok[3]))
                cur.kin_e.append(float(tok[4]))
                cur.lx.append(float(tok[5]))
            except ValueError:
                continue
    return stages


# ── Reporting ─────────────────────────────────────────────────────


def text_summary(lammps: list[Stage], atv: list[Stage]) -> None:
    n_stages = max(len(lammps), len(atv))
    for i in range(n_stages):
        L = lammps[i] if i < len(lammps) else None
        A = atv[i]    if i < len(atv)    else None
        print(f"\n── Stage {i + 1} ─────────────────────────────────")
        print(f"  LAMMPS   : {L.label if L else '(no log)'}"
              f"  -- {len(L.step) if L else 0} samples")
        print(f"  Atomiverse: {A.label if A else '(no log)'}"
              f"  -- {len(A.step) if A else 0} samples")

        if L:
            t_mean = sum(L.temp) / len(L.temp)
            p_mean = sum(L.press) / len(L.press)
            print(f"     LAMMPS  ⟨T⟩={t_mean:9.2f} K   ⟨P⟩={p_mean:13.3e} atm")
        if A:
            t_mean = sum(A.temp) / len(A.temp)
            p_mean = sum(A.press) / len(A.press)
            print(f"     Atomivse⟨T⟩={t_mean:9.2f} K   ⟨P⟩={p_mean:13.3e} atm")


def plot_overlay(lammps: list[Stage], atv: list[Stage], out: Path) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[compare] matplotlib not installed — skipping plot.")
        return

    n_stages = max(len(lammps), len(atv))
    fig, axes = plt.subplots(3, n_stages, figsize=(6 * n_stages, 9),
                              sharex="col", squeeze=False)

    for i in range(n_stages):
        L = lammps[i] if i < len(lammps) else None
        A = atv[i]    if i < len(atv)    else None
        col = axes[:, i]
        col[0].set_title(f"Stage {i + 1}")
        if L:
            col[0].plot(L.step, L.temp, label="LAMMPS", color="C0")
            col[1].plot(L.step, L.press, label="LAMMPS", color="C0")
        if A:
            col[0].plot(A.step, A.temp, label="Atomiverse", color="C1")
            col[1].plot(A.step, A.press, label="Atomiverse", color="C1")
            col[2].plot(A.step, A.lx, label="Atomiverse Lx", color="C1")
        col[0].set_ylabel("T [K]")
        col[1].set_ylabel("P [atm]")
        col[2].set_ylabel("Lx [Å]")
        col[2].set_xlabel("step")
        for ax in col:
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=8)
        # Symlog for press if it spans many orders of magnitude.
        if A and any(abs(x) > 1e8 for x in A.press):
            col[1].set_yscale("symlog", linthresh=1e3)
        if L and any(abs(x) > 1e8 for x in L.press):
            col[1].set_yscale("symlog", linthresh=1e3)

    fig.tight_layout()
    fig.savefig(out, dpi=120)
    print(f"[compare] wrote {out}")


# ── main ──────────────────────────────────────────────────────────


def main() -> int:
    here = Path(__file__).parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--lammps",
                    default=str(here / "log.PE_nc10_cl1000.txt"))
    ap.add_argument("--atomiverse",
                    default=str(here / "log.atomiverse.txt"))
    ap.add_argument("--out",
                    default=str(here / "compare.png"))
    args = ap.parse_args()

    lammps = parse_lammps(Path(args.lammps))
    atv    = parse_atomiverse(Path(args.atomiverse))

    print(f"LAMMPS    : {args.lammps}  -> {len(lammps)} stages")
    print(f"Atomiverse: {args.atomiverse}  -> {len(atv)} stages")

    text_summary(lammps, atv)
    plot_overlay(lammps, atv, Path(args.out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
