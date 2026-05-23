#!/usr/bin/env python3
"""Canonicalize raw MSU-PE engine output into the downstream canonical CSV schema.

This is the Phase D.1 canonicalization step for
``.claude/specs/msu-pe-phase-d-01-reference-runs.md``. It is pure stdlib Python
(``csv`` + ``argparse`` only) so it runs on a fresh install, is side-effect-free
beyond writing the requested output file, and is idempotent: running it twice on
the same input yields byte-identical output.

"Canonicalization" = converting a raw engine output file into the schema the
downstream chain links (``msu-pe-align-03`` / ``-d-03`` / ``-d-04``) consume:

* thermo        — ``step,T,P,E_pot,E_kin,Pxx,Pyy,Pzz``
* stress-strain — ``strain,sigma_xx``

Engine / kind matrix
--------------------

* ``lammps stress-strain`` — the input is the upstream paper deck's
  space-separated, ``#``-commented ``def1`` print file (columns
  ``strain pxx_GPa pyy_GPa pzz_GPa lx ly lz temp``, the deck's
  ``-pxx/10000*1.01325`` already converts atm→GPa and flips the sign to
  tension-positive). Output ``strain,sigma_xx`` is def1 col0, col1.
* ``lammps thermo`` — the online ``lammps.in`` ``fix print`` already streams the
  canonical 8-column schema. This mode validates the header / column order and
  passes the rows through unchanged; a non-canonical header fails loudly.
* ``atomiverse thermo`` — the benchmark's ``thermo.csv`` is remapped exactly as
  the inline ``awk`` in ``tools/run_msu_atomiverse.sh`` did: ``$1,$4,$5,$6,$7``
  become ``step,T,P,E_pot,E_kin`` and the isotropic ``P`` ($5) is replicated
  into ``Pxx/Pyy/Pzz``. (Pre msu-pe-cuda-components-01: Phase G.2 stand-in. Real per-axis
  pressure is an umbrella-phase fix) and is preserved verbatim here.
* ``atomiverse stress-strain`` — the benchmark's ``stress_strain.csv`` is
  remapped as the inline ``awk`` did: ``$3 → strain``, ``$5 → sigma_xx``.

Behaviour
---------

* ``#``-comment lines and blank lines in the input are skipped.
* An empty or header-only input yields a header-only output, not a crash.
* A malformed row (wrong column count) raises :class:`CanonicalizeError`
  naming the offending line; the CLI turns that into a non-zero exit.

CLI
---

::

    canonicalize_msu_output.py <engine> <kind> <in> <out>
        engine ∈ {lammps, atomiverse}
        kind   ∈ {thermo, stress-strain}

Exit codes
----------

* ``0`` — conversion succeeded.
* ``1`` — a malformed input row or a non-canonical header was found.
* ``2`` — a bad CLI argument (handled by argparse) or a missing input file.

@see .claude/specs/msu-pe-phase-d-01-reference-runs.md
"""

import argparse
import csv
import os
import re
import sys


#: Canonical thermo CSV schema consumed by the downstream chain links.
#: Per CLAUDE.md project units: T [K], P [bar], E_pot/E_kin [eV],
#: Pxx/Pyy/Pzz [bar]. Atomiverse already streams in these units;
#: LAMMPS ``units real`` is converted on canonicalisation.
THERMO_HEADER = ["step", "T", "P", "E_pot", "E_kin", "Pxx", "Pyy", "Pzz"]

#: Canonical stress-strain CSV schema consumed by msu-pe-phase-d-04.
STRESS_HEADER = ["strain", "sigma_xx"]

#: LAMMPS ``units real`` → canonical unit conversions.
#:   pressure: 1 atm     = 1.01325 bar
#:   energy:   1 kcal/mol = 1/23.060547830619026 eV
_ATM_TO_BAR = 1.01325
_KCAL_PER_MOL_TO_EV = 1.0 / 23.060547830619026


class CanonicalizeError(Exception):
    """A raw engine output file could not be canonicalized.

    Raised on a malformed input row (wrong column count) or a non-canonical
    LAMMPS thermo header. The message names the offending line so the failure
    is loud, not silent.
    """


# ── input parsing ────────────────────────────────────────────────────────────


def _iter_data_lines(path):
    """Yield ``(line_number, raw_line)`` for non-blank, non-``#``-comment lines.

    ``line_number`` is 1-based and counts every physical line of the file so the
    number in an error message points at the real line. Blank lines and lines
    whose first non-whitespace character is ``#`` are skipped.

    @raises CanonicalizeError  When the input file does not exist.
    """
    if not os.path.isfile(path):
        raise CanonicalizeError("input file not found: {}".format(path))
    with open(path, newline="") as f:
        for lineno, raw in enumerate(f, start=1):
            stripped = raw.strip()
            if not stripped or stripped.startswith("#"):
                continue
            yield lineno, raw


def _split_ws(raw):
    """Split a whitespace-separated row into its fields (no empty fields)."""
    return raw.split()


def _split_csv(raw):
    """Split a comma-separated row into its fields, trimming surrounding space."""
    return [cell.strip() for cell in next(csv.reader([raw]))]


# ── per-engine / per-kind converters ─────────────────────────────────────────


def _lammps_stress_strain(path):
    """``lammps stress-strain``: def1 print → ``strain,sigma_xx``.

    The input is the upstream paper deck's space-separated, ``#``-commented
    ``def1`` file with 8 columns ``strain pxx_GPa pyy_GPa pzz_GPa lx ly lz
    temp``. The canonical row is ``(col0, col1)`` — ``sigma_xx`` is the deck's
    already-GPa, already-tension-positive ``pxx``.

    The online ``lammps.in`` ``fix ss5`` already streams the canonical
    comma-separated ``strain,sigma_xx`` schema. When the first data line is
    exactly that header this mode validates the column order and passes the
    rows through unchanged, so the online and offline LAMMPS stress-strain
    paths share one tested canonicalization implementation.
    """
    rows = []
    expected = 8
    canonical_mode = False
    header_seen = False
    for lineno, raw in _iter_data_lines(path):
        if not header_seen:
            header_seen = True
            if _split_csv(raw) == STRESS_HEADER:
                canonical_mode = True
                continue
            # else: not the online CSV — fall through, treat as a def1 row.
        if canonical_mode:
            fields = _split_csv(raw)
            if len(fields) != len(STRESS_HEADER):
                raise CanonicalizeError(
                    "malformed stress-strain row at line {}: expected {} "
                    "columns, got {} -- {!r}".format(
                        lineno, len(STRESS_HEADER), len(fields), raw.strip()))
            rows.append(fields)
            continue
        fields = _split_ws(raw)
        if len(fields) != expected:
            raise CanonicalizeError(
                "malformed def1 row at line {}: expected {} columns, got {} "
                "-- {!r}".format(lineno, expected, len(fields), raw.strip()))
        rows.append([fields[0], fields[1]])
    return STRESS_HEADER, rows


def _lammps_thermo(path):
    """``lammps thermo``: convert ``units real`` to canonical units + dedupe.

    The online ``lammps.in`` ``fix print`` emits the canonical 8-column
    header ``step,T,P,E_pot,E_kin,Pxx,Pyy,Pzz`` but the values are in
    ``units real`` (T [K], P [atm], E [kcal/mol]). This pass converts
    P / Pxx / Pyy / Pzz from atm to bar (×1.01325) and E_pot / E_kin from
    kcal/mol to eV (÷23.060547830619026), so the canonical CSV matches the
    Atomiverse-emitted unit set and ``validate_polyethylene.py`` can compare
    ensemble means in a single unit.

    Consecutive rows with the same ``step`` value (LAMMPS sometimes emits a
    duplicate at stage transitions where one ``fix print`` rolls off and
    another rolls on at the same timestep) are deduplicated — only the first
    is kept. This satisfies the strict-monotonic-step contract that
    ``tests/cpu/md/test_msu_reference_inputs.cpp`` enforces.

    A non-canonical header still fails loudly; a malformed data row still
    raises :class:`CanonicalizeError` naming the line.
    """
    header = None
    rows = []
    expected = len(THERMO_HEADER)
    last_step = None
    for lineno, raw in _iter_data_lines(path):
        fields = _split_csv(raw)
        if header is None:
            header = fields
            if header != THERMO_HEADER:
                raise CanonicalizeError(
                    "non-canonical LAMMPS thermo header at line {}: expected "
                    "{}, got {}".format(lineno, THERMO_HEADER, header))
            continue
        if len(fields) != expected:
            raise CanonicalizeError(
                "malformed thermo row at line {}: expected {} columns, got {} "
                "-- {!r}".format(lineno, expected, len(fields), raw.strip()))
        if fields[0] == last_step:
            continue  # skip consecutive-same-step duplicate
        last_step = fields[0]
        try:
            step = fields[0]
            t_K = fields[1]                                    # K — unchanged
            p_bar = repr(float(fields[2]) * _ATM_TO_BAR)
            e_pot_eV = repr(float(fields[3]) * _KCAL_PER_MOL_TO_EV)
            e_kin_eV = repr(float(fields[4]) * _KCAL_PER_MOL_TO_EV)
            pxx_bar = repr(float(fields[5]) * _ATM_TO_BAR)
            pyy_bar = repr(float(fields[6]) * _ATM_TO_BAR)
            pzz_bar = repr(float(fields[7]) * _ATM_TO_BAR)
        except ValueError as exc:
            raise CanonicalizeError(
                "non-numeric thermo cell at line {}: {!r} -- {}".format(
                    lineno, raw.strip(), exc))
        rows.append([step, t_K, p_bar, e_pot_eV, e_kin_eV,
                     pxx_bar, pyy_bar, pzz_bar])
    return THERMO_HEADER, rows


def _atomiverse_thermo(path):
    """``atomiverse thermo``: benchmark ``thermo.csv`` → canonical 8 columns.

    Reads the Phase G.2 + G.3 17-column extended schema produced by
    ``benchmarks/md/polyethylene_tension.cpp``::

        step,stage,T_target_K,T_K,P_bar,Pxx,Pyy,Pzz,E_pot_eV,E_pair_eV,
        E_bond_eV,E_angle_eV,E_dihedral_eV,E_kin_eV,box_x,box_y,box_z

    and emits the canonical 8-column schema ``step,T,P,E_pot,E_kin,Pxx,Pyy,Pzz``
    by passthrough — no more isotropic-P replication into Pxx/Pyy/Pzz (that
    was the Phase G.2 stand-in, retired by msu-pe-cuda-components-01).
    The 4 per-term energy columns (E_pair / E_bond / E_angle / E_dihedral)
    are read by the benchmark output but not propagated to the canonical
    schema in this link; downstream consumers that need the decomposition
    read the raw thermo.csv directly.
    """
    header = None
    rows = []
    min_cols = 14  # need at least up to E_kin_eV (column index 13)
    for lineno, raw in _iter_data_lines(path):
        fields = _split_csv(raw)
        if header is None:
            header = fields  # benchmark header, not validated column-by-column
            continue
        if len(fields) < min_cols:
            raise CanonicalizeError(
                "malformed atomiverse thermo row at line {}: expected at least "
                "{} columns, got {} -- {!r}".format(
                    lineno, min_cols, len(fields), raw.strip()))
        step = fields[0]
        t = fields[3]       # T_K
        p = fields[4]       # P_bar
        pxx = fields[5]     # Pxx (real per-axis virial, bar)
        pyy = fields[6]     # Pyy
        pzz = fields[7]     # Pzz
        e_pot = fields[8]   # E_pot_eV
        e_kin = fields[13]  # E_kin_eV
        rows.append([step, t, p, e_pot, e_kin, pxx, pyy, pzz])
    return THERMO_HEADER, rows


def _atomiverse_stress_strain(path):
    """``atomiverse stress-strain``: benchmark ``stress_strain.csv`` → 2 columns.

    Reproduces the inline ``awk`` of ``tools/run_msu_atomiverse.sh``:
    ``$3 → strain``, ``$5 → sigma_xx`` (the ``NF >= 5`` guard).
    """
    header = None
    rows = []
    min_cols = 5  # awk guard: NF >= 5
    for lineno, raw in _iter_data_lines(path):
        fields = _split_csv(raw)
        if header is None:
            header = fields
            continue
        if len(fields) < min_cols:
            raise CanonicalizeError(
                "malformed atomiverse stress-strain row at line {}: expected "
                "at least {} columns, got {} -- {!r}".format(
                    lineno, min_cols, len(fields), raw.strip()))
        rows.append([fields[2], fields[4]])
    return STRESS_HEADER, rows


#: Dispatch table: (engine, kind) → converter returning ``(header, rows)``.
_CONVERTERS = {
    ("lammps", "stress-strain"): _lammps_stress_strain,
    ("lammps", "thermo"): _lammps_thermo,
    ("atomiverse", "thermo"): _atomiverse_thermo,
    ("atomiverse", "stress-strain"): _atomiverse_stress_strain,
}


# ── msu-pe-cuda-components-08 — backend tag ──────────────────────────────────

#: Backend values recognised by the CLI.
_BACKENDS = ("cpu", "cuda", "auto")


def infer_backend_from_path(path):
    """Infer the backend (cpu / cuda) from a raw-output directory path.

    Convention (msu-pe-cuda-components-08 § Design):

    * ``runs/atomiverse_cuda_*`` (any segment of the path) → ``cuda``
    * ``runs/lammps_*_gpu_*``                              → ``cuda``
    * anything else (incl. pre-link-08 CPU-era paths)      → ``cpu``

    The mapping is intentionally permissive: pre-existing CPU-era runs
    without a backend tag in their path map to ``cpu`` so back-compat is
    preserved (no breaking change for the existing
    ``thermo_atomiverse_1r.csv`` / ``thermo_lammps_1r.csv`` artefacts).

    @param path  Any file path. Only the path segments are inspected; the
                 file content is not read.
    @return  ``"cuda"`` or ``"cpu"``.
    """
    norm = path.replace("\\", "/")
    # atomiverse_cuda_* — explicit CUDA tag in any segment.
    if re.search(r"(^|/)atomiverse_cuda(_|/)", norm):
        return "cuda"
    # lammps_*_gpu_* — GPU mode tag in any segment.
    if re.search(r"(^|/)lammps_[^/]*_gpu(_|/)", norm):
        return "cuda"
    return "cpu"


def _detect_kind_from_path(path):
    """Best-effort kind inference from a raw-output filename.

    Used by the flag-based CLI form when ``--kind`` is omitted:

    * filenames containing ``stress_strain`` (or ``stress-strain``) → ``stress-strain``
    * filenames containing ``thermo``                               → ``thermo``
    * otherwise                                                     → ``None`` (CLI demands explicit ``--kind``)
    """
    name = os.path.basename(path).lower()
    if "stress_strain" in name or "stress-strain" in name:
        return "stress-strain"
    if "thermo" in name:
        return "thermo"
    return None


# ── public API ───────────────────────────────────────────────────────────────


def canonicalize(engine, kind, in_path, out_path, metadata=None):
    """Convert a raw engine output file into a canonical-schema CSV.

    @param engine    ``"lammps"`` or ``"atomiverse"``.
    @param kind      ``"thermo"`` or ``"stress-strain"``.
    @param in_path   Raw engine output file (def1 print / benchmark CSV /
                     online LAMMPS canonical CSV).
    @param out_path  Canonical CSV to write. An empty / header-only input
                     yields a header-only output (metadata comment line, if
                     requested, still emitted).
    @param metadata  Optional ``dict`` of provenance fields written as a
                     leading ``#`` comment line so downstream readers
                     (``_iter_data_lines`` skips ``#`` lines verbatim) can
                     surface engine / backend / seed without re-parsing the
                     run directory. Keys are emitted in deterministic order:
                     ``engine``, ``backend``, ``seed``, then everything else
                     alphabetically. ``None`` (default) keeps the legacy
                     header-only output untouched for back-compat.
    @raises CanonicalizeError  On an unknown ``(engine, kind)`` pair, a missing
                               input file, a malformed input row, or a
                               non-canonical LAMMPS thermo header.
    """
    converter = _CONVERTERS.get((engine, kind))
    if converter is None:
        raise CanonicalizeError(
            "unknown engine/kind combination: {!r}/{!r}".format(engine, kind))
    header, rows = converter(in_path)
    with open(out_path, "w", newline="") as f:
        if metadata:
            # Deterministic field order: engine, backend, seed first; then
            # any remaining keys alphabetically. Comma-separated, leading
            # ``#`` so the canonical reader skips it transparently.
            ordered_keys = []
            for k in ("engine", "backend", "seed"):
                if k in metadata:
                    ordered_keys.append(k)
            for k in sorted(metadata):
                if k not in ordered_keys:
                    ordered_keys.append(k)
            parts = ["{}={}".format(k, metadata[k]) for k in ordered_keys]
            f.write("# " + ",".join(parts) + "\n")
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def _backend_tagged_output_name(engine, backend, kind, runs):
    """Filename for the canonical output under the flag-based CLI form.

    Pattern (msu-pe-cuda-components-08 § Design):

    * ``thermo``        → ``thermo_<engine>_<backend>_<runs>r.csv``
    * ``stress-strain`` → ``stress_strain_<engine>_<backend>_<runs>r.csv``

    Back-compat note: pre-link-08 names ``thermo_atomiverse_1r.csv`` (no
    backend tag) are *not* re-emitted by this helper — they continue to be
    produced by the legacy positional CLI form / by the older runner
    scripts. The new tagged names live alongside them; readers that want
    backend-discriminated artefacts opt in by using the new names.
    """
    stem = "thermo" if kind == "thermo" else "stress_strain"
    return "{}_{}_{}_{}r.csv".format(stem, engine, backend, runs)


# ── CLI ──────────────────────────────────────────────────────────────────────


def _is_legacy_positional_form(argv):
    """Heuristic: legacy form is ``<engine> <kind> <in> <out>``.

    Detected by the first non-flag argument being one of the engine choices
    AND there being at least four non-flag tokens. The new flag-based form
    always has ``--engine`` somewhere and at most one positional (the input
    file).
    """
    positionals = [a for a in argv if not a.startswith("-")]
    if len(positionals) < 4:
        return False
    return positionals[0] in ("lammps", "atomiverse")


def main(argv=None):
    """CLI entry point. Two supported invocation forms:

    Legacy (positional)::

        canonicalize.py <engine> <kind> <input> <output>

    Flag-based (msu-pe-cuda-components-08)::

        canonicalize.py <input> --engine {lammps,atomiverse}
                                [--kind {thermo,stress-strain}]
                                [--backend {cpu,cuda,auto}]
                                [--seed <S>]
                                [--runs <N>]
                                [--output-dir <D>]
                                [--output <FILE>]

    ``--backend auto`` (default in the flag form) infers backend from the
    input path via :func:`infer_backend_from_path`. ``--kind`` is inferred
    from the filename when omitted (``thermo`` / ``stress_strain``). When
    ``--output-dir`` is given, the output filename is
    ``thermo_<engine>_<backend>_<runs>r.csv`` (or ``stress_strain_…``); pass
    ``--output`` to override with a concrete path.

    See the module docstring for the exit-code contract.
    """
    if argv is None:
        argv = sys.argv[1:]

    # Legacy positional form: <engine> <kind> <in> <out>.
    if _is_legacy_positional_form(argv):
        p = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter,
        )
        p.add_argument("engine", choices=["lammps", "atomiverse"],
                       help="Source engine of the raw output file.")
        p.add_argument("kind", choices=["thermo", "stress-strain"],
                       help="Kind of canonical schema to emit.")
        p.add_argument("input", help="Raw engine output file to canonicalize.")
        p.add_argument("output", help="Canonical-schema CSV file to write.")
        args = p.parse_args(argv)
        try:
            canonicalize(args.engine, args.kind, args.input, args.output)
        except CanonicalizeError as exc:
            print("canonicalize_msu_output: {}".format(exc), file=sys.stderr)
            return 1
        return 0

    # Flag-based form: <input> --engine ... [--backend ...] [--output-dir ...]
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("input", help="Raw engine output file to canonicalize.")
    p.add_argument("--engine", required=True,
                   choices=["lammps", "atomiverse"],
                   help="Source engine of the raw output file.")
    p.add_argument("--kind", choices=["thermo", "stress-strain"],
                   default=None,
                   help="Kind of canonical schema to emit "
                        "(default: inferred from the input filename).")
    p.add_argument("--backend", choices=list(_BACKENDS), default="auto",
                   help="Backend tag for the metadata comment + filename "
                        "suffix (default: auto, infer from input path).")
    p.add_argument("--seed", default=None,
                   help="Provenance seed; written into the metadata comment.")
    p.add_argument("--runs", default="1",
                   help="Number of ranks / runs; used in the tagged "
                        "filename `<engine>_<backend>_<runs>r.csv` "
                        "(default: 1).")
    p.add_argument("--output-dir", default=None,
                   help="Directory under which the tagged filename is "
                        "written. Mutually exclusive with --output.")
    p.add_argument("--output", default=None,
                   help="Concrete output path (overrides --output-dir + "
                        "tagged-filename logic).")
    args = p.parse_args(argv)

    if args.output_dir and args.output:
        print("canonicalize_msu_output: --output-dir and --output are "
              "mutually exclusive", file=sys.stderr)
        return 2

    kind = args.kind or _detect_kind_from_path(args.input)
    if kind is None:
        print("canonicalize_msu_output: could not infer --kind from "
              "filename {!r}; pass --kind explicitly".format(args.input),
              file=sys.stderr)
        return 2

    backend = args.backend
    if backend == "auto":
        backend = infer_backend_from_path(args.input)

    if args.output:
        out_path = args.output
    else:
        out_dir = args.output_dir or os.path.dirname(args.input) or "."
        out_path = os.path.join(
            out_dir,
            _backend_tagged_output_name(args.engine, backend, kind, args.runs))

    metadata = {"engine": args.engine, "backend": backend}
    if args.seed is not None:
        metadata["seed"] = args.seed

    try:
        canonicalize(args.engine, kind, args.input, out_path,
                     metadata=metadata)
    except CanonicalizeError as exc:
        print("canonicalize_msu_output: {}".format(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
