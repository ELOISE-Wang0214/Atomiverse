"""``python -m atomiverse`` argparse entry — ``install`` / ``status``.

This module must **not** import ``atomiverse._binding`` or ``_atomiverse``
at import time (pure-Python wheel / soft-extension contract).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional, Sequence

from atomiverse.cli.config import (
    MODULE_TO_DEFINE,
    ConfigError,
    InstallConfig,
    apply_inline_modules,
    load_toml_config,
    merge_unbound_cli,
    validate_mode,
)
from atomiverse.cli.install import run_install
from atomiverse.cli.source import resolve_source
from atomiverse.cli.status import run_status

# Product-module CLI flags (excluding the special --global / --hf).
_MODULE_FLAGS: tuple[str, ...] = tuple(
    m for m in MODULE_TO_DEFINE if m != "global"
)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m atomiverse",
        description=(
            "Atomiverse on-demand native install. "
            "Pure-Python wheels ship without the C++ extension; "
            "use `install` to compile modules into site-packages."
        ),
    )
    parser.add_argument(
        "--version",
        action="store_true",
        help="print the installed atomiverse package version and exit",
    )

    sub = parser.add_subparsers(dest="command")

    # ── install ──────────────────────────────────────────────────────────
    inst = sub.add_parser(
        "install",
        help="build/install the native _atomiverse extension for selected modules",
    )
    inst.add_argument(
        "--config",
        metavar="PATH.toml",
        default=None,
        help="unified TOML config (XOR with inline --md/--pw/… flags)",
    )
    for name in _MODULE_FLAGS:
        inst.add_argument(
            f"--{name}",
            nargs="?",
            const="",
            default=None,
            metavar="k:v,...",
            help=f"enable module '{name}' with optional inline options",
        )
    inst.add_argument(
        "--global",
        dest="global_opts",
        nargs="?",
        const="",
        default=None,
        metavar="k:v,...",
        help="cross-module options (real_storage, mpi, …)",
    )
    inst.add_argument(
        "--hf",
        action="store_true",
        help="(rejected) HF/SCF bindings are not available yet",
    )
    inst.add_argument(
        "--src",
        default=None,
        metavar="PATH",
        help="local Atomiverse checkout (python/pyproject version must match)",
    )
    inst.add_argument(
        "--url",
        default=None,
        metavar="URL",
        help="source tarball/URL override (default: GitHub tag of installed version)",
    )
    inst.add_argument(
        "-D",
        dest="defines",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="unbound CMake cache entry (bound keys are rejected)",
    )
    inst.add_argument(
        "--dry-run",
        action="store_true",
        help="resolve modules/defines and print the pip command without building",
    )

    # ── status ───────────────────────────────────────────────────────────
    sub.add_parser(
        "status",
        help="show enabled modules and define fingerprint from last install",
    )

    return parser


def _module_flags_from_args(args: argparse.Namespace) -> dict:
    return {name: getattr(args, name.replace("-", "_"), None) for name in _MODULE_FLAGS}


def _cmd_install(args: argparse.Namespace) -> int:
    module_flags = _module_flags_from_args(args)
    try:
        validate_mode(
            config_path=args.config,
            module_flags=module_flags,
            global_spec=args.global_opts,
            hf=bool(args.hf),
        )
        if args.config:
            cfg = load_toml_config(Path(args.config))
        else:
            cfg = InstallConfig()
            apply_inline_modules(cfg, module_flags, global_spec=args.global_opts)
        if args.defines:
            merge_unbound_cli(cfg, args.defines)
        source = resolve_source(src=args.src, url=args.url)
    except ConfigError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    try:
        return run_install(cfg, source, dry_run=bool(args.dry_run))
    except ConfigError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


def _cmd_status(_args: argparse.Namespace) -> int:
    return run_status()


def main(argv: Optional[Sequence[str]] = None) -> int:
    """CLI entry point.

    @param argv: Arguments without the program name (default: ``sys.argv[1:]``).
    @return Process exit code.
    """
    argv_list: List[str] = list(sys.argv[1:] if argv is None else argv)
    parser = _build_parser()

    # Bare invocation: friendly guidance (pure-wheel contract / pure_import).
    if not argv_list:
        try:
            from atomiverse import __version__
        except Exception:  # pragma: no cover
            __version__ = "unknown"
        print(f"atomiverse {__version__}")
        print()
        print("On-demand native install:")
        print("  python -m atomiverse install --md")
        print("  python -m atomiverse install --config path.toml")
        print("  python -m atomiverse status")
        print()
        print("Run `python -m atomiverse --help` for full options.")
        return 0

    # `python -m atomiverse --version` (package version — not an install flag)
    if argv_list == ["--version"] or argv_list == ["-v"]:
        from atomiverse.cli.source import package_version

        try:
            print(package_version())
        except ConfigError:
            from atomiverse import __version__

            print(__version__)
        return 0

    args = parser.parse_args(argv_list)

    if getattr(args, "version", False) and not args.command:
        from atomiverse.cli.source import package_version

        try:
            print(package_version())
        except ConfigError:
            from atomiverse import __version__

            print(__version__)
        return 0

    if args.command == "install":
        return _cmd_install(args)
    if args.command == "status":
        return _cmd_status(args)

    parser.print_help()
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
