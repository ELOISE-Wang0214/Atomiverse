"""``python -m atomiverse`` Typer entry — scheme B install CLI.

Modules are **positional names** (not ``--md k:v``)::

    python -m atomiverse install md --mpi
    python -m atomiverse install md pw --real-storage double
    python -m atomiverse install --config path.toml

``--config`` is mutually exclusive with CLI modules/options.

This module must **not** import ``atomiverse._binding`` or ``_atomiverse``
at import time (pure-Python wheel / soft-extension contract).
"""

from __future__ import annotations

import sys
from enum import Enum
from pathlib import Path
from typing import List, Optional, Sequence

import typer

from atomiverse.cli.config import (
    MODULE_TO_DEFINE,
    ConfigError,
    InstallConfig,
    apply_cli_modules,
    load_toml_config,
    merge_unbound_cli,
    validate_mode,
)
from atomiverse.cli.install import run_install
from atomiverse.cli.source import resolve_source
from atomiverse.cli.status import run_status

app = typer.Typer(
    name="atomiverse",
    help=(
        "Atomiverse on-demand native install. "
        "Pure-Python wheels ship without the C++ extension; "
        "use `install` to compile selected modules into site-packages."
    ),
    add_completion=False,
    no_args_is_help=False,
    invoke_without_command=True,
)


class RealPrecision(str, Enum):
    """Storage / work precision tier (CMake ``ATV_REAL_*``)."""

    float = "float"
    double = "double"


class IdSize(str, Enum):
    """Integer-ID width (CMake ``ATV_ID_SIZE``)."""

    bigbig = "bigbig"
    smallbig = "smallbig"


def _print_banner() -> None:
    try:
        from atomiverse import __version__
    except Exception:  # pragma: no cover
        __version__ = "unknown"
    print(f"atomiverse {__version__}")
    print()
    print("On-demand native install (scheme B — no k:v strings):")
    print("  python -m atomiverse install md")
    print("  python -m atomiverse install md pw --mpi --real-storage double")
    print("  python -m atomiverse install --config path.toml")
    print("  python -m atomiverse status")
    print()
    print("Run `python -m atomiverse --help` for full options.")


@app.callback(invoke_without_command=True)
def _root(
    ctx: typer.Context,
    version: bool = typer.Option(
        False,
        "--version",
        help="print the installed atomiverse package version and exit",
    ),
) -> None:
    """Atomiverse CLI root."""
    if version:
        from atomiverse.cli.source import package_version

        try:
            typer.echo(package_version())
        except ConfigError:
            from atomiverse import __version__

            typer.echo(__version__)
        raise typer.Exit(0)
    if ctx.invoked_subcommand is None:
        _print_banner()
        raise typer.Exit(0)


@app.command("install")
def install_cmd(
    modules: Optional[List[str]] = typer.Argument(
        None,
        help=(
            "product modules to enable (scheme B): md, minimize, mc, evaluate, "
            "vib, path, pw, gpw, am1, cuda, kspace, potentials, neighbor"
        ),
        metavar="MODULE",
    ),
    config: Optional[Path] = typer.Option(
        None,
        "--config",
        help="unified TOML config (XOR with CLI modules/options)",
        exists=False,
        dir_okay=False,
        readable=True,
    ),
    mpi: Optional[bool] = typer.Option(
        None,
        "--mpi/--no-mpi",
        help="global: ATV_USE_MPI",
    ),
    real_storage: Optional[RealPrecision] = typer.Option(
        None,
        "--real-storage",
        help="global: ATV_REAL_STORAGE (float|double)",
    ),
    real_work: Optional[RealPrecision] = typer.Option(
        None,
        "--real-work",
        help="global: ATV_REAL_WORK (float|double)",
    ),
    id_size: Optional[IdSize] = typer.Option(
        None,
        "--id-size",
        help="global: ATV_ID_SIZE (bigbig|smallbig)",
    ),
    hf: bool = typer.Option(
        False,
        "--hf",
        help="(rejected) HF/SCF bindings are not available yet",
    ),
    src: Optional[Path] = typer.Option(
        None,
        "--src",
        help="local Atomiverse checkout (python/pyproject version must match)",
        exists=False,
        file_okay=False,
        dir_okay=True,
    ),
    url: Optional[str] = typer.Option(
        None,
        "--url",
        help="source tarball/URL override (default: GitHub tag of installed version)",
    ),
    define: Optional[List[str]] = typer.Option(
        None,
        "-D",
        "--define",
        help="unbound CMake cache entry KEY=VALUE (bound keys are rejected)",
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        help="resolve modules/defines and print the pip command without building",
    ),
) -> None:
    """Build/install the native ``_atomiverse`` extension for selected modules."""
    mod_list = list(modules or [])
    # Reject unknown module tokens early with a clear message.
    for m in mod_list:
        if m not in MODULE_TO_DEFINE and m != "hf":
            typer.echo(
                f"error: unknown module {m!r}; known: {sorted(MODULE_TO_DEFINE)}",
                err=True,
            )
            raise typer.Exit(2)

    cli_options = {
        "mpi": mpi,
        "real_storage": real_storage.value if real_storage is not None else None,
        "real_work": real_work.value if real_work is not None else None,
        "id_size": id_size.value if id_size is not None else None,
    }
    has_cli_options = any(v is not None for v in cli_options.values())

    try:
        validate_mode(
            config_path=str(config) if config is not None else None,
            modules=mod_list,
            hf=hf,
            has_cli_options=has_cli_options,
        )
        if config is not None:
            cfg = load_toml_config(config)
        else:
            cfg = InstallConfig()
            apply_cli_modules(cfg, mod_list, options=cli_options)
        if define:
            merge_unbound_cli(cfg, define)
        source = resolve_source(
            src=str(src) if src is not None else None,
            url=url,
        )
    except ConfigError as exc:
        typer.echo(f"error: {exc}", err=True)
        raise typer.Exit(2) from exc

    try:
        code = run_install(cfg, source, dry_run=dry_run)
    except ConfigError as exc:
        typer.echo(f"error: {exc}", err=True)
        raise typer.Exit(2) from exc
    if code != 0:
        raise typer.Exit(code)


@app.command("status")
def status_cmd() -> None:
    """Show enabled modules and define fingerprint from the last install."""
    code = run_status()
    if code != 0:
        raise typer.Exit(code)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """CLI entry point compatible with tests and ``python -m atomiverse``.

    @param argv: Arguments without the program name (default: ``sys.argv[1:]``).
    @return Process exit code.
    """
    args = list(sys.argv[1:] if argv is None else argv)
    try:
        # standalone_mode=False: Click may *return* the exit code or raise Exit.
        result = app(args=args, standalone_mode=False)
        if isinstance(result, int):
            return result
        return 0
    except typer.Exit as exc:
        code = getattr(exc, "exit_code", None)
        if code is None:
            code = getattr(exc, "code", 0)
        return int(code) if code is not None else 0
    except SystemExit as exc:  # pragma: no cover
        code = exc.code
        if code is None:
            return 0
        return int(code) if isinstance(code, int) else 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
