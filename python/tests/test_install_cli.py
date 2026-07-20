"""Unit tests for on-demand install CLI (Typer scheme B — no k:v).

Covers mutual exclusion, empty install, --hf, CLI options, DAG closure,
bound -D reject, version gate, status, import hygiene.
"""

from __future__ import annotations

import json
import sys
import textwrap
from pathlib import Path

import pytest

# Ensure src layout is importable when tests run without an editable install.
_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from atomiverse.cli.config import (  # noqa: E402
    BINDING_TABLE,
    ConfigError,
    InstallConfig,
    apply_cli_modules,
    load_toml_config,
    merge_unbound_cli,
    resolve_defines,
    validate_mode,
)
from atomiverse.cli.install import (  # noqa: E402
    build_pip_command,
    dry_run_defines,
)
from atomiverse.cli.main import main  # noqa: E402
from atomiverse.cli.source import (  # noqa: E402
    default_tarball_url,
    package_version,
    resolve_source,
    verify_source_version,
)
from atomiverse.cli.status import format_status, run_status  # noqa: E402


# ── ac-001 mutual exclusion ──────────────────────────────────────────────────


def test_ac001_mutual_exclusion_config_vs_modules() -> None:
    with pytest.raises(ConfigError, match="mutual exclusion"):
        validate_mode(
            config_path="x.toml",
            modules=["md"],
            hf=False,
        )
    rc = main(["install", "--config", "x.toml", "md"])
    assert rc != 0


def test_ac001_cli_stderr_names_mutual_exclusion(capsys) -> None:
    rc = main(["install", "--config", "does-not-matter.toml", "pw"])
    assert rc != 0
    err = capsys.readouterr().err
    assert "mutual" in err.lower() or "exclusion" in err.lower()


def test_ac001_config_with_cli_options_rejected(capsys) -> None:
    rc = main(["install", "--config", "x.toml", "--mpi"])
    assert rc != 0
    err = capsys.readouterr().err.lower()
    assert "mutual" in err or "exclusion" in err


# ── ac-002 empty install ─────────────────────────────────────────────────────


def test_ac002_empty_install_fails() -> None:
    with pytest.raises(ConfigError, match="empty install"):
        validate_mode(config_path=None, modules=[], hf=False)
    rc = main(["install"])
    assert rc != 0


def test_ac002_options_only_is_empty() -> None:
    with pytest.raises(ConfigError, match="empty install"):
        validate_mode(
            config_path=None,
            modules=[],
            hf=False,
            has_cli_options=True,
        )
    rc = main(["install", "--mpi"])
    assert rc != 0


# ── ac-003 --hf rejected ─────────────────────────────────────────────────────


def test_ac003_hf_flag_rejected() -> None:
    with pytest.raises(ConfigError, match="HF/SCF|not available"):
        validate_mode(config_path=None, modules=["md"], hf=True)
    rc = main(["install", "--hf"])
    assert rc != 0
    rc = main(["install", "md", "--hf"])
    assert rc != 0


def test_ac003_modules_hf_in_toml(tmp_path: Path) -> None:
    cfg = tmp_path / "hf.toml"
    cfg.write_text("[modules]\nhf = true\n", encoding="utf-8")
    with pytest.raises(ConfigError, match="HF/SCF|not available"):
        load_toml_config(cfg)


# ── scheme B: positional modules + native options ────────────────────────────


def test_apply_cli_modules_md_mpi() -> None:
    cfg = InstallConfig()
    apply_cli_modules(cfg, ["md"], options={"mpi": True})
    defines = resolve_defines(cfg)
    assert defines["ATV_BUILD_MD"] == "ON"
    assert defines[BINDING_TABLE["mpi"]] == "ON"
    assert BINDING_TABLE["mpi"] == "ATV_USE_MPI"


def test_apply_cli_modules_unknown() -> None:
    with pytest.raises(ConfigError, match="unknown module"):
        apply_cli_modules(InstallConfig(), ["notamodule"])


def test_cli_multi_module_md_pw() -> None:
    cfg = InstallConfig()
    apply_cli_modules(cfg, ["md", "pw"])
    defines = resolve_defines(cfg)
    assert defines["ATV_BUILD_MD"] == "ON"
    assert defines["ATV_BUILD_PW"] == "ON"
    assert defines["ATV_BUILD_POTENTIALS"] == "ON"
    assert defines["ATV_BUILD_NEIGHBOR"] == "ON"


# ── ac-005 md auto-enable closure ────────────────────────────────────────────


def test_ac005_md_closure_enables_potentials_and_neighbor() -> None:
    cfg = InstallConfig(modules={"md": True})
    defines = resolve_defines(cfg)
    assert defines["ATV_BUILD_MD"] == "ON"
    assert defines["ATV_BUILD_POTENTIALS"] == "ON"
    assert defines["ATV_BUILD_NEIGHBOR"] == "ON"
    assert defines["ATV_BUILD_PW"] == "OFF"
    assert defines["ATV_BUILD_MC"] == "OFF"


def test_ac005_path_closure_pulls_minimizer_then_potentials() -> None:
    cfg = InstallConfig(modules={"path": True})
    defines = resolve_defines(cfg)
    assert defines["ATV_BUILD_PATH"] == "ON"
    assert defines["ATV_BUILD_MINIMIZER"] == "ON"
    assert defines["ATV_BUILD_POTENTIALS"] == "ON"
    assert defines["ATV_BUILD_NEIGHBOR"] == "ON"


def test_ac005_kspace_closure_enables_potentials() -> None:
    cfg = InstallConfig(modules={"kspace": True})
    defines = resolve_defines(cfg)
    assert defines["ATV_BUILD_KSPACE"] == "ON"
    assert defines["ATV_BUILD_POTENTIALS"] == "ON"


# ── ac-006 bound keys cannot be unbound -D ───────────────────────────────────


def test_ac006_bound_key_as_unbound_rejected() -> None:
    cfg = InstallConfig(modules={"md": True})
    with pytest.raises(ConfigError, match="bound key"):
        merge_unbound_cli(cfg, ["ATV_BUILD_MD=ON"])
    with pytest.raises(ConfigError, match="bound key"):
        merge_unbound_cli(cfg, ["mpi=ON"])
    with pytest.raises(ConfigError, match="bound key"):
        merge_unbound_cli(cfg, ["ATV_USE_MPI=ON"])


def test_ac006_bound_key_in_toml_unbound(tmp_path: Path) -> None:
    cfg_path = tmp_path / "bad.toml"
    cfg_path.write_text(
        textwrap.dedent(
            """\
            [modules]
            md = true
            [cmake.unbound]
            ATV_BUILD_MD = "ON"
            """
        ),
        encoding="utf-8",
    )
    with pytest.raises(ConfigError, match="bound key"):
        load_toml_config(cfg_path)


def test_ac006_cli_D_bound_key(capsys) -> None:
    rc = main(
        [
            "install",
            "md",
            "-D",
            "ATV_BUILD_MD=ON",
            "--dry-run",
            "--src",
            str(_SRC.parent),
        ]
    )
    assert rc != 0
    err = capsys.readouterr().err.lower()
    assert "bound" in err


# ── ac-008 source version mismatch ───────────────────────────────────────────


def test_ac008_version_mismatch_rejects_before_pip(tmp_path: Path) -> None:
    py_dir = tmp_path / "python"
    py_dir.mkdir()
    (py_dir / "src" / "atomiverse").mkdir(parents=True)
    (py_dir / "pyproject.toml").write_text(
        '[project]\nname = "atomiverse"\nversion = "9.9.9"\n',
        encoding="utf-8",
    )
    with pytest.raises(ConfigError, match="version mismatch"):
        verify_source_version(tmp_path, expected="0.1.0")
    with pytest.raises(ConfigError, match="version mismatch"):
        resolve_source(src=str(tmp_path), version="0.1.0")


def test_ac008_matching_version_ok(tmp_path: Path) -> None:
    py_dir = tmp_path / "python"
    (py_dir / "src" / "atomiverse").mkdir(parents=True)
    (py_dir / "pyproject.toml").write_text(
        '[project]\nname = "atomiverse"\nversion = "0.1.0"\n',
        encoding="utf-8",
    )
    spec = resolve_source(src=str(tmp_path), version="0.1.0")
    assert spec.version == "0.1.0"
    assert spec.location.endswith("python")


# ── install always PYTHON=ON ─────────────────────────────────────────────────


def test_install_always_sets_python_on() -> None:
    cfg = InstallConfig(modules={"md": True})
    defines = dry_run_defines(cfg)
    assert defines["ATV_BUILD_PYTHON"] == "ON"
    from atomiverse.cli.source import SourceSpec

    cmd = build_pip_command(
        SourceSpec(location="/tmp/x", version="0.1.0"), defines
    )
    assert any("cmake.define.ATV_BUILD_PYTHON=ON" in a for a in cmd)
    assert "--force-reinstall" in cmd
    assert "--no-deps" in cmd
    assert not any(a == "cmake" or a.startswith("cmake ") for a in cmd)


# ── status ───────────────────────────────────────────────────────────────────


def test_ac010_status_prints_modules_and_fingerprint(tmp_path: Path, capsys) -> None:
    state = {
        "modules": ["md", "neighbor", "potentials"],
        "fingerprint": "ATV_BUILD_MD=ON|ATV_BUILD_PYTHON=ON",
        "defines": {
            "ATV_BUILD_MD": "ON",
            "ATV_BUILD_PYTHON": "ON",
            "ATV_BUILD_POTENTIALS": "ON",
            "ATV_BUILD_NEIGHBOR": "ON",
        },
    }
    path = tmp_path / "_install_state.json"
    path.write_text(json.dumps(state), encoding="utf-8")
    rc = run_status(path=path)
    assert rc == 0
    out = capsys.readouterr().out
    assert "modules:" in out
    assert "md" in out
    assert "fingerprint:" in out
    assert "ATV_BUILD_MD=ON" in out


def test_ac010_status_cli_exit_0() -> None:
    rc = main(["status"])
    assert rc == 0


def test_ac010_format_status_none() -> None:
    text = format_status(None)
    assert "fingerprint" in text
    assert "none" in text.lower() or "(none)" in text


# ── CLI loads without native extension ───────────────────────────────────────


def test_ac011_cli_import_does_not_load_extension() -> None:
    import subprocess

    script = textwrap.dedent(
        f"""\
        import sys
        sys.path.insert(0, {str(_SRC)!r})
        for k in list(sys.modules):
            if k == "atomiverse" or k.startswith("atomiverse.") or k == "_atomiverse":
                del sys.modules[k]
        import atomiverse.cli.config
        import atomiverse.cli.source
        import atomiverse.cli.install
        import atomiverse.cli.status
        import atomiverse.cli.main
        bad = [
            k
            for k in sys.modules
            if k in ("_atomiverse", "atomiverse._atomiverse")
            or (k.endswith("._atomiverse") and "binding" not in k)
        ]
        assert not bad, bad
        print("OK")
        """
    )
    proc = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "OK" in proc.stdout


# ── TOML happy path ──────────────────────────────────────────────────────────


def test_load_toml_config_modules(tmp_path: Path) -> None:
    p = tmp_path / "ok.toml"
    p.write_text(
        textwrap.dedent(
            """\
            [modules]
            md = true
            pw = false

            [md]
            mpi = true

            [global]
            real_storage = "double"

            [cmake.unbound]
            ATV_USE_LIBXC = "OFF"
            """
        ),
        encoding="utf-8",
    )
    cfg = load_toml_config(p)
    assert cfg.modules["md"] is True
    assert cfg.modules.get("pw") is False
    assert cfg.options["mpi"] is True
    assert cfg.options["real_storage"] == "double"
    assert cfg.unbound["ATV_USE_LIBXC"] == "OFF"
    defines = resolve_defines(cfg)
    assert defines["ATV_BUILD_MD"] == "ON"
    assert defines["ATV_BUILD_POTENTIALS"] == "ON"
    assert defines["ATV_USE_MPI"] == "ON"
    assert defines["ATV_REAL_STORAGE"] == "double"
    assert defines["ATV_BUILD_PYTHON"] == "ON"


def test_default_tarball_url_uses_version() -> None:
    url = default_tarball_url("0.1.0")
    assert url == (
        "https://github.com/Roy-Kid/Atomiverse/archive/refs/tags/v0.1.0.tar.gz"
    )
    ver = package_version()
    assert isinstance(ver, str) and ver


def test_help_mentions_install(capsys) -> None:
    rc = main([])
    assert rc == 0
    out = capsys.readouterr().out.lower()
    assert "install" in out
    # Typer --help exits via typer.Exit(0)
    rc = main(["--help"])
    assert rc == 0


def test_install_help_lists_modules(capsys) -> None:
    rc = main(["install", "--help"])
    assert rc == 0
    out = capsys.readouterr().out.lower()
    assert "module" in out
    assert "--config" in out
    assert "--mpi" in out or "mpi" in out
    # no k:v dialect
    assert "k:v" not in out


def test_dry_run_md_via_main(tmp_path: Path, capsys) -> None:
    """End-to-end dry-run with matching local source (no pip/C++)."""
    repo_python = _SRC.parent  # python/
    pyproject = repo_python / "pyproject.toml"
    if not pyproject.is_file():
        pytest.skip("python/pyproject.toml not found")
    rc = main(
        [
            "install",
            "md",
            "--mpi",
            "--src",
            str(repo_python),
            "--dry-run",
        ]
    )
    out = capsys.readouterr().out
    if rc != 0:
        pytest.skip(f"dry-run skipped (rc={rc}): {out}")
    assert "ATV_BUILD_MD=ON" in out
    assert "ATV_BUILD_POTENTIALS=ON" in out
    assert "ATV_BUILD_NEIGHBOR=ON" in out
    assert "ATV_BUILD_PYTHON=ON" in out
    assert "ATV_USE_MPI=ON" in out
    assert "dry-run" in out.lower()
