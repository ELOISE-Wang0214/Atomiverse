"""Install configuration: binding table, TOML/inline parse, auto-enable DAG.

Maps user-facing module flags and friendly option keys onto the CMake
``ATV_*`` defines declared in ``cmake/AtvProductModules.cmake`` (phase 01).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Optional, Set

# ── Module flags (CLI / [modules] section) → product ATV_BUILD_* ──────────
# Names match the user-facing CLI surface. ``am1`` maps to SEMIEMPIRICAL;
# ``minimize`` maps to MINIMIZER. ``neighbor`` is optional-explicit (also
# pulled in by the auto-enable closure).

MODULE_TO_DEFINE: Dict[str, str] = {
    "md": "ATV_BUILD_MD",
    "minimize": "ATV_BUILD_MINIMIZER",
    "mc": "ATV_BUILD_MC",
    "evaluate": "ATV_BUILD_EVALUATE",
    "vib": "ATV_BUILD_VIB",
    "path": "ATV_BUILD_PATH",
    "pw": "ATV_BUILD_PW",
    "gpw": "ATV_BUILD_GPW",
    "am1": "ATV_BUILD_SEMIEMPIRICAL",
    "cuda": "ATV_BUILD_CUDA",
    "kspace": "ATV_BUILD_KSPACE",
    "potentials": "ATV_BUILD_POTENTIALS",
    "neighbor": "ATV_BUILD_NEIGHBOR",
}

# Product modules that are set ON/OFF explicitly on every resolve so monorepo
# CMake defaults (mostly ON) do not leak into an on-demand install.
PRODUCT_MODULES: tuple[str, ...] = tuple(MODULE_TO_DEFINE.keys())

# Auto-enable closure — same one-way DAG as cmake/AtvProductModules.cmake:
#   MD, MC, MINIMIZER → POTENTIALS + NEIGHBOR
#   PATH → MINIMIZER  (then inherits minimizer → potentials+neighbor)
#   KSPACE → POTENTIALS
AUTO_ENABLE: Dict[str, tuple[str, ...]] = {
    "md": ("potentials", "neighbor"),
    "mc": ("potentials", "neighbor"),
    "minimize": ("potentials", "neighbor"),
    "path": ("minimize",),
    "kspace": ("potentials",),
}

# Friendly option keys → ATV_* CMake cache names (01 names).
# Cross-module options live under [global] / --global only.
# Module flags themselves are also bound (cannot pass as -D / unbound).
BINDING_TABLE: Dict[str, str] = {
    # module enable flags
    **MODULE_TO_DEFINE,
    # common options
    "mpi": "ATV_USE_MPI",
    "real_storage": "ATV_REAL_STORAGE",
    "real_work": "ATV_REAL_WORK",
    "id_size": "ATV_ID_SIZE",
}

# Reverse: any ATV_* name that is bound (cannot appear under cmake.unbound / -D).
BOUND_DEFINES: Set[str] = set(BINDING_TABLE.values()) | {"ATV_BUILD_PYTHON"}

_BOOL_TRUE = frozenset({"yes", "true", "on", "1"})
_BOOL_FALSE = frozenset({"no", "false", "off", "0"})

HF_REJECT_MSG = (
    "HF/SCF Python bindings are not available yet. "
    "Remove --hf / modules.hf until SCF bindings ship."
)


class ConfigError(ValueError):
    """User-facing configuration error (fail-loud)."""


@dataclass
class InstallConfig:
    """Resolved install request before define mapping.

    @param modules: Module name → enabled bool (only explicitly requested /
        auto-enabled modules need be present; resolve fills the rest).
    @param options: Bound friendly-key → value (bool or str).
    @param unbound: Unbound CMake cache key → value (stringified).
    """

    modules: Dict[str, bool] = field(default_factory=dict)
    options: Dict[str, Any] = field(default_factory=dict)
    unbound: Dict[str, str] = field(default_factory=dict)


def parse_bool(raw: str) -> Optional[bool]:
    """Parse a bool synonym; return None if not a recognised bool token."""
    key = raw.strip().lower()
    if key in _BOOL_TRUE:
        return True
    if key in _BOOL_FALSE:
        return False
    return None


def parse_inline_kv(spec: str) -> Dict[str, Any]:
    """Parse ``k:v,k:v`` inline options.

    Items are comma-separated; the **first** ``:`` splits key from value.
    Bool synonyms ``yes/no/true/false/on/off/1/0`` become ``bool``; other
    values stay as strings.

    @param spec: Inline option string (may be empty).
    @return Mapping of friendly key → value.
    @raises ConfigError: on empty key or missing colon.
    """
    out: Dict[str, Any] = {}
    text = (spec or "").strip()
    if not text:
        return out
    for item in text.split(","):
        piece = item.strip()
        if not piece:
            continue
        if ":" not in piece:
            raise ConfigError(
                f"inline option {piece!r} missing key:value (first-colon split)"
            )
        key, _, value = piece.partition(":")
        key = key.strip()
        value = value.strip()
        if not key:
            raise ConfigError(f"inline option has empty key in {piece!r}")
        b = parse_bool(value)
        out[key] = b if b is not None else value
    return out


def _load_toml(path: Path) -> Mapping[str, Any]:
    try:
        import tomllib  # py>=3.11
    except ModuleNotFoundError:  # pragma: no cover - exercised on 3.9/3.10
        try:
            import tomli as tomllib  # type: ignore[no-redef]
        except ModuleNotFoundError as exc:
            raise ConfigError(
                "TOML config requires Python 3.11+ (tomllib) or the 'tomli' package"
            ) from exc
    with path.open("rb") as fh:
        return tomllib.load(fh)


def load_toml_config(path: Path | str) -> InstallConfig:
    """Load unified TOML config (``[modules]``, per-module, ``[global]``, unbound).

    @param path: Path to a ``.toml`` file.
    @return Parsed :class:`InstallConfig`.
    @raises ConfigError: on unknown sections/keys, hf enable, bound-as-unbound.
    """
    p = Path(path)
    if not p.is_file():
        raise ConfigError(f"config file not found: {p}")
    data = dict(_load_toml(p))
    cfg = InstallConfig()

    modules_tbl = data.pop("modules", None) or {}
    if not isinstance(modules_tbl, Mapping):
        raise ConfigError("[modules] must be a table")
    for name, enabled in modules_tbl.items():
        name = str(name)
        if name == "hf":
            if enabled:
                raise ConfigError(HF_REJECT_MSG)
            continue
        if name == "global":
            raise ConfigError("'global' is not a product module; use [global]")
        if name not in MODULE_TO_DEFINE:
            raise ConfigError(f"unknown module in [modules]: {name!r}")
        cfg.modules[name] = bool(enabled)

    global_tbl = data.pop("global", None) or {}
    if not isinstance(global_tbl, Mapping):
        raise ConfigError("[global] must be a table")
    for key, value in global_tbl.items():
        _set_bound_option(cfg, str(key), value, where="[global]")

    cmake_tbl = data.pop("cmake", None) or {}
    if not isinstance(cmake_tbl, Mapping):
        raise ConfigError("[cmake] must be a table")
    unbound_tbl = cmake_tbl.get("unbound", {}) if cmake_tbl else {}
    # allow only the unbound sub-table
    extra_cmake = set(cmake_tbl.keys()) - {"unbound"}
    if extra_cmake:
        raise ConfigError(
            f"[cmake] only supports [cmake.unbound], not {sorted(extra_cmake)}"
        )
    if not isinstance(unbound_tbl, Mapping):
        raise ConfigError("[cmake.unbound] must be a table")
    for key, value in unbound_tbl.items():
        _set_unbound(cfg, str(key), value)

    # remaining top-level tables are per-module option sections
    for section, body in data.items():
        if section in ("modules", "global", "cmake"):
            continue
        if section == "hf":
            raise ConfigError(HF_REJECT_MSG)
        if section not in MODULE_TO_DEFINE:
            raise ConfigError(f"unknown config section [{section}]")
        if not isinstance(body, Mapping):
            raise ConfigError(f"[{section}] must be a table")
        # presence of a per-module table does not alone enable the module;
        # [modules] does. But if the module is listed with options, require
        # it be enabled (or enable it implicitly when the table is non-empty
        # only if modules already said so). Spec: modules table enables.
        for key, value in body.items():
            _set_bound_option(cfg, str(key), value, where=f"[{section}]")

    return cfg


def _stringify_cmake_value(value: Any) -> str:
    if isinstance(value, bool):
        return "ON" if value else "OFF"
    return str(value)


def _set_bound_option(
    cfg: InstallConfig, key: str, value: Any, *, where: str
) -> None:
    if key not in BINDING_TABLE:
        raise ConfigError(
            f"unknown bound key {key!r} in {where}; "
            f"known keys: {sorted(BINDING_TABLE)}"
        )
    # module names used as options are not valid under [global]/module]
    if key in MODULE_TO_DEFINE:
        raise ConfigError(
            f"{key!r} is a module flag, not an option key (set it under [modules])"
        )
    if isinstance(value, str):
        b = parse_bool(value)
        cfg.options[key] = b if b is not None else value
    else:
        cfg.options[key] = value


def _set_unbound(cfg: InstallConfig, key: str, value: Any) -> None:
    # Reject friendly bound keys and their ATV_* targets.
    if key in BINDING_TABLE or key in BOUND_DEFINES:
        raise ConfigError(
            f"bound key {key!r} cannot be passed as cmake.unbound / -D; "
            "use the friendly module flag or [global] option instead"
        )
    # Also reject if someone passes the friendly name that maps to a define
    for friendly, define in BINDING_TABLE.items():
        if key == define or key == friendly:
            raise ConfigError(
                f"bound key {key!r} cannot be passed as cmake.unbound / -D"
            )
    cfg.unbound[key] = _stringify_cmake_value(value)


def apply_inline_modules(
    cfg: InstallConfig,
    module_flags: Mapping[str, Optional[str]],
    global_spec: Optional[str] = None,
) -> InstallConfig:
    """Merge inline ``--md [k:v]`` / ``--global k:v`` into *cfg*.

    Presence of a module flag enables that module (even with empty k:v).

    @param cfg: Config to mutate (returned for chaining).
    @param module_flags: Module name → inline k:v string or None (absent).
    @param global_spec: Optional ``--global`` inline k:v string.
    """
    for name, spec in module_flags.items():
        if spec is None:
            continue
        if name == "hf":
            raise ConfigError(HF_REJECT_MSG)
        if name not in MODULE_TO_DEFINE:
            raise ConfigError(f"unknown module flag --{name}")
        cfg.modules[name] = True
        for key, value in parse_inline_kv(spec).items():
            _set_bound_option(cfg, key, value, where=f"--{name}")
    if global_spec is not None:
        for key, value in parse_inline_kv(global_spec).items():
            _set_bound_option(cfg, key, value, where="--global")
    return cfg


def apply_auto_enable(modules: MutableMapping[str, bool]) -> Set[str]:
    """Close *modules* under the product-module DAG (01 semantics).

    @param modules: Module enable map (mutated in place).
    @return Set of module names that are enabled after closure.
    """
    # iterate to fixed point (path → minimize → potentials+neighbor)
    changed = True
    while changed:
        changed = False
        enabled = [m for m, on in modules.items() if on]
        for m in enabled:
            for dep in AUTO_ENABLE.get(m, ()):
                if not modules.get(dep, False):
                    modules[dep] = True
                    changed = True
    return {m for m, on in modules.items() if on}


def resolve_defines(cfg: InstallConfig) -> Dict[str, str]:
    """Map an :class:`InstallConfig` to CMake ``ATV_*`` define strings.

    Always forces ``ATV_BUILD_PYTHON=ON``. Every product module is set
    explicitly ON or OFF so monorepo defaults do not leak.

    @param cfg: Parsed install configuration.
    @return Ordered-stable dict of define name → ``ON``/``OFF``/value.
    @raises ConfigError: empty install (no modules enabled).
    """
    modules: Dict[str, bool] = {m: False for m in PRODUCT_MODULES}
    for name, on in cfg.modules.items():
        if name == "hf" and on:
            raise ConfigError(HF_REJECT_MSG)
        if name not in MODULE_TO_DEFINE:
            raise ConfigError(f"unknown module {name!r}")
        modules[name] = bool(on)

    enabled = apply_auto_enable(modules)
    if not enabled:
        raise ConfigError(
            "empty install: enable at least one module "
            "(e.g. --md, --pw, or [modules] in --config)"
        )

    defines: Dict[str, str] = {}
    for name in PRODUCT_MODULES:
        defines[MODULE_TO_DEFINE[name]] = "ON" if modules.get(name) else "OFF"

    defines["ATV_BUILD_PYTHON"] = "ON"

    for key, value in cfg.options.items():
        if key not in BINDING_TABLE:
            raise ConfigError(f"unknown bound option {key!r}")
        define = BINDING_TABLE[key]
        defines[define] = _stringify_cmake_value(value)

    for key, value in cfg.unbound.items():
        if key in BOUND_DEFINES or key in BINDING_TABLE:
            raise ConfigError(
                f"bound key {key!r} cannot be passed as cmake.unbound / -D"
            )
        defines[key] = value

    return defines


def enabled_modules_from_defines(defines: Mapping[str, str]) -> List[str]:
    """List product modules whose define is ON (stable sorted)."""
    on: List[str] = []
    for name, define in MODULE_TO_DEFINE.items():
        if defines.get(define, "OFF").upper() in ("ON", "TRUE", "1"):
            on.append(name)
    return on


def define_fingerprint(defines: Mapping[str, str]) -> str:
    """Stable fingerprint of a define map (sorted ``k=v`` joined)."""
    parts = [f"{k}={defines[k]}" for k in sorted(defines)]
    return "|".join(parts)


def any_module_flag_present(module_flags: Mapping[str, Optional[str]]) -> bool:
    """True if any module flag was supplied on the CLI (value not None)."""
    return any(v is not None for v in module_flags.values())


def validate_mode(
    *,
    config_path: Optional[str],
    module_flags: Mapping[str, Optional[str]],
    global_spec: Optional[str],
    hf: bool,
) -> None:
    """Enforce XOR of ``--config`` vs inline flags; reject ``--hf`` and empty.

    @raises ConfigError: on mutual exclusion, empty install, or --hf.
    """
    if hf:
        raise ConfigError(HF_REJECT_MSG)

    inline = any_module_flag_present(module_flags) or global_spec is not None
    if config_path and inline:
        raise ConfigError(
            "mutual exclusion: use either --config PATH.toml or inline "
            "module flags (--md/--pw/…/--global), not both"
        )
    if not config_path and not any_module_flag_present(module_flags):
        # --global alone is not an install
        raise ConfigError(
            "empty install: pass --config PATH.toml or at least one module "
            "flag (e.g. --md, --pw)"
        )


def merge_unbound_cli(cfg: InstallConfig, pairs: Iterable[str]) -> None:
    """Apply ``-D KEY=VALUE`` / ``--define KEY=VALUE`` CLI unbound pairs."""
    for raw in pairs:
        if "=" not in raw:
            raise ConfigError(f"-D expects KEY=VALUE, got {raw!r}")
        key, _, value = raw.partition("=")
        key = key.strip()
        # strip optional leading -D
        if key.startswith("-D"):
            key = key[2:]
        _set_unbound(cfg, key, value.strip())
