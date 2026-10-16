"""Internal shim: soft-imports the C++ extension module as ``_cpp``.

Two valid import paths when the extension is present:
  - ``atomiverse._atomiverse`` — installed wheel (relative import)
  - ``_atomiverse``           — development (bare .so on PYTHONPATH)

When the extension is absent the package still imports:
``extension_available()`` is False and attribute access on the ``_cpp``
sentinel (or an explicit ``require_cpp()`` call) raises an install-oriented
``ImportError``. An extension that is present but fails to load (undefined
symbol, missing shared library) is a different failure: the package still
imports — so the CLI can rebuild it — but every compute access re-raises the
loader's own error rather than claiming nothing is installed.

Important: never use ``from atomiverse import _atomiverse`` here — that
re-enters package ``__init__`` while ``_binding`` is still loading and
creates a circular import.
"""

from __future__ import annotations

import importlib
from types import ModuleType
from typing import Any, Optional

__all__ = [
    "_cpp",
    "extension_available",
    "extension_load_error",
    "require_cpp",
    "INSTALL_HINT",
]

INSTALL_HINT = (
    "The Atomiverse C++ extension (_atomiverse) is not installed.\n"
    "Build it from a source checkout (git clone --recursive):\n"
    "  pip install 'atomiverse[native]'\n"
    "  python -m atomiverse config md --src /path/to/Atomiverse\n"
    "  python -m atomiverse compile\n"
    "  python -m atomiverse install"
)


class _MissingExtension:
    """Sentinel standing in for ``_atomiverse`` when the extension is absent.

    Any attribute access raises ``ImportError`` with the install hint so
    compute paths fail loud instead of with an opaque ``AttributeError``.
    """

    __slots__ = ()

    def __getattr__(self, name: str) -> Any:
        require_cpp()

    def __repr__(self) -> str:
        return "<atomiverse._binding: C++ extension not installed>"


_EXTENSION_NAMES = ("atomiverse._atomiverse", "_atomiverse")


def _try_import_extension() -> tuple[Optional[ModuleType], Optional[ImportError]]:
    """Return ``(module, None)``, ``(None, None)`` if absent, or
    ``(None, error)`` if an extension was found but failed to load."""
    # import_module, not `from . import`: inside the half-initialized package
    # the latter reports an absent submodule as a plain ImportError ("cannot
    # import name"), which would read as a broken extension.
    try:
        mod = importlib.import_module("atomiverse._atomiverse")

        return mod, None
    except ModuleNotFoundError as exc:
        if exc.name not in _EXTENSION_NAMES:
            return None, exc
    except ImportError as exc:
        return None, exc
    try:
        import _atomiverse as mod  # type: ignore[import-untyped]

        return mod, None
    except ModuleNotFoundError as exc:
        if exc.name not in _EXTENSION_NAMES:
            return None, exc
        return None, None
    except ImportError as exc:
        return None, exc


_ext, _load_error = _try_import_extension()
_cpp: Any = _ext if _ext is not None else _MissingExtension()


def extension_available() -> bool:
    """Return True when the ``_atomiverse`` C++ extension is importable."""
    return _ext is not None


def extension_load_error() -> Optional[ImportError]:
    """Return the loader error of an extension that is present but broken.

    @return The ``ImportError`` raised while loading ``_atomiverse``, or None
        when it loaded or is simply not installed.
    """
    return _load_error


def require_cpp() -> ModuleType:
    """Return the C++ extension module, or raise an install-oriented ImportError.

    @return The loaded ``_atomiverse`` module.
    @raises ImportError: when the extension has not been built/installed, or
        (chained from the loader's error) when it is present but fails to load.
    """
    if _ext is not None:
        return _ext
    if _load_error is not None:
        raise ImportError(
            f"The Atomiverse C++ extension failed to load: {_load_error}"
        ) from _load_error
    raise ImportError(INSTALL_HINT)
