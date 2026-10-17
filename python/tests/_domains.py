"""Which product modules the installed extension carries.

A native wheel holds only the modules it was configured with, so a test that
needs a domain says so and is skipped, naming the flag, on an extension built
without it (or on a pure-Python install).
"""

from __future__ import annotations

import pytest

from atomiverse._binding import _cpp, extension_available

# A symbol each domain's binding file always registers when it is built.
_PROBES: dict[str, str] = {
    "ATV_BUILD_POTENTIALS": "PairLJCut",
    "ATV_BUILD_MD": "MD",
    "ATV_BUILD_MINIMIZER": "Minimize",
    "ATV_BUILD_MC": "MC",
    "ATV_BUILD_VIB": "Hessian",
    "ATV_BUILD_SEMIEMPIRICAL": "AM1",
    "ATV_BUILD_HF": "RHF",
    "ATV_BUILD_PW": "PW",
    "ATV_BUILD_KSPACE": "KSpacePME",
}


def built(*flags: str) -> bool:
    """True when the extension is installed and carries every ``flags`` domain."""
    return extension_available() and all(hasattr(_cpp, _PROBES[f]) for f in flags)


def needs(*flags: str) -> pytest.MarkDecorator:
    """Skip mark for a test, class or parameter that needs ``flags``."""
    return pytest.mark.skipif(not built(*flags), reason="needs " + " and ".join(flags))


def skip_module_unless(*flags: str) -> None:
    """Skip the calling test file unless the extension carries ``flags``.

    Call it before the file's compute imports, which fail without the domain.
    """
    if not built(*flags):
        pytest.skip("needs " + " and ".join(flags), allow_module_level=True)
