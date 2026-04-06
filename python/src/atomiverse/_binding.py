"""Internal shim: imports the C++ extension module as ``_cpp``.

Two valid import paths:
  - ``atomiverse._atomiverse`` — installed wheel
  - ``_atomiverse``           — development (bare .so on PYTHONPATH)

If BOTH fail, the C++ extension was not built. Raise immediately.
"""

try:
    from atomiverse import _atomiverse as _cpp  # installed wheel
except ImportError:
    try:
        import _atomiverse as _cpp  # type: ignore[import-untyped]
    except ImportError:
        raise ImportError(
            "Cannot import the _atomiverse C++ extension.\n"
            "Build it with:\n"
            "  cmake -B build -DATV_BUILD_PYTHON=ON && cmake --build build\n"
            "Then either install the wheel or add the build directory to PYTHONPATH."
        ) from None

__all__ = ["_cpp"]
