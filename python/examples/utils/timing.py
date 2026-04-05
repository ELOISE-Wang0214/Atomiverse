"""Timing decorator and context manager for Atomiverse examples."""

from __future__ import annotations

import functools
import time
from typing import Callable


def timed(fn: Callable | None = None, *, label: str | None = None):
    """Decorator that prints wall-clock time after function returns.

    Usage::

        @timed
        def heavy():
            ...

        @timed(label="PySCF RHF")
        def run_pyscf():
            ...
    """
    def decorator(func: Callable) -> Callable:
        name = label or func.__name__

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            t0 = time.perf_counter()
            result = func(*args, **kwargs)
            elapsed = time.perf_counter() - t0
            wrapper.elapsed = elapsed
            wrapper.elapsed_ms = elapsed * 1e3
            return result

        wrapper.elapsed = 0.0
        wrapper.elapsed_ms = 0.0
        wrapper.label = name
        return wrapper

    if fn is not None:
        return decorator(fn)
    return decorator
