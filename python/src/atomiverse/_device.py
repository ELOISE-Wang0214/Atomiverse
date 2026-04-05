"""Device handlers for molpy integration.

Registers Atomiverse compute backends into molpy's device registry,
enabling `frame.to("cuda")` via the standard molpy protocol.
"""

import numpy as np


class CpuDeviceHandler:
    """Handles CPU-resident frames via Atomiverse cpu::Frame."""

    def run(self, frame, sim, n_steps):
        """Run simulation on a CPU frame using Atomiverse backend."""
        from atomiverse import run_md
        return run_md(frame, sim, n_steps)


def register_devices():
    """Register Atomiverse device handlers into molpy.

    Call this once at import time. Safe to call multiple times.
    """
    try:
        import molpy as mp
        if hasattr(mp, "register_device"):
            mp.register_device("atomiverse_cpu", CpuDeviceHandler())
    except ImportError:
        pass
