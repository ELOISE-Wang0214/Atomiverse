"""Device handlers for molpy integration.

Registers Atomiverse compute backends into molpy's device registry,
enabling ``frame.to("cuda")`` / ``frame.to("cpu")`` via the molpy protocol.

This module is NOT auto-imported. Users opt in explicitly::

    import atomiverse
    atomiverse.register_molpy()   # raises if molpy is missing or incompatible

Or::

    from atomiverse._device import register_devices
    register_devices()

Note: Atomiverse drivers operate directly on a ``molrs.Frame``. These
handlers exist only for the optional molpy device-registry protocol.
"""

from __future__ import annotations


class CpuDeviceHandler:
    """Placeholder CPU device handler for the molpy device registry.

    Atomiverse drivers consume a ``molrs.Frame`` directly via
    ``MD.run(frame, n_steps)`` — no device handler is needed for the
    normal path. This handler is a stub for the optional molpy
    ``frame.to("cpu")`` protocol.
    """

    def from_cpu(self, frame: object) -> object:
        """molpy Frame(cpu) → Atomiverse CPU execution path."""
        raise NotImplementedError(
            "CpuDeviceHandler.from_cpu() is not yet implemented. "
            "Pass a molrs.Frame to the Atomiverse driver directly."
        )

    def to_cpu(self, frame: object) -> object:
        """Atomiverse CPU representation → molpy Frame(cpu)."""
        raise NotImplementedError(
            "CpuDeviceHandler.to_cpu() is not yet implemented. "
            "Pass a molrs.Frame to the Atomiverse driver directly."
        )


class CudaDeviceHandler:
    """Placeholder CUDA device handler for the molpy device registry."""

    def from_cpu(self, frame: object) -> object:
        """molpy Frame(cpu) → Frame(cuda) via H2D transfer."""
        raise NotImplementedError(
            "CudaDeviceHandler.from_cpu() is not yet implemented. "
            "Requires ATV_BUILD_CUDA=ON."
        )

    def to_cpu(self, frame: object) -> object:
        """molpy Frame(cuda) → Frame(cpu) via D2H transfer."""
        raise NotImplementedError(
            "CudaDeviceHandler.to_cpu() is not yet implemented. "
            "Requires ATV_BUILD_CUDA=ON."
        )


def register_devices() -> None:
    """Register Atomiverse device handlers into molpy's device registry.

    Raises
    ------
    ImportError
        If molpy is not installed.
    AttributeError
        If molpy.Frame does not support the device registry protocol
        (``register_device`` classmethod).
    """
    try:
        import molpy as mp
    except ImportError:
        raise ImportError(
            "molpy is required for device registration.\n"
            "Install it with: pip install molpy"
        ) from None

    if not hasattr(mp.Frame, "register_device"):
        raise AttributeError(
            "molpy.Frame does not have a 'register_device' classmethod.\n"
            "The device registry protocol is not yet implemented in your "
            "version of molpy. Update molpy or implement the protocol:\n"
            "  Frame._device_registry: ClassVar[dict]\n"
            "  Frame.register_device(name, handler)\n"
            "  Frame.to(device) -> Frame"
        )

    mp.Frame.register_device("cpu", CpuDeviceHandler())
    mp.Frame.register_device("cuda", CudaDeviceHandler())
