"""Device handlers for molpy integration.

Registers Atomiverse compute backends into molpy's device registry,
enabling ``frame.to("cuda")`` / ``frame.to("cpu")`` via the molpy protocol.

This module is NOT auto-imported. Users opt in explicitly::

    import atomiverse
    atomiverse.register_molpy()   # raises if molpy is missing or incompatible

Or::

    from atomiverse._device import register_devices
    register_devices()
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    pass  # mp.Frame type used in signatures below


class CpuDeviceHandler:
    """Transfers mp.Frame data to/from Atomiverse cpu::Frame."""

    def from_cpu(self, frame: object) -> object:
        """mp.Frame(cpu) → internal CpuFrame, run, write back."""
        # TODO: implement when molpy Frame protocol is ready
        raise NotImplementedError(
            "CpuDeviceHandler.from_cpu() is not yet implemented. "
            "Use av.Frame for direct Atomiverse access."
        )

    def to_cpu(self, frame: object) -> object:
        """Internal representation → mp.Frame(cpu)."""
        raise NotImplementedError(
            "CpuDeviceHandler.to_cpu() is not yet implemented. "
            "Use av.Frame for direct Atomiverse access."
        )


class CudaDeviceHandler:
    """Transfers mp.Frame data to/from Atomiverse cuda::Frame via DLPack."""

    def from_cpu(self, frame: object) -> object:
        """mp.Frame(cpu) → mp.Frame(cuda) via H2D transfer."""
        raise NotImplementedError(
            "CudaDeviceHandler.from_cpu() is not yet implemented. "
            "Requires ATV_BUILD_CUDA=ON."
        )

    def to_cpu(self, frame: object) -> object:
        """mp.Frame(cuda) → mp.Frame(cpu) via D2H transfer."""
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
