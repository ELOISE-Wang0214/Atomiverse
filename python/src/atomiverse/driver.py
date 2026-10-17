"""Driver entry points — one-shot evaluations that need no loop of their own.

Canonical deep path: ``atomiverse.driver.single_point``; top-level
``atomiverse.single_point`` is the same function. Its C++ counterpart is
``cpu::single_point`` (``cpu/evaluate/single_point.h``), the helper that
runs one pass of the ``Evaluate`` driver preset.
"""

from atomiverse._binding import _cpp, has_domain, missing_domain_message, require_cpp


def single_point(frame: object, evaluator: object, device: str = "cpu") -> None:
    """Compute energy (and any evaluator-specific properties) for one fixed
    geometry — a *single point*, as opposed to a trajectory or an
    optimization.

    Energy lands in the frame's typed ``meta`` map [eV]; charge-capable
    evaluators (AM1) write per-atom charges [e] back into the same
    ``molpy.Frame``. The evaluator is borrowed — inspect
    ``evaluator.last_result()`` afterwards for solver-internal quantities.

    Two dispatch paths, one call shape. A component-style evaluator such as
    ``AM1Evaluator`` is run by the C++ ``cpu::single_point`` helper, which
    drives one pass of the ``Evaluate`` driver over it. A Gaussian SCF method
    (``RHF`` / ``UHF`` / ``RKS`` / ``UKS``) is already a driver in its own
    right — it owns the self-consistent-field loop — so it is dispatched
    straight to ``evaluator.run(frame)``. Nothing is returned either way;
    results are read off the frame or off the evaluator.

    Parameters
    ----------
    frame : molpy.Frame
        System state (``atoms`` block with ``x, y, z, atomic_number``; the
        ``atomic_number`` column [uint] carries Z). Results land as typed
        metadata read back as plain values, for example
        ``frame.meta["energy"]`` is a ``float`` [eV].
    evaluator : atomiverse.AM1Evaluator or atomiverse.RHF/UHF/RKS/UKS
        The evaluator component. ``AM1Evaluator`` runs through the C++
        ``Evaluate`` pass; a Gaussian SCF method owns its own SCF loop, so
        it is driven by ``evaluator.run(frame)`` instead (read the Hartree
        energy / loop state off the method afterwards with ``energy()`` /
        ``converged()`` / ``iterations()``).

        Forces differ between the two paths. The ``Evaluate`` pass requests
        gradients, so ``AM1Evaluator`` publishes fx/fy/fz [eV/Å] here as
        well; ``AM1Evaluator.forces(frame)`` is the entry that also *returns*
        them as an (N, 3) array. A Gaussian SCF ``run`` computes no forces —
        call that method's own ``forces(frame)`` for those.
    device : str, optional
        Execution backend, PyTorch-style. v1 accepts only ``"cpu"``; any
        other value raises ``ValueError`` before any compute.

    Raises
    ------
    ValueError
        If *device* is anything but ``"cpu"``, or (SCF path) if the method's
        spin multiplicity cannot split this frame's electron count.
    TypeError
        If *evaluator* is not a Gaussian SCF method and the extension was
        built without the ``am1`` module, which binds the ``Evaluate`` pass.
    RuntimeError
        If the frame is missing a required column or carries a wrong dtype.
    """
    if device != "cpu":
        raise ValueError(
            f"device={device!r}: CUDA bindings are not yet available; "
            "only device='cpu' is supported"
        )
    require_cpp()
    # A Gaussian SCF method is its own driver loop, not an Evaluate
    # component. Imported here, and only on an `hf` build, so this module
    # imports on any module selection.
    if has_domain("hf"):
        from atomiverse import scf

        if isinstance(evaluator, (scf.RHF, scf.UHF, scf.RKS, scf.UKS)):
            evaluator.run(frame)
            return
    if not has_domain("am1"):
        raise TypeError(
            f"single_point cannot run {type(evaluator).__name__}: "
            + missing_domain_message("the Evaluate pass of single_point", "am1")
        )
    _cpp.single_point(frame, evaluator)


__all__ = ["single_point"]
