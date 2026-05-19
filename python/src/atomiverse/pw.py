"""Plane-wave DFT for periodic systems.

Provides ``pw_scf`` for total energy calculation using
plane-wave basis, norm-conserving pseudopotentials, and LDA.

Example::

    import atomiverse as av

    # Build a structure
    frame = av.Frame(2)
    frame.atoms.type[:] = [14, 14]  # Si
    frame.atoms.mass[:] = [28.0855, 28.0855]
    frame.atoms.x[:] = [0.0, 1.3575]
    frame.atoms.y[:] = [0.0, 1.3575]
    frame.atoms.z[:] = [0.0, 1.3575]
    # ... set frame.box for FCC cell

    result = av.pw_scf(frame, "Si.upf", ecut_ry=30.0, k_mesh=4)
    print(result.energy)  # Ry
"""

from __future__ import annotations

from atomiverse._binding import _cpp


class PWResult:
    """Result of a plane-wave DFT calculation.

    Attributes
    ----------
    energy : float
        Total energy in Rydberg.
    converged : bool
        Whether the SCF loop converged within max_iter.
    iterations : int
        Number of SCF iterations taken.
    """

    def __init__(self, cpp_result: object) -> None:
        self._cpp = cpp_result

    @property
    def energy(self) -> float:
        return float(self._cpp.energy)

    @property
    def converged(self) -> bool:
        return bool(self._cpp.converged)

    @property
    def iterations(self) -> int:
        return int(self._cpp.iterations)

    def __repr__(self) -> str:
        conv = "converged" if self.converged else "NOT converged"
        return (
            f"PWResult({conv}, E={self.energy:.6f} Ry, "
            f"iters={self.iterations})"
        )


def pw_scf(
    frame: object,
    upf_path: str = "",
    ecut_ry: float = 30.0,
    k_mesh: int = 4,
    tol: float = 1e-6,
    max_iter: int = 50,
) -> PWResult:
    """Run plane-wave DFT total energy for a periodic system.

    Uses a Monkhorst-Pack k-point mesh and LDA exchange-correlation.
    The Frame provides atom positions, types (Z), and simulation cell.

    Parameters
    ----------
    frame : Frame
        Input structure with positions, atom types (Z), and cell.
    upf_path : str
        Path to UPF pseudopotential file. Empty = toy local potential.
    ecut_ry : float
        Plane-wave energy cutoff in Rydberg. Default 30.0.
    k_mesh : int
        Monkhorst-Pack mesh size per dimension (e.g. 4 = 4x4x4).
    tol : float
        SCF convergence tolerance on density change.
    max_iter : int
        Maximum SCF iterations.

    Returns
    -------
    PWResult
        Result with ``.energy`` (Ry), ``.converged``, ``.iterations``.
    """
    return PWResult(
        _cpp.pw_scf(frame, upf_path, ecut_ry, k_mesh, tol, max_iter)
    )
