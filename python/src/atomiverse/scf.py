"""Unified SCF solver — accepts ``mp.Frame``, converts internally.

Replaces the old per-method classes (RHF, UHF, UKS) with a single
``SCF`` class backed by the unified C++ ``scf_solve`` engine.
"""

from __future__ import annotations

import numpy as np
import molpy as mp

from atomiverse._binding import _cpp
from atomiverse._convert import mp_to_cpu


# ═══════════════════════════════════════════════════════════════════════
# Integrals
# ═══════════════════════════════════════════════════════════════════════

class Integrals:
    """Molecular integrals (S, T, V, Hcore, ERI) with property access.

    Users do not construct this directly — it is available as
    ``solver.integrals`` after creating an :class:`SCF` object.
    """

    def __init__(self, cpp_ip: object) -> None:
        self._cpp = cpp_ip

    @property
    def S(self) -> np.ndarray:
        return np.asarray(self._cpp.get_S())

    @property
    def T(self) -> np.ndarray:
        return np.asarray(self._cpp.get_T())

    @property
    def V(self) -> np.ndarray:
        return np.asarray(self._cpp.get_V())

    @property
    def Hcore(self) -> np.ndarray:
        return np.asarray(self._cpp.get_Hcore())

    @property
    def ERI(self) -> np.ndarray:
        return np.asarray(self._cpp.get_ERI())

    @property
    def V_nuc(self) -> float:
        return float(self._cpp.V_nuc)

    @property
    def n_basis(self) -> int:
        return int(self._cpp.n_basis)

    @property
    def n_electrons(self) -> int:
        return int(self._cpp.n_electrons)

    @property
    def n_occ(self) -> int:
        return int(self._cpp.n_occ)

    def __repr__(self) -> str:
        return (
            f"Integrals(n_basis={self.n_basis}, "
            f"n_electrons={self.n_electrons})"
        )


# ═══════════════════════════════════════════════════════════════════════
# Result classes
# ═══════════════════════════════════════════════════════════════════════

class SCFResult:
    """Result of a unified SCF calculation (RHF/UHF/UKS)."""

    def __init__(self, cpp_result: object) -> None:
        self._cpp = cpp_result

    @property
    def converged(self) -> bool:
        return bool(self._cpp.converged)

    @property
    def iterations(self) -> int:
        return int(self._cpp.iterations)

    @property
    def energy(self) -> float:
        """Total energy (electronic + nuclear) [Hartree]."""
        return float(self._cpp.energy)

    @property
    def residual(self) -> float:
        return float(self._cpp.residual)

    @property
    def one_electron_energy(self) -> float:
        return float(self._cpp.one_electron_energy)

    @property
    def coulomb_energy(self) -> float:
        return float(self._cpp.coulomb_energy)

    @property
    def exchange_energy(self) -> float:
        return float(self._cpp.exchange_energy)

    @property
    def xc_energy(self) -> float:
        return float(self._cpp.xc_energy)

    @property
    def nuclear_repulsion(self) -> float:
        return float(self._cpp.nuclear_repulsion)

    @property
    def n_alpha(self) -> int:
        return int(self._cpp.n_alpha)

    @property
    def n_beta(self) -> int:
        return int(self._cpp.n_beta)

    @property
    def spin_contamination(self) -> float:
        return float(self._cpp.spin_contamination)

    @property
    def eigenvalues_alpha(self) -> np.ndarray:
        return np.asarray(self._cpp.get_eigenvalues_alpha())

    @property
    def eigenvalues_beta(self) -> np.ndarray:
        return np.asarray(self._cpp.get_eigenvalues_beta())

    @property
    def eigenvalues(self) -> np.ndarray:
        """Orbital energies [Hartree]. Alpha channel (same as eigenvalues_alpha)."""
        return self.eigenvalues_alpha

    @property
    def density_matrix_alpha(self) -> np.ndarray:
        return np.asarray(self._cpp.get_P_alpha())

    @property
    def density_matrix_beta(self) -> np.ndarray:
        return np.asarray(self._cpp.get_P_beta())

    @property
    def density_matrix(self) -> np.ndarray:
        """Total density matrix P_alpha + P_beta."""
        return self.density_matrix_alpha + self.density_matrix_beta

    @property
    def coefficients_alpha(self) -> np.ndarray:
        return np.asarray(self._cpp.get_C_alpha())

    @property
    def coefficients_beta(self) -> np.ndarray:
        return np.asarray(self._cpp.get_C_beta())

    @property
    def trace(self) -> list[dict]:
        return list(self._cpp.trace)

    @property
    def integrated_n_alpha(self) -> float:
        return float(self._cpp.integrated_n_alpha)

    @property
    def integrated_n_beta(self) -> float:
        return float(self._cpp.integrated_n_beta)

    def __repr__(self) -> str:
        conv = "converged" if self.converged else "NOT converged"
        return (
            f"SCFResult({conv}, E={self.energy:.8f} Ha, "
            f"iters={self.iterations})"
        )


class GradientResult:
    """Result of a gradient computation."""

    def __init__(self, cpp_result: object) -> None:
        self._cpp = cpp_result

    @property
    def forces(self) -> np.ndarray:
        """Forces [n_atoms, 3] in eV/Angstrom."""
        return np.asarray(self._cpp.get_forces())

    @property
    def energy(self) -> float:
        """Total energy [eV]."""
        return float(self._cpp.energy)

    def __repr__(self) -> str:
        f = self.forces
        return f"GradientResult(E={self.energy:.6f} eV, n_atoms={f.shape[0]})"


# ═══════════════════════════════════════════════════════════════════════
# XC functional
# ═══════════════════════════════════════════════════════════════════════

_XC_MAP: dict[str, object] = {
    "LSDA": _cpp.XCType.LSDA,
    "lsda": _cpp.XCType.LSDA,
    "PBE": _cpp.XCType.PBE,
    "pbe": _cpp.XCType.PBE,
}

# Method string to (SpinType, MeanFieldType) mapping
_METHOD_MAP: dict[str, tuple[str, str]] = {
    "rhf": ("Restricted", "HF"),
    "uhf": ("Unrestricted", "HF"),
    "rks": ("Restricted", "KS"),
    "uks": ("Unrestricted", "KS"),
}


# ═══════════════════════════════════════════════════════════════════════
# Unified SCF solver
# ═══════════════════════════════════════════════════════════════════════

def _build_integrals(frame: mp.Frame, basis: str, charge: int) -> object:
    """Build integrals from mp.Frame."""
    if basis != "sto-3g":
        raise ValueError(
            f"Only 'sto-3g' basis is currently supported, got {basis!r}"
        )
    cf = mp_to_cpu(frame)
    return _cpp.build_integrals_sto3g(cf, charge)


class SCF:
    """Unified SCF solver — RHF, UHF, RKS, UKS.

    Parameters
    ----------
    frame : mp.Frame
        Molecular geometry. ``atoms.type`` must be atomic numbers (1=H, ..., 8=O).
    basis : str
        Basis set name (currently only ``"sto-3g"``).
    method : str
        SCF method: ``"rhf"``, ``"uhf"``, ``"rks"``, ``"uks"``.
    xc : str or None
        XC functional (``"LSDA"``, ``"PBE"``). Required for KS methods.
    charge : int
        Net molecular charge.
    multiplicity : int
        Spin multiplicity 2S+1 (1=singlet, 2=doublet, 3=triplet).

    Example::

        import molpy as mp
        import atomiverse as av

        frame = mp.Frame({"atoms": mp.Block({
            "x": [0.0, 0.96], "y": [0, 0.59], "z": [0, 0],
            "type": [8, 1], "mass": [15.999, 1.008],
        })})
        frame.box = mp.Box(np.eye(3) * 100.0)

        scf = av.SCF(frame, method="rhf")
        result = scf.solve()
        print(result.energy)
    """

    def __init__(
        self,
        frame: mp.Frame,
        basis: str = "sto-3g",
        method: str = "rhf",
        xc: str | None = None,
        charge: int = 0,
        multiplicity: int = 1,
        **kwargs: float,
    ) -> None:
        method_key = method.lower()
        if method_key not in _METHOD_MAP:
            raise ValueError(
                f"Unknown method {method!r}. "
                f"Supported: {', '.join(_METHOD_MAP)}"
            )

        spin_str, mf_str = _METHOD_MAP[method_key]
        is_ks = mf_str == "KS"

        if is_ks and xc is None:
            xc = "LSDA"
        if is_ks:
            xc_enum = _XC_MAP.get(xc)  # type: ignore[arg-type]
            if xc_enum is None:
                valid = ", ".join(f'"{k}"' for k in _XC_MAP if k == k.upper())
                raise ValueError(
                    f"Unknown XC functional {xc!r}. Supported: {valid}"
                )

        self._frame = frame
        self._ip_cpp = _build_integrals(frame, basis, charge)
        self.integrals = Integrals(self._ip_cpp)
        self._charge = charge

        n_el = self.integrals.n_electrons
        spin = multiplicity - 1
        if (n_el + spin) % 2 != 0:
            raise ValueError(
                f"Incompatible n_electrons={n_el} and "
                f"multiplicity={multiplicity}"
            )
        n_alpha = (n_el + spin) // 2
        n_beta = (n_el - spin) // 2

        # Build SCFConfig
        config = _cpp.SCFConfig()
        config.spin = getattr(_cpp.SpinType, spin_str)
        config.method = getattr(_cpp.MeanFieldType, mf_str)
        config.n_alpha = n_alpha
        config.n_beta = n_beta
        if is_ks:
            config.xc_type = xc_enum  # type: ignore[possibly-undefined]
            config.n_radial = int(kwargs.get("n_radial", 50))
            config.n_angular = int(kwargs.get("n_angular", 110))
        self._config = config

    def solve(self, tol: float = 1e-8, max_iter: int = 100) -> SCFResult:
        """Run the SCF procedure."""
        self._config.tol = tol
        self._config.max_iter = max_iter
        result = _cpp.scf_solve(self._ip_cpp, self._config)
        self._result = result
        return SCFResult(result)

    def gradient(self) -> GradientResult:
        """Compute the analytic nuclear gradient.

        Must call ``solve()`` first.
        """
        if not hasattr(self, "_result"):
            raise RuntimeError("Call solve() before gradient()")
        cf = mp_to_cpu(self._frame)
        grad = _cpp.scf_gradient(cf, self._ip_cpp, self._result, self._config)
        return GradientResult(grad)

    def __repr__(self) -> str:
        ip = self.integrals
        return (
            f"SCF(n_basis={ip.n_basis}, n_electrons={ip.n_electrons})"
        )


# ===================================================================
# Convenience wrappers (backward compatible API)
# ===================================================================


class RHF(SCF):
    """Restricted Hartree-Fock solver.

    Convenience wrapper around :class:`SCF` with ``method="rhf"``.

    Parameters
    ----------
    frame : mp.Frame
        Molecular geometry. ``atoms.type`` must be atomic numbers.
    basis : str
        Basis set name (currently only ``"sto-3g"``).
    charge : int
        Net molecular charge.

    Example::

        rhf = av.RHF(frame, basis="sto-3g")
        result = rhf.solve(tol=1e-10)
        print(result.energy)
    """

    def __init__(
        self,
        frame: mp.Frame,
        basis: str = "sto-3g",
        charge: int = 0,
    ) -> None:
        super().__init__(frame, basis=basis, method="rhf", charge=charge)

    def __repr__(self) -> str:
        ip = self.integrals
        return f"RHF(n_basis={ip.n_basis}, n_electrons={ip.n_electrons})"


class UHF(SCF):
    """Unrestricted Hartree-Fock solver.

    Convenience wrapper around :class:`SCF` with ``method="uhf"``.

    Parameters
    ----------
    frame : mp.Frame
        Molecular geometry.
    basis : str
        Basis set name.
    charge : int
        Net molecular charge.
    multiplicity : int
        Spin multiplicity 2S+1 (1=singlet, 2=doublet, 3=triplet).

    Example::

        uhf = av.UHF(frame, multiplicity=2)
        result = uhf.solve()
    """

    def __init__(
        self,
        frame: mp.Frame,
        basis: str = "sto-3g",
        charge: int = 0,
        multiplicity: int = 1,
    ) -> None:
        super().__init__(
            frame, basis=basis, method="uhf",
            charge=charge, multiplicity=multiplicity,
        )

    def __repr__(self) -> str:
        ip = self.integrals
        return (
            f"UHF(n_alpha={self._config.n_alpha}, "
            f"n_beta={self._config.n_beta}, n_basis={ip.n_basis})"
        )


class UKS(SCF):
    """Unrestricted Kohn-Sham (DFT) solver.

    Convenience wrapper around :class:`SCF` with ``method="uks"``.

    Parameters
    ----------
    frame : mp.Frame
        Molecular geometry.
    basis : str
        Basis set name.
    xc : str
        Exchange-correlation functional (``"LSDA"`` or ``"PBE"``).
    charge : int
        Net molecular charge.
    multiplicity : int
        Spin multiplicity 2S+1.
    n_radial : int
        Number of radial grid points.
    n_angular : int
        Number of angular grid points.

    Example::

        dft = av.UKS(frame, xc="PBE", multiplicity=1)
        result = dft.solve()
        print(result.xc_energy)
    """

    def __init__(
        self,
        frame: mp.Frame,
        basis: str = "sto-3g",
        xc: str = "LSDA",
        charge: int = 0,
        multiplicity: int = 1,
        n_radial: int = 50,
        n_angular: int = 110,
    ) -> None:
        super().__init__(
            frame, basis=basis, method="uks", xc=xc,
            charge=charge, multiplicity=multiplicity,
            n_radial=n_radial, n_angular=n_angular,
        )

    def __repr__(self) -> str:
        ip = self.integrals
        return (
            f"UKS(n_alpha={self._config.n_alpha}, "
            f"n_beta={self._config.n_beta}, n_basis={ip.n_basis})"
        )


# Backward-compatible result aliases
RHFResult = SCFResult
UHFResult = SCFResult
UKSResult = SCFResult
