"""Tests for plane-wave DFT Python interface."""

import pytest


@pytest.fixture
def has_pw():
    """Check if PW bindings are available."""
    try:
        from atomiverse._binding import _cpp
        _cpp.PWResult
        return True
    except AttributeError:
        return False


@pytest.mark.skipif(
    "not has_pw",
    reason="ATV_BUILD_PW not enabled in this build",
)
class TestPW:
    """Tests that require the PW backend."""

    def test_solve_si_import(self):
        from atomiverse import solve_si, PWResult

        assert callable(solve_si)
        assert PWResult is not None

    def test_solve_si_toy_potential(self):
        """solve_si with toy local potential returns finite energy."""
        from atomiverse import solve_si

        result = solve_si(ecut_ry=30.0, k_mesh=4, tol=1e-4, max_iter=5)
        assert isinstance(result.energy, float)
        assert isinstance(result.converged, bool)
        assert isinstance(result.iterations, int)
        assert result.iterations > 0

    def test_solve_si_defaults(self):
        """solve_si works with all default arguments."""
        from atomiverse import solve_si

        result = solve_si()
        assert result.energy is not None

    def test_pwresult_repr(self):
        """PWResult.__repr__ is informative."""
        from atomiverse import solve_si

        result = solve_si(ecut_ry=30.0, k_mesh=4, tol=1e-4, max_iter=3)
        s = repr(result)
        assert "PWResult" in s
        assert "Ry" in s


class TestPWNoBackend:
    """Tests that work even without PW compiled."""

    def test_pw_module_importable(self):
        """pw.py should always be importable (wrapper, not binding)."""
        from atomiverse import pw

        assert hasattr(pw, "PWResult")
        assert hasattr(pw, "solve_si")
        assert pw.PWResult is not None
        assert callable(pw.solve_si)
