"""Two-phase benchmark: run PySCF first, then Atomiverse with 2x timeout."""

import importlib.util
import signal
import contextlib
import sys
import time
from pathlib import Path

# Load 04_aspirin_opt.py (numeric prefix prevents normal import)
_spec = importlib.util.spec_from_file_location(
    "aspirin_opt", Path(__file__).parent / "04_aspirin_opt.py"
)
_mod = importlib.util.module_from_spec(_spec)
sys.modules["aspirin_opt"] = _mod
_spec.loader.exec_module(_mod)

make_perturbed_coords = _mod.make_perturbed_coords
run_pyscf = _mod.run_pyscf
run_atomiverse = _mod.run_atomiverse
print_summary = _mod.print_summary

MAX_STEPS = 5
FTOL = 0.05
MEMORY = 5
PYSCF_SCF_TOL = 1e-8


@contextlib.contextmanager
def timeout_ctx(seconds: float):
    def _handler(signum, frame):
        raise TimeoutError(f"Atomiverse exceeded {seconds:.1f}s timeout")
    old = signal.signal(signal.SIGALRM, _handler)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old)


def main() -> None:
    types, ref_coords, perturbed_coords = make_perturbed_coords()

    # Phase 1: PySCF
    print(f"[Phase 1] Running PySCF ({MAX_STEPS} steps)...")
    pyscf = run_pyscf(types, perturbed_coords, MAX_STEPS, FTOL, PYSCF_SCF_TOL)
    print(f"  PySCF done in {pyscf.elapsed_s:.2f}s  (steps={pyscf.iterations})")

    atv_timeout = pyscf.elapsed_s * 2.0
    print(f"\n[Phase 2] Running Atomiverse (timeout={atv_timeout:.1f}s)...")
    t0 = time.perf_counter()
    try:
        with timeout_ctx(atv_timeout):
            atv = run_atomiverse(types, perturbed_coords, MAX_STEPS, FTOL, MEMORY)
    except TimeoutError as e:
        elapsed = time.perf_counter() - t0
        print(f"  TIMEOUT after {elapsed:.2f}s: {e}")
        print(f"  Atomiverse is >{2:.0f}x slower than PySCF for this workload.")
        return

    print_summary(ref_coords, perturbed_coords, atv, pyscf)


if __name__ == "__main__":
    main()
