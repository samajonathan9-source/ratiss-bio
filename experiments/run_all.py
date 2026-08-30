"""Lanceur global — Expérience n°114, les 7 phases enchaînées."""

import subprocess
import sys
import time
from pathlib import Path

PHASES = [f"experiments/phase{i}_{n}.py" for i, n in enumerate(
    ["invariant", "eth", "bridges", "synthesis", "decoherence", "guard", "verdict"], 1)]

if __name__ == "__main__":
    root = Path(__file__).resolve().parent.parent
    t0 = time.time()
    for p in PHASES:
        print(f"\n{'=' * 60}\n>>> {p}\n{'=' * 60}")
        r = subprocess.run([sys.executable, str(root / p)], cwd=root)
        if r.returncode != 0:
            print(f"ÉCHEC à {p}")
            sys.exit(1)
    print(f"\n✅ Expérience n°114 complète en {time.time() - t0:.1f}s")
