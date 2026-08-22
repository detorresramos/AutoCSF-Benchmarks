#!/usr/bin/env python3
"""Render the nine bound-validation plots and method-comparison plot."""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
subprocess.run([sys.executable, str(ROOT / "experiments/bound_validation/plot.py")], check=True)
subprocess.run([sys.executable, str(ROOT / "experiments/synthetic_comparison/plot.py")], check=True)
actual = len(list((ROOT / "results/bound_validation/figures/paper").glob("*.png")))
actual += int((ROOT / "results/synthetic_comparison/method-comparison.png").exists())
if actual != 10:
    raise SystemExit(f"expected 10 plots, found {actual}")
print(f"Rendered {actual} paper plots beneath results/")
