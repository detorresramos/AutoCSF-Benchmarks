#!/usr/bin/env python3
"""Render the ten accepted-paper plots into results/figures."""
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
subprocess.run([sys.executable, str(ROOT / "theory_validation/make_plots.py")], check=True)
subprocess.run([sys.executable, str(ROOT / "method_comparison/make_plot.py")], check=True)
destination = ROOT / "results/figures"
destination.mkdir(parents=True, exist_ok=True)
source = ROOT / "theory_validation/figures/combined"
for path in source.glob("*.png"):
    shutil.copy2(path, destination / path.name)
actual = len(list(destination.glob("*.png")))
if actual != 10:
    raise SystemExit(f"expected 10 plots, found {actual}")
print(f"Rendered {actual} plots in {destination}")
