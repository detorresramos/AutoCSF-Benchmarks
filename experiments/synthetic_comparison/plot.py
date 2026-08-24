#!/usr/bin/env python3
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
rows = list(csv.DictReader((ROOT / "results/synthetic_comparison/data.csv").open()))
series = defaultdict(list)
for row in rows:
    series[(row["distribution"], row["method"])].append((float(row["alpha"]), float(row["bits_saved"])))
styles = {
    "autocsf": ("AutoCSF", "tab:blue", "-"),
    "bcsf": ("BCSF", "tab:orange", "--"),
    "hkp": ("HKP", "tab:green", "-."),
    "vlburr": ("GFT (VL-BuRR)", "tab:red", ":o"),
}
fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
for ax, dist, title in zip(axes, ("uniform_100", "zipfian"), ("Uniform-100", "Zipfian")):
    for method, (label, color, style) in styles.items():
        points = sorted(series[(dist, method)])
        ax.plot([x for x, _ in points], [y for _, y in points], style, color=color, label=label, linewidth=2)
    ax.axhline(0, color="black", linestyle=":", linewidth=.8)
    ax.set_title(title)
    ax.set_xlabel(r"$\alpha$")
    ax.grid(alpha=.25)
axes[0].set_ylabel("Bits/key saved vs no filter")
axes[0].legend()
fig.tight_layout()
out = ROOT / "results/synthetic_comparison/method-comparison.png"
out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(out, dpi=180, bbox_inches="tight")
