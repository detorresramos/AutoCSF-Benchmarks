#!/usr/bin/env python3
"""Plot the VL-BuRR filter's empirical bits-per-key-saved vs alpha, in the style
of the companion theory-validation `plot_alpha_sweep`.

"Bits per key saved" = NoFilter storage - Opt storage, i.e. how much the LSF
VL-BuRR filter trick (FilterLengthStrategyOpt) saves over the no-filter CSF.
Positive = filter helps; negative = the dead zone where it hurts.

(Only the VL-BuRR empirical curve is drawn; an AutoCSF curve can be added later.)
"""
import argparse, re, statistics
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RESULT_RE = re.compile(r"(\w+)=(\S+)")
NAME_RE   = re.compile(r"acsf_([a-z0-9]+)_p(\d+)_s(\d+)")

DIST_DISPLAY = {
    "unique": "Unique",
    "zipfian": "Zipfian",
    "uniform10": "Uniform-10",
    "uniform100": "Uniform-100",
    "geometric": "Geometric",
    "twovalue": "Two-value",
}


def parse(path: Path):
    storage = defaultdict(list)  # (dist, alpha, variant) -> bits/key
    with path.open() as f:
        for line in f:
            if not line.startswith("RESULT"):
                continue
            kv = dict(RESULT_RE.findall(line))
            m = NAME_RE.match(kv.get("dataset_name", ""))
            if not m:
                continue
            dist, alpha = m.group(1), int(m.group(2)) / 100
            variant = "Opt" if "_Opt" in kv["storage_name"] else "NoFilter"
            storage[(dist, alpha, variant)].append(float(kv["storage_bits"]))
    return storage


def series(storage, dist):
    alphas = sorted({k[1] for k in storage if k[0] == dist})
    out = []
    for a in alphas:
        nf = storage.get((dist, a, "NoFilter"))
        op = storage.get((dist, a, "Opt"))
        if not (nf and op):
            continue
        saved = statistics.mean(nf) - statistics.mean(op)  # bits/key saved by filter
        sd = statistics.stdev([n - o for n, o in zip(nf, op)]) if len(nf) > 1 else 0.0
        out.append((a, saved, sd))
    return out


def plot_alpha_sweep_vlburr(dist, rows, out_path: Path):
    fig, ax = plt.subplots(figsize=(8, 5))
    alphas = [r[0] for r in rows]
    saved  = [r[1] for r in rows]
    sd     = [r[2] for r in rows]

    ax.plot(alphas, saved, "r.-", linewidth=1.5, markersize=8,
            label="VL-BuRR (empirical)")
    ax.axhline(y=0, color="gray", linestyle="-", alpha=0.5, linewidth=1)

    ax.set_xlabel(r"$\alpha$")
    ax.set_ylabel("Bits per key saved")
    ax.set_xlim(0.5, 1.0)
    ax.set_title(f"VL-BuRR filter — {DIST_DISPLAY.get(dist, dist)}")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out_path), dpi=300, bbox_inches="tight")
    fig.savefig(str(out_path.with_suffix(".pdf")), bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("log", type=Path, help="sweep.log from run_sweep.sh")
    ap.add_argument("--out-dir", type=Path, default=Path("figures/alpha_sweep"))
    args = ap.parse_args()

    storage = parse(args.log)
    for dist in sorted({k[0] for k in storage}):
        rows = series(storage, dist)
        out = args.out_dir / f"alpha_sweep_vlburr_{dist}.png"
        plot_alpha_sweep_vlburr(dist, rows, out)
        print(f"Saved: {out} (+ .pdf)")
        for a, s, sd in rows:
            print(f"  alpha={a:.2f}  bits_saved={s:+.4f} +/- {sd:.4f}")


if __name__ == "__main__":
    main()
