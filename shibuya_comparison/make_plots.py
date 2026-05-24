"""Plot generation for prior methods comparison. Reads JSON data, renders plots.

Usage:
    python shibuya_comparison/make_plots.py
"""
import csv
import json
import os
import sys
from collections import defaultdict

import matplotlib.pyplot as plt

_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_dir, ".."))

FIGURES_DIR = os.path.join(_dir, "figures")
DATA_DIR = os.path.join(FIGURES_DIR, "data")
DISTS = ["unique", "zipfian", "uniform_100"]
DIST_LABELS = {"unique": "Unique", "zipfian": "Zipfian", "uniform_100": "Uniform-100"}

# Consistent styling across all methods
METHODS = [
    {
        "key": "our",
        "label": "AutoCSF",
        "color": "tab:blue",
        "linestyle": "-",
        "linewidth": 2,
    },
    {
        "key": "shib",
        "label": "BCSF",
        "color": "tab:orange",
        "linestyle": "--",
        "linewidth": 2,
    },
    {
        "key": "hkp",
        "label": "HKP",
        "color": "tab:green",
        "linestyle": "-.",
        "linewidth": 2,
    },
]


def load_json(path):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def plot_bits_per_key(data):
    N = data["N"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)
    for ax, dist in zip(axes, DISTS):
        rows = data["distributions"][dist]
        alphas = [r["alpha"] for r in rows]

        ax.plot(
            alphas,
            [r["baseline_bpk"] for r in rows],
            "-",
            color="tab:gray",
            linewidth=1.5,
            label="No filter",
        )

        for method in METHODS:
            ax.plot(
                alphas,
                [r[f"{method['key']}_bpk"] for r in rows],
                method["linestyle"],
                color=method["color"],
                linewidth=method["linewidth"],
                label=method["label"],
            )

        ax.set_xlabel(r"$\alpha$")
        ax.set_title(DIST_LABELS[dist])
        ax.legend(loc="best")
        ax.grid(True, alpha=0.3)

    axes[0].set_ylabel("Bits per key")
    fig.suptitle(
        f"Measured bits/key: Bloom filter recommendation comparison (N={N:,})",
        fontsize=14,
        y=1.02,
    )
    plt.tight_layout()
    return fig


def plot_bits_per_key_saved(data):
    N = data["N"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)
    for ax, dist in zip(axes, DISTS):
        rows = data["distributions"][dist]
        alphas = [r["alpha"] for r in rows]

        for method in METHODS:
            ax.plot(
                alphas,
                [r["baseline_bpk"] - r[f"{method['key']}_bpk"] for r in rows],
                method["linestyle"],
                color=method["color"],
                linewidth=method["linewidth"],
                label=method["label"],
            )

        ax.axhline(0, color="black", linewidth=0.5, linestyle=":")
        ax.set_xlabel(r"$\alpha$")
        ax.set_title(DIST_LABELS[dist])
        ax.legend(loc="best")
        ax.grid(True, alpha=0.3)

    axes[0].set_ylabel("Bits/key saved vs no filter")
    fig.suptitle(
        f"Measured bits/key saved: Bloom filter recommendation comparison (N={N:,})",
        fontsize=14,
        y=1.02,
    )
    plt.tight_layout()
    return fig


VLBURR_DOWNLOADS = os.path.expanduser("~/Downloads")
VLBURR_STYLE = {
    "label": "VL-BuRR",
    "color": "tab:red",
    "linestyle": ":",
    "linewidth": 2,
    "marker": "o",
    "markersize": 5,
}

# Panels of the combined VL-BuRR plot. Each entry: dict with json/csv keys + display info.
VLBURR_PANELS = [
    {
        "json_dist_key": "uniform_100",
        "csv_dist": "uniform100",
        "csv_path": os.environ.get(
            "VLBURR_CSV_UNIFORM", os.path.join(VLBURR_DOWNLOADS, "sweep_raw.csv")
        ),
        "n1m_json": "prior_methods_comparison_n1m_uniform.json",
        "title": "Uniform-100",
    },
    {
        "json_dist_key": "zipfian",
        "csv_dist": "zipfian",
        "csv_path": os.environ.get(
            "VLBURR_CSV_ZIPFIAN", os.path.join(VLBURR_DOWNLOADS, "sweep_zipfian_raw.csv")
        ),
        "n1m_json": "prior_methods_comparison_n1m_zipfian.json",
        "title": "Zipfian",
    },
]


def load_vlburr(path, csv_dist):
    """Return {alpha: {"nofilter_bpk": mean, "vlburr_bpk": mean}} for matching rows."""
    if not os.path.exists(path):
        return None
    by_alpha_method = defaultdict(lambda: defaultdict(list))
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            if row["distribution"] != csv_dist:
                continue
            alpha = float(row["alpha"])
            by_alpha_method[alpha][row["method"]].append(float(row["storage_bits"]))

    out = {}
    for alpha, methods in by_alpha_method.items():
        if "nofilter" not in methods or "vlburr" not in methods:
            continue
        out[alpha] = {
            "nofilter_bpk": sum(methods["nofilter"]) / len(methods["nofilter"]),
            "vlburr_bpk": sum(methods["vlburr"]) / len(methods["vlburr"]),
        }
    return dict(sorted(out.items()))


def _draw_vlburr_panel(ax, data, vlburr, dist_key, title, show_legend):
    rows = data["distributions"][dist_key]
    alphas = [r["alpha"] for r in rows]

    for method in METHODS:
        ax.plot(
            alphas,
            [r["baseline_bpk"] - r[f"{method['key']}_bpk"] for r in rows],
            method["linestyle"],
            color=method["color"],
            linewidth=method["linewidth"],
            label=method["label"],
        )

    if vlburr:
        vl_alphas = list(vlburr.keys())
        vl_saved = [v["nofilter_bpk"] - v["vlburr_bpk"] for v in vlburr.values()]
        ax.plot(
            vl_alphas,
            vl_saved,
            VLBURR_STYLE["linestyle"],
            color=VLBURR_STYLE["color"],
            linewidth=VLBURR_STYLE["linewidth"],
            marker=VLBURR_STYLE["marker"],
            markersize=VLBURR_STYLE["markersize"],
            label=VLBURR_STYLE["label"],
        )

    ax.axhline(0, color="black", linewidth=0.5, linestyle=":")
    ax.set_xlabel(r"$\alpha$")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    if show_legend:
        ax.legend(loc="best")


def plot_bits_per_key_saved_vlburr(panels, fallback_data):
    """Render a multi-panel figure of bits/key saved vs no filter, with VL-BuRR overlay.

    `panels` is a list of (panel_cfg, panel_data, vlburr) tuples. `fallback_data`
    is used only to report missing-N=1M info; not plotted directly.
    """
    n_panels = len(panels)
    fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels, 4.5), sharey=True, squeeze=False)
    for i, (cfg, data, vlburr) in enumerate(panels):
        _draw_vlburr_panel(axes[0, i], data, vlburr, cfg["json_dist_key"], cfg["title"], show_legend=(i == 0))
    axes[0, 0].set_ylabel("Bits/key saved vs no filter")
    plt.tight_layout()
    return fig


def main():
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "legend.fontsize": 9,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "lines.linewidth": 1.5,
        }
    )

    path = os.path.join(DATA_DIR, "prior_methods_comparison.json")
    data = load_json(path)
    if data is None:
        print(f"No data found at {path}. Run run_experiments.py first.")
        return

    os.makedirs(FIGURES_DIR, exist_ok=True)
    for name, plot_fn in [
        ("bloom_bits_per_key", plot_bits_per_key),
        ("bloom_bits_per_key_saved", plot_bits_per_key_saved),
    ]:
        fig = plot_fn(data)
        out = os.path.join(FIGURES_DIR, f"{name}.png")
        fig.savefig(out, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {out}")

    panels = []
    for cfg in VLBURR_PANELS:
        vlburr = load_vlburr(cfg["csv_path"], cfg["csv_dist"])
        if vlburr is None:
            print(f"  Skipping {cfg['title']} VL-BuRR panel — no CSV at {cfg['csv_path']}")
            continue
        n1m_path = os.path.join(DATA_DIR, cfg["n1m_json"])
        panel_data = load_json(n1m_path) or data
        if panel_data is data:
            print(
                f"  Note: using N={data['N']:,} data for {cfg['title']} VL-BuRR panel "
                f"(no {cfg['n1m_json']})"
            )
        panels.append((cfg, panel_data, vlburr))

    if panels:
        fig = plot_bits_per_key_saved_vlburr(panels, data)
        out = os.path.join(FIGURES_DIR, "bloom_bits_per_key_saved_vlburr.png")
        fig.savefig(out, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {out}")


if __name__ == "__main__":
    main()
