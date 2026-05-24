"""Generate paper-quality plots and tables for AutoCSF (VLDB 2026).

Reads existing benchmark data from baselines/figures/data/.
No rerunning needed — uses JSON files from previous benchmark runs.

Generates:
    baselines/figures/paper_pareto.png
    baselines/figures/paper_memory_vs_alpha.png
    baselines/figures/tables.md

Usage:
    python baselines/paper_plots.py
"""

import glob
import json
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

_dir = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIR = os.path.join(_dir, "figures")
DATA_DIR = os.path.join(FIGURES_DIR, "data")

# ── Layout ────────────────────────────────────────────────────────────────────

DISTS = ["uniform_100", "zipfian", "unique"]
DIST_LABEL = {"uniform_100": "Uniform-100", "zipfian": "Zipfian", "unique": "Unique"}
ALPHAS = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]
HIGHLIGHT_ALPHAS = {0.5, 0.8, 0.95}
TABLE_ALPHAS = [0.5, 0.8, 0.95]
N = 100_000

GENOMICS_FILES = [
    ("baselines_ecoli_sakai_k15_binary_fuse.json", "E. coli"),
    ("baselines_srr10211353_k15_binary_fuse.json", "SRR"),
    ("baselines_celegans_k15_binary_fuse.json", "C. elegans"),
]

# ── Method styling (order = legend order) ─────────────────────────────────────

METHODS = [
    "csf_filter_optimal_binary_fuse",
    "csf_filter_shibuya_bloom",
    "cpp_hash_table",
    "java_mph",
]

STYLES = {
    "csf_filter_optimal_binary_fuse": dict(
        label="AutoCSF", color="#2166ac", marker="o", lw=2.0, ms=7
    ),
    "csf_filter_shibuya_bloom": dict(
        label="BCSF (Shibuya)", color="#67a9cf", marker="^", lw=1.5, ms=6
    ),
    "cpp_hash_table": dict(
        label="C++ Hash Table", color="#d6604d", marker="D", lw=1.5, ms=6
    ),
    "java_mph": dict(
        label="MPH Table", color="#762a83", marker="P", lw=1.5, ms=6
    ),
}

# ── Helpers ───────────────────────────────────────────────────────────────────


def load_synthetic():
    """Load all synthetic experiment files. Returns dict[(dist, alpha)] -> data."""
    exps = {}
    for dist in DISTS:
        for alpha in ALPHAS:
            path = os.path.join(
                DATA_DIR, f"baselines_n{N}_a{alpha}_{dist}_binary_fuse.json"
            )
            if os.path.exists(path):
                with open(path) as f:
                    exps[(dist, alpha)] = json.load(f)
    return exps


def _pick_best_learned_tokenizer(base_fname, label):
    """Load all tokenizer variants for a genomics dataset, merge best learned CSF.

    Returns (label, merged_data) where the learned CSF result is the one with
    lowest total bpk across all tokenizer variants.
    """
    base_path = os.path.join(DATA_DIR, base_fname)
    if not os.path.exists(base_path):
        return None

    with open(base_path) as f:
        base_data = json.load(f)

    # Find tokenizer variant files: baselines_{dataset}_binary_fuse_{tok}.json
    stem = base_fname.replace(".json", "")  # e.g. baselines_ecoli_sakai_k15_binary_fuse
    variant_files = glob.glob(os.path.join(DATA_DIR, f"{stem}_*.json"))

    learned_method = "lsf_ours_filtered-huffman_opt"
    best_learned = find(base_data, learned_method)
    best_bpk = float("inf")
    if best_learned:
        best_bpk = mem_bpk(best_learned, base_data["dataset"]["N"])

    for vf in variant_files:
        with open(vf) as f:
            vdata = json.load(f)
        # Tokenizer variants have method names like lsf_lsf_learned_kmer_ordinal
        # but results_to_json always uses "lsf_ours_filtered-huffman_opt"
        vr = find(vdata, learned_method)
        if not vr:
            continue
        vbpk = mem_bpk(vr, vdata["dataset"]["N"])
        if vbpk < best_bpk:
            best_bpk = vbpk
            best_learned = vr

    # Replace the learned CSF result in base_data with the best one
    if best_learned:
        base_data["results"] = [
            r for r in base_data["results"] if r["method"] != learned_method
        ]
        base_data["results"].append(best_learned)

    return (label, base_data)


def load_genomics():
    """Load genomics experiment files, picking best tokenizer for learned CSF.

    Returns list of (label, data) tuples.
    """
    out = []
    for fname, label in GENOMICS_FILES:
        result = _pick_best_learned_tokenizer(fname, label)
        if result:
            out.append(result)
        else:
            print(f"  Warning: {fname} not found, skipping")
    return out


def mem_bpk(result, n):
    """Extract memory in bits per key."""
    m = result["memory"]
    raw = m.get("theoretical") or m.get("serialized") or m.get("serialized_bytes")
    return raw * 8 / n


def lat_ns(result):
    """Extract inference latency in ns."""
    return result["inference_ns"]["mean"]


def find(exp, method):
    return next((r for r in exp["results"] if r["method"] == method), None)


def _collect(exps, dist, method, alphas=ALPHAS):
    """Collect (alpha, memory_bpk, latency_ns) tuples for a method+dist."""
    out = []
    for a in alphas:
        e = exps.get((dist, a))
        if not e:
            continue
        r = find(e, method)
        if not r:
            continue
        out.append((a, mem_bpk(r, e["dataset"]["N"]), lat_ns(r)))
    return out


def _format_n(n):
    if n >= 1_000_000:
        return f"{n // 1_000_000}M"
    elif n >= 1_000:
        return f"{n // 1_000}k"
    return str(n)


# ── Shared legend ─────────────────────────────────────────────────────────────


def _add_legend(fig, methods):
    handles = [
        Line2D(
            [0],
            [0],
            marker=STYLES[m]["marker"],
            color=STYLES[m]["color"],
            linewidth=STYLES[m]["lw"],
            markersize=5,
            label=STYLES[m]["label"],
        )
        for m in methods
        if m in STYLES
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=min(len(handles), 7),
        fontsize=8,
        bbox_to_anchor=(0.5, -0.02),
        frameon=True,
        edgecolor="#cccccc",
    )


# ── Figure 1: Pareto (Memory vs Latency) — 1×4 with genomics ────────────────


def plot_pareto(exps, genomics):
    fig, axes = plt.subplots(1, 4, figsize=(14, 4.0))

    # --- Synthetic subplots (1-3) ---
    for i, dist in enumerate(DISTS):
        ax = axes[i]
        for mkey in METHODS:
            s = STYLES[mkey]
            pts = _collect(exps, dist, mkey)
            if not pts:
                continue
            alps, ms, ls = zip(*pts)

            # Trajectory line
            if len(ms) > 1:
                ax.plot(
                    ms, ls,
                    color=s["color"], lw=s["lw"] * 0.7, alpha=0.5, zorder=2,
                )

            # Points — larger at highlight alphas
            for a, m, l in zip(alps, ms, ls):
                sz = s["ms"] if a in HIGHLIGHT_ALPHAS else s["ms"] * 0.6
                ax.scatter(
                    m, l,
                    s=sz**2, marker=s["marker"], color=s["color"],
                    edgecolors="white", linewidths=0.4, zorder=3,
                )

            # Annotate alpha at endpoints for AutoCSF only
            if mkey == "csf_filter_optimal_binary_fuse" and len(alps) >= 2:
                ax.annotate(
                    r"$\alpha\!=\!0.5$",
                    (ms[0], ls[0]),
                    fontsize=6, color=s["color"],
                    textcoords="offset points", xytext=(5, 3),
                )
                ax.annotate(
                    r"$\alpha\!=\!0.99$",
                    (ms[-1], ls[-1]),
                    fontsize=6, color=s["color"],
                    textcoords="offset points", xytext=(5, -2),
                )

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Memory (bits/key)")
        ax.set_title(DIST_LABEL[dist])
        ax.grid(True, alpha=0.15, which="both")
        if i == 0:
            ax.set_ylabel("Avg. inference latency (ns)")

    # --- Genomics subplot (4th) ---
    ax = axes[3]
    for mkey in METHODS:
        s = STYLES[mkey]
        for g_label, g_data in genomics:
            r = find(g_data, mkey)
            if not r:
                continue
            n = g_data["dataset"]["N"]
            m = mem_bpk(r, n)
            l = lat_ns(r)
            ax.scatter(
                m, l,
                s=s["ms"] ** 2, marker=s["marker"], color=s["color"],
                edgecolors="white", linewidths=0.4, zorder=3,
            )
            ax.annotate(
                f"${g_label}$",
                (m, l),
                fontsize=5.5,
                textcoords="offset points", xytext=(4, 2),
                color=s["color"], alpha=0.7,
            )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Memory (bits/key)")
    ax.set_title("Genomics")
    ax.grid(True, alpha=0.15, which="both")

    _add_legend(fig, METHODS)
    plt.tight_layout(rect=[0, 0.07, 1, 1])
    return fig


# ── Figure 2: Memory vs Alpha ────────────────────────────────────────────────


def plot_memory_vs_alpha(exps):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), sharey=False)

    for i, dist in enumerate(DISTS):
        ax = axes[i]
        for mkey in METHODS:
            s = STYLES[mkey]
            pts = _collect(exps, dist, mkey)
            if not pts:
                continue
            alps, ms, _ = zip(*pts)
            ax.plot(
                alps, ms,
                color=s["color"], marker=s["marker"],
                markersize=s["ms"], linewidth=s["lw"],
            )

        ax.set_xlabel(r"$\alpha$ (majority fraction)")
        ax.set_title(DIST_LABEL[dist])
        ax.grid(True, alpha=0.15)
        if i == 0:
            ax.set_ylabel("Memory (bits/key)")

    _add_legend(fig, METHODS)
    plt.tight_layout(rect=[0, 0.07, 1, 1])
    return fig


# ── Tables (Markdown) ────────────────────────────────────────────────────────


def generate_synthetic_table(exps):
    """Combined markdown table: methods × (alpha → memory, latency, construction time)."""
    lines = []
    for dist in DISTS:
        lines.append(f"\n### {DIST_LABEL[dist]}\n")

        header = "| Method |"
        sep = "| --- |"
        for a in TABLE_ALPHAS:
            header += f" bpk (α={a}) | ns | build (s) |"
            sep += " --- | --- | --- |"
        lines.append(header)
        lines.append(sep)

        for mkey in METHODS:
            s = STYLES[mkey]
            row = f"| {s['label']} |"
            for a in TABLE_ALPHAS:
                e = exps.get((dist, a))
                if not e:
                    row += " --- | --- | --- |"
                    continue
                r = find(e, mkey)
                if not r:
                    row += " --- | --- | --- |"
                    continue
                m = mem_bpk(r, e["dataset"]["N"])
                l = lat_ns(r)
                ct = r["construction_time_s"]
                row += f" {m:.1f} | {l:.0f} | {ct:.3f} |"
            lines.append(row)
        lines.append("")
    return "\n".join(lines)


def generate_genomics_table(genomics):
    """Genomics markdown table: methods × datasets → memory, latency, construction time."""
    if not genomics:
        return ""

    lines = []
    lines.append("\n### Genomics\n")

    header = "| Method |"
    sep = "| --- |"
    for g_label, g_data in genomics:
        ds = g_data["dataset"]
        n_str = _format_n(ds["N"])
        alpha_str = f"{ds['alpha']:.2f}"
        header += f" bpk ({g_label}, n={n_str}, α={alpha_str}) | ns | build (s) |"
        sep += " --- | --- | --- |"
    lines.append(header)
    lines.append(sep)

    for mkey in METHODS:
        s = STYLES[mkey]
        row = f"| {s['label']} |"
        for g_label, g_data in genomics:
            r = find(g_data, mkey)
            if not r:
                row += " --- | --- | --- |"
                continue
            n = g_data["dataset"]["N"]
            m = mem_bpk(r, n)
            l = lat_ns(r)
            ct = r["construction_time_s"]
            row += f" {m:.2f} | {l:.0f} | {ct:.1f} |"
        lines.append(row)
    lines.append("")
    return "\n".join(lines)


# ── Main ──────────────────────────────────────────────────────────────────────


def main():
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 9,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "legend.fontsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
        }
    )

    print("Loading data...")
    exps = load_synthetic()
    print(f"  Loaded {len(exps)} synthetic experiments")
    genomics = load_genomics()
    print(f"  Loaded {len(genomics)} genomics datasets")

    os.makedirs(FIGURES_DIR, exist_ok=True)

    # --- Plots ---
    print("\nFigure 1: Pareto frontier (log-log, 1x4)...")
    fig = plot_pareto(exps, genomics)
    path = os.path.join(FIGURES_DIR, "paper_pareto.png")
    fig.savefig(path)
    print(f"  {path}")
    plt.close(fig)

    print("\nFigure 2: Memory vs alpha...")
    fig = plot_memory_vs_alpha(exps)
    path = os.path.join(FIGURES_DIR, "paper_memory_vs_alpha.png")
    fig.savefig(path)
    print(f"  {path}")
    plt.close(fig)

    # --- Tables ---
    print("\nGenerating tables.md...")
    md = "# Baseline Results\n"
    md += generate_synthetic_table(exps)
    md += "\n"
    md += generate_genomics_table(genomics)

    tables_path = os.path.join(FIGURES_DIR, "tables.md")
    with open(tables_path, "w") as f:
        f.write(md)
    print(f"  {tables_path}")

    # Print to stdout too
    print(md)


if __name__ == "__main__":
    main()
