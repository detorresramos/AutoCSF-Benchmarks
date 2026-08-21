#!/usr/bin/env python3
"""Export sweep.log to tidy CSVs for cross-method comparison plots.

Writes two files:
  <out>_raw.csv      one row per (dataset, variant) RESULT line (per seed).
  <out>_summary.csv  one row per (distribution, alpha, method): seed means/stds
                     of storage bits, query ns, and bits-per-key saved vs the
                     no-filter baseline. Add more methods (e.g. autocsf) as new
                     rows with the same columns to overlay them.

Methods: 'nofilter' (baseline CSF) and 'vlburr' (the LSF FilterLengthStrategyOpt
filter trick = the Opt variant). bits_saved = nofilter_bits - method_bits.
"""
import argparse, csv, re, statistics
from collections import defaultdict
from pathlib import Path

RESULT_RE = re.compile(r"(\w+)=(\S+)")
NAME_RE   = re.compile(r"acsf_([a-z0-9]+)_p(\d+)_s(\d+)")
VARIANT_METHOD = {"NoFilter": "nofilter", "Opt": "vlburr"}


def parse(path: Path):
    rows = []  # raw per-(dataset, variant)
    with path.open() as f:
        for line in f:
            if not line.startswith("RESULT"):
                continue
            kv = dict(RESULT_RE.findall(line))
            m = NAME_RE.match(kv.get("dataset_name", ""))
            if not m:
                continue
            variant = "Opt" if "_Opt" in kv["storage_name"] else "NoFilter"
            rows.append({
                "distribution": m.group(1),
                "alpha": int(m.group(2)) / 100,
                "seed": int(m.group(3)),
                "method": VARIANT_METHOD[variant],
                "entropy_bits": float(kv["entropy"]),
                "storage_bits": float(kv["storage_bits"]),
                "query_ns": float(kv["query_nanos"]),
                "construct_ms": float(kv["construct_ms"]),
            })
    return rows


def mean_std(xs):
    return statistics.mean(xs), (statistics.stdev(xs) if len(xs) > 1 else 0.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("log", type=Path, help="sweep.log from run_sweep.sh")
    ap.add_argument("--out", type=Path, default=Path("results/sweep"),
                    help="Output path prefix (writes <out>_raw.csv and <out>_summary.csv).")
    args = ap.parse_args()

    rows = parse(args.log)
    args.out.parent.mkdir(parents=True, exist_ok=True)

    # --- raw per-seed CSV ---
    raw_path = args.out.with_name(args.out.name + "_raw.csv")
    fields = ["distribution", "alpha", "seed", "method", "entropy_bits",
              "storage_bits", "query_ns", "construct_ms"]
    with raw_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in sorted(rows, key=lambda r: (r["distribution"], r["alpha"], r["seed"], r["method"])):
            w.writerow(r)

    # --- aggregated summary CSV (per distribution, alpha, method) ---
    # index seed-level values for bits_saved (paired against nofilter per seed)
    bits = defaultdict(dict)   # (dist, alpha, seed) -> {method: storage_bits}
    grp  = defaultdict(lambda: defaultdict(list))  # (dist, alpha, method) -> metric -> [vals]
    for r in rows:
        key = (r["distribution"], r["alpha"], r["method"])
        grp[key]["storage_bits"].append(r["storage_bits"])
        grp[key]["query_ns"].append(r["query_ns"])
        grp[key]["construct_ms"].append(r["construct_ms"])
        grp[key]["entropy_bits"].append(r["entropy_bits"])
        bits[(r["distribution"], r["alpha"], r["seed"])][r["method"]] = r["storage_bits"]

    summ_path = args.out.with_name(args.out.name + "_summary.csv")
    with summ_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["distribution", "alpha", "method", "n_seeds", "entropy_bits",
                    "storage_bits_mean", "storage_bits_std",
                    "bits_saved_mean", "bits_saved_std",
                    "query_ns_mean", "query_ns_std",
                    "construct_ms_mean"])
        for (dist, alpha, method) in sorted(grp):
            g = grp[(dist, alpha, method)]
            sb_m, sb_s = mean_std(g["storage_bits"])
            q_m, q_s = mean_std(g["query_ns"])
            # per-seed bits saved vs nofilter baseline
            saved = [b["nofilter"] - b[method]
                     for b in (bits[(dist, alpha, s)] for s in
                               {r["seed"] for r in rows
                                if r["distribution"] == dist and r["alpha"] == alpha})
                     if "nofilter" in b and method in b]
            sv_m, sv_s = mean_std(saved) if saved else (0.0, 0.0)
            w.writerow([dist, alpha, method, len(g["storage_bits"]),
                        f"{statistics.mean(g['entropy_bits']):.6f}",
                        f"{sb_m:.6f}", f"{sb_s:.6f}",
                        f"{sv_m:.6f}", f"{sv_s:.6f}",
                        f"{q_m:.4f}", f"{q_s:.4f}",
                        f"{statistics.mean(g['construct_ms']):.4f}"])

    print(f"Wrote {raw_path} ({len(rows)} rows)")
    print(f"Wrote {summ_path}")


if __name__ == "__main__":
    main()
