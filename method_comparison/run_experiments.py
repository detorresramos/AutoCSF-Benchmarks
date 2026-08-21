#!/usr/bin/env python3
"""Reproduce the accepted-paper four-method synthetic comparison."""

import argparse
import csv
from pathlib import Path
import re
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from baselines.methods import CSFFilter
from shared.data_gen import gen_alpha_values, gen_keys

ALPHAS = [0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.99]
DISTS = ["uniform_100", "zipfian"]
LOCAL = {
    "autocsf": lambda: CSFFilter("bloom", "optimal"),
    "bcsf": lambda: CSFFilter("bloom", "shibuya"),
    "hkp": lambda: CSFFilter("bloom", "hkp"),
}
KV = re.compile(r"(\w+)=(\S+)")


def serialized_bpk(method, keys, values):
    structure = method.construct(keys, values)
    return method.measure_memory_from_structure(structure)["serialized"] * 8 / len(keys)


def run_local(n):
    records = []
    keys = gen_keys(n)
    for dist in DISTS:
        for alpha in ALPHAS:
            values = gen_alpha_values(n, alpha, seed=42, minority_dist=dist)
            baseline = CSFFilter("bloom", "hkp")
            baseline.construct = lambda k, v: __import__("carameldb").Caramel(k, v, prefilter=None, verbose=False)
            baseline_bpk = serialized_bpk(baseline, keys, values)
            for name, factory in LOCAL.items():
                bpk = serialized_bpk(factory(), keys, values)
                records.append({"distribution": dist, "alpha": alpha, "method": name, "bits_saved": baseline_bpk - bpk})
    return records


def run_vlburr(n, work):
    binary = ROOT / "data" / "cache" / "bin" / "ribbon_learned_bench"
    if not binary.exists():
        subprocess.run([str(ROOT / "vlburr" / "build.sh")], check=True)
    lrdata = work / "lrdata"
    subprocess.run([sys.executable, str(ROOT / "vlburr/scripts/gen_datasets.py"), "--out", str(lrdata), "--n", str(n), "--minority-dist", *DISTS, "--p-max", *map(str, ALPHAS), "--seeds", "0"], check=True)
    records = []
    for labels in sorted(lrdata.glob("*_y.lrbin")):
        name = labels.name.removesuffix("_y.lrbin")
        run = subprocess.run([str(binary), "-r", str(lrdata) + "/", "-c", "CSF", "-d", name], capture_output=True, text=True, check=True)
        variants = {}
        for line in run.stdout.splitlines():
            if line.startswith("RESULT"):
                values = dict(KV.findall(line))
                variants["opt" if "_Opt" in values["storage_name"] else "plain"] = float(values["storage_bits"])
        match = re.match(r"acsf_(uniform100|zipfian)_p(\d+)_s0", name)
        if not match or set(variants) != {"plain", "opt"}:
            raise RuntimeError(f"incomplete VL-BuRR output for {name}")
        dist = "uniform_100" if match.group(1) == "uniform100" else "zipfian"
        records.append({"distribution": dist, "alpha": int(match.group(2)) / 100, "method": "vlburr", "bits_saved": variants["plain"] - variants["opt"]})
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=1_000_000)
    parser.add_argument("--output", type=Path, default=ROOT / "method_comparison/data.csv")
    parser.add_argument("--component", choices=("all", "local", "vlburr"), default="all")
    parser.add_argument("--append", action="store_true")
    args = parser.parse_args()
    records = []
    if args.component in ("all", "local"):
        records.extend(run_local(args.n))
    if args.component in ("all", "vlburr"):
        records.extend(run_vlburr(args.n, ROOT / "data/cache/method-comparison"))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    mode = "a" if args.append else "w"
    needs_header = not args.append or not args.output.exists()
    with args.output.open(mode, newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["distribution", "alpha", "method", "bits_saved"])
        if needs_header:
            writer.writeheader()
        writer.writerows(sorted(records, key=lambda r: (r["distribution"], r["method"], r["alpha"])))


if __name__ == "__main__":
    main()
