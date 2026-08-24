#!/usr/bin/env python3
"""AutoCSF filter-vs-no-filter cost table on one real genomics dataset.

The camera-ready savings table reports space only. This reports the other two
axes -- construction time and query latency -- for AutoCSF alone, so the
comparison is between AutoCSF's filtered index and the plain CSF it was built
from. Both sides use the same native harness and the same CSF implementation,
so the numbers are directly comparable; nothing here crosses method boundaries.

Construction is deterministic and expensive, so it is timed once. A query batch
is 10k random lookups and is noisy at that scale, so it is repeated and reported
as a median with the observed range.
"""

import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from experiments.genomics_comparison.run import (
    dataset_profile, environment, plain_table,
)
from methods.decision_rules import select_filter

BINARY = ROOT / "data/cache/bin/caramel_bench"


def measure(dataset, kind, bpe, hashes, builds, batches):
    if not BINARY.exists():
        subprocess.run([str(ROOT / "experiments/genomics_comparison/build_native.sh")], check=True)
    result = subprocess.run(
        [str(BINARY), str(plain_table(dataset)), kind, str(bpe), str(hashes),
         str(builds), str(batches)],
        capture_output=True, text=True, check=True,
    )
    return json.loads(result.stdout.strip().splitlines()[-1])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", default="rice")
    parser.add_argument("--builds", type=int, default=1)
    parser.add_argument("--query-batches", type=int, default=11)
    parser.add_argument("--output", type=Path, default=ROOT / "results/genomics_comparison")
    args = parser.parse_args()

    manifest, stats = dataset_profile(args.dataset)
    n = manifest["records"]
    selected = select_filter("autocsf", stats)
    if selected is None:
        raise SystemExit(f"AutoCSF declines to filter on {args.dataset}; "
                         "there is no filtered arm to time")
    bpe, hashes = selected["bloom_bits_per_element"], selected["bloom_num_hashes"]

    arms = {
        "no filter": measure(args.dataset, "none", 0, 0, args.builds, args.query_batches),
        "filter": measure(args.dataset, "bloom", bpe, hashes, args.builds, args.query_batches),
    }
    for arm in arms.values():
        arm["bits_per_key"] = arm["serialized_bytes"] * 8 / n

    payload = {
        "schema_version": 1,
        "environment": environment(),
        "dataset": args.dataset,
        "n": n,
        "alpha": manifest["alpha"],
        "method": "autocsf",
        "selected_filter": {"bloom_bits_per_element": bpe, "bloom_num_hashes": hashes},
        "builds": args.builds,
        "query_batches": args.query_batches,
        "arms": arms,
    }
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / f"latency-{args.dataset}.json").write_text(json.dumps(payload, indent=2) + "\n")

    saved = arms["no filter"]["bits_per_key"] - arms["filter"]["bits_per_key"]
    caption = (f"AutoCSF on {args.dataset} (N={n:,}, alpha={manifest['alpha']:.3f}), "
               f"Bloom {bpe} bits/element, {hashes} hashes, saving {saved:.4f} bits/key. "
               f"Construction timed once; query latency is the median of "
               f"{args.query_batches} batches of 10,000 random lookups.")

    md = [caption, "",
          "| | Construction (s) | Query (ns) |",
          "|---|---:|---:|"]
    for label in ("no filter", "filter"):
        arm = arms[label]
        md.append(f"| {label} | {arm['build_seconds']:.2f} | {arm['query_ns']:.1f} |")
    (args.output / f"latency-{args.dataset}.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))
    for label, arm in arms.items():
        print(f"  {label:10s} query range {arm['query_ns_min']:.1f}-{arm['query_ns_max']:.1f} ns, "
              f"{arm['bits_per_key']:.4f} bits/key")


if __name__ == "__main__":
    main()
