#!/usr/bin/env python3
"""Parse a sweep.log produced by run_sweep.sh, aggregate per
(minority_dist, alpha, variant), print a summary table per distribution, and
write one 2-panel figure (storage + query latency) per distribution.
"""
import argparse, re, statistics, sys
from collections import defaultdict
from pathlib import Path


RESULT_RE = re.compile(r'(\w+)=(\S+)')
NAME_RE   = re.compile(r'acsf_([a-z0-9]+)_p(\d+)_s(\d+)')


def parse(path: Path):
    storage = defaultdict(list)   # (dist, alpha, variant) -> bits/key
    query   = defaultdict(list)   # (dist, alpha, variant) -> ns/query
    construct = defaultdict(list)
    entropy = defaultdict(list)   # (dist, alpha) -> empirical entropy (bits)
    with path.open() as f:
        for line in f:
            if not line.startswith("RESULT"):
                continue
            kv = dict(RESULT_RE.findall(line))
            m = NAME_RE.match(kv.get('dataset_name', ''))
            if not m:
                continue
            dist = m.group(1)
            alpha = int(m.group(2)) / 100
            variant = "Opt" if "_Opt" in kv['storage_name'] else "NoFilter"
            storage[(dist, alpha, variant)].append(float(kv['storage_bits']))
            query[(dist, alpha, variant)].append(float(kv['query_nanos']))
            construct[(dist, alpha, variant)].append(float(kv['construct_ms']))
            entropy[(dist, alpha)].append(float(kv['entropy']))
    return storage, query, construct, entropy


def summarize(storage, query, construct, entropy):
    """Returns {dist: [rows]} with rows sorted by alpha."""
    by_dist = defaultdict(list)
    dists = sorted({k[0] for k in storage})
    for dist in dists:
        alphas = sorted({k[1] for k in storage if k[0] == dist})
        for alpha in alphas:
            bn = storage[(dist, alpha, "NoFilter")]
            bo = storage[(dist, alpha, "Opt")]
            qn = query[(dist, alpha, "NoFilter")]
            qo = query[(dist, alpha, "Opt")]
            H = entropy[(dist, alpha)]
            if not (bn and bo):
                continue
            by_dist[dist].append((
                alpha, statistics.mean(H),
                statistics.mean(bn), statistics.stdev(bn) if len(bn) > 1 else 0.0,
                statistics.mean(bo), statistics.stdev(bo) if len(bo) > 1 else 0.0,
                statistics.mean(qn), statistics.stdev(qn) if len(qn) > 1 else 0.0,
                statistics.mean(qo), statistics.stdev(qo) if len(qo) > 1 else 0.0,
            ))
    return by_dist


def print_table(dist, rows) -> None:
    print(f"\n=== minority_dist = {dist} ===")
    h = ("alpha", "H(bits)", "NoFilter bits", "Opt bits", "Δ bits",
         "NoFilter ns", "Opt ns", "Δ ns")
    print(f"{h[0]:>6} {h[1]:>8} {h[2]:>14} {h[3]:>10} {h[4]:>9} {h[5]:>13} {h[6]:>9} {h[7]:>9}")
    print("-" * 90)
    for r in rows:
        p, H, bnm, _bns, bom, _bos, qnm, _qns, qom, _qos = r
        verdict = "filter idle" if abs(bom-bnm) < 0.001 else ("HELPS" if bom < bnm else "HURTS")
        print(f"{p:>6.2f} {H:>8.3f} {bnm:>14.4f} {bom:>10.4f} {bom-bnm:>+9.4f}"
              f" {qnm:>13.2f} {qom:>9.2f} {qom-qnm:>+9.2f}   {verdict}")


def plot(dist, rows, out_path: Path) -> None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    alphas = [r[0] for r in rows]
    H      = [r[1] for r in rows]
    bn_m   = [r[2] for r in rows]; bn_s = [r[3] for r in rows]
    bo_m   = [r[4] for r in rows]; bo_s = [r[5] for r in rows]
    qn_m   = [r[6] for r in rows]; qn_s = [r[7] for r in rows]
    qo_m   = [r[8] for r in rows]; qo_s = [r[9] for r in rows]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
    ax1.plot(alphas, H, 'k--', label='empirical entropy', alpha=0.6, linewidth=1.5)
    ax1.errorbar(alphas, bn_m, yerr=bn_s, marker='o', label='NoFilter',
                 color='#1f77b4', capsize=4, linewidth=2)
    ax1.errorbar(alphas, bo_m, yerr=bo_s, marker='s', label='Opt (filter)',
                 color='#d62728', capsize=4, linewidth=2)
    ax1.axvspan(2/3, 0.8, alpha=0.12, color='gray', label='filter dead zone')
    ax1.set_xlabel(r'$\alpha$ (majority-value fraction)')
    ax1.set_ylabel('storage bits/key')
    ax1.set_title(f'Storage: filter vs no-filter\n(minority_dist = {dist})')
    ax1.legend(loc='upper right', framealpha=0.95)
    ax1.grid(alpha=0.3)

    ax2.errorbar(alphas, qn_m, yerr=qn_s, marker='o', label='NoFilter',
                 color='#1f77b4', capsize=4, linewidth=2)
    ax2.errorbar(alphas, qo_m, yerr=qo_s, marker='s', label='Opt (filter)',
                 color='#d62728', capsize=4, linewidth=2)
    ax2.set_xlabel(r'$\alpha$')
    ax2.set_ylabel('query time (ns)')
    ax2.set_title(f'Query latency: filter vs no-filter\n(minority_dist = {dist})')
    ax2.legend(loc='upper right', framealpha=0.95)
    ax2.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(str(out_path), dpi=150, bbox_inches='tight')
    plt.savefig(str(out_path.with_suffix('.pdf')), bbox_inches='tight')


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("log", type=Path, help="Path to sweep.log from run_sweep.sh")
    ap.add_argument("--out-dir", type=Path, default=Path("figures"),
                    help="Directory for per-distribution figures.")
    args = ap.parse_args()

    storage, query, construct, entropy = parse(args.log)
    if not storage:
        print(f"error: no RESULT lines parsed from {args.log}", file=sys.stderr)
        return 1
    by_dist = summarize(storage, query, construct, entropy)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for dist, rows in by_dist.items():
        print_table(dist, rows)
        out = args.out_dir / f"sweep_{dist}.png"
        plot(dist, rows, out)
        print(f"Figure written: {out} (+ .pdf)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
