#!/usr/bin/env python3
import json
import re
import statistics
import sys

KV = re.compile(r"(\w+)=(\S+)")


def parse(lines, n):
    variants = {"opt": [], "plain": []}
    for line in lines:
        if not line.startswith("RESULT"):
            continue
        row = dict(KV.findall(line))
        name = row.get("storage_name", "")
        if name.endswith("_Opt"):
            variants["opt"].append(row)
        elif name.endswith("_No"):
            variants["plain"].append(row)
    if not variants["opt"] or not variants["plain"]:
        raise ValueError("VL-BuRR log must contain both _Opt and _No RESULT lines")
    opt_bits = statistics.median(float(row["storage_bits"]) for row in variants["opt"])
    plain_bits = statistics.median(float(row["storage_bits"]) for row in variants["plain"])
    return {
        "serialized_bytes": round(opt_bits * n / 8),
        "plain_serialized_bytes": round(plain_bits * n / 8),
        "bits_per_key": opt_bits,
        "baseline_bits_per_key": plain_bits,
        "bits_saved_vs_plain": plain_bits - opt_bits,
        "build_seconds": statistics.median(float(row["construct_ms"]) for row in variants["opt"]) / 1000,
        "query_ns": statistics.median(float(row["query_nanos"]) for row in variants["opt"]),
        "parameters": {
            "implementation": "Filtered-HuffmanCSF_Opt",
            "baseline_implementation": "Filtered-HuffmanCSF_No",
            "reported_storage_bits": opt_bits,
            "baseline_storage_bits": plain_bits,
        },
        "repetitions": len(variants["opt"]),
    }


if __name__ == "__main__":
    try:
        result = parse(sys.stdin, int(sys.argv[1]))
    except ValueError as error:
        raise SystemExit(str(error)) from error
    print("JSON " + json.dumps(result, sort_keys=True))
