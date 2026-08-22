#!/usr/bin/env python3
"""Convert the public genomics TSV format to LSF's frequency-model input."""

import argparse
import json
from pathlib import Path
import struct
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("table", type=Path)
    parser.add_argument("output", type=Path, help="output prefix, without _X/_y")
    parser.add_argument("--manifest", type=Path,
                        help="manifest containing records and histogram (avoids an extra scan)")
    args = parser.parse_args()
    if args.manifest:
        manifest = json.loads(args.manifest.read_text())
        count = int(manifest["records"])
        distinct = sorted(int(value) for value in manifest["histogram"])
    else:
        count, seen = 0, set()
        process = subprocess.Popen(["zstd", "-dc", str(args.table)], stdout=subprocess.PIPE, text=True)
        assert process.stdout is not None
        for line in process.stdout:
            _, value = line.rstrip("\n").split("\t")
            count += 1
            seen.add(int(value))
        if process.wait():
            raise SystemExit("zstd failed")
        distinct = sorted(seen)
    classes = {value: index for index, value in enumerate(distinct)}
    if len(classes) > 65535:
        raise SystemExit(f"LSF uint16 label limit exceeded: {len(classes)} classes")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.with_name(args.output.name + "_X.lrbin").open("wb") as x_handle:
        x_handle.write(struct.pack("<QQ", count, 1))
        zero_chunk = struct.pack("<f", 0.0) * min(count, 1_000_000)
        remaining = count
        while remaining:
            chunk_records = min(remaining, 1_000_000)
            x_handle.write(zero_chunk[: chunk_records * 4])
            remaining -= chunk_records
    with args.output.with_name(args.output.name + "_y.lrbin").open("wb") as y_handle:
        y_handle.write(struct.pack("<H", len(classes)))
        process = subprocess.Popen(["zstd", "-dc", str(args.table)], stdout=subprocess.PIPE, text=True)
        assert process.stdout is not None
        buffer = bytearray()
        observed = 0
        for line in process.stdout:
            _, value = line.rstrip("\n").split("\t")
            buffer += struct.pack("<H", classes[int(value)])
            observed += 1
            if len(buffer) >= 2_000_000:
                y_handle.write(buffer)
                buffer.clear()
        y_handle.write(buffer)
        if process.wait():
            raise SystemExit("zstd failed")
        if observed != count:
            raise SystemExit(f"manifest says {count} rows but table has {observed}")
    print(f"{count} rows, {len(classes)} labels")


if __name__ == "__main__":
    main()
