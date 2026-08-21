#!/usr/bin/env python3
"""Download, generate, and validate AutoCSF genomics count tables."""

import argparse
import collections
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = Path(__file__).with_name("sources.json")
DATA_ROOT = ROOT / "data"


def config():
    with CONFIG_PATH.open() as handle:
        return json.load(handle)


def names(args, cfg):
    if args.dataset:
        if args.dataset not in cfg["datasets"]:
            raise SystemExit(f"unknown dataset: {args.dataset}")
        return [args.dataset]
    return list(cfg["datasets"])


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def decompressed_sha256(path):
    """Hash canonical TSV bytes, independent of the zstd encoder version."""
    digest = hashlib.sha256()
    process = subprocess.Popen(["zstd", "-dc", str(path)], stdout=subprocess.PIPE)
    assert process.stdout is not None
    for block in iter(lambda: process.stdout.read(8 * 1024 * 1024), b""):
        digest.update(block)
    if process.wait():
        raise RuntimeError(f"zstd failed while hashing {path}")
    return digest.hexdigest()


def download(args, cfg):
    try:
        from huggingface_hub import snapshot_download
    except ImportError as exc:
        raise SystemExit("install requirements.txt before downloading") from exc
    repo = os.environ.get("AUTOCSF_DATASET_REPO", cfg["dataset_repo"])
    allow = []
    for name in names(args, cfg):
        allow.extend([f"processed/{name}_k15.tsv.zst", f"manifests/{name}_k15.json"])
        if args.include_raw:
            allow.extend(cfg["datasets"][name]["raw"])
    snapshot_download(
        repo_id=repo,
        repo_type="dataset",
        local_dir=DATA_ROOT,
        allow_patterns=allow,
    )


def fetch(args, cfg):
    """Download original archives from the authoritative public sources."""
    for name in names(args, cfg):
        entry = cfg["datasets"][name]
        if len(entry.get("urls", [])) != len(entry["raw"]):
            raise SystemExit(f"{name}: raw paths and URLs do not match")
        for relative, url in zip(entry["raw"], entry["urls"]):
            destination = DATA_ROOT / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            if destination.exists() and not args.force:
                print(f"exists: {destination}")
                continue
            temporary = destination.with_suffix(destination.suffix + ".part")
            print(f"download: {url}")
            urllib.request.urlretrieve(url, temporary)
            temporary.replace(destination)
            print(f"saved: {destination} ({sha256(destination)})")


def stage(args, cfg):
    destination = args.output.resolve()
    for directory in ("raw", "processed", "manifests"):
        source = DATA_ROOT / directory
        if source.exists():
            shutil.copytree(source, destination / directory, dirs_exist_ok=True)
    shutil.copy2(Path(__file__).with_name("DATASET_CARD.md"), destination / "README.md")
    shutil.copy2(CONFIG_PATH, destination / "sources.json")
    checksums = []
    for path in sorted(p for p in destination.rglob("*") if p.is_file() and p.name != "SHA256SUMS"):
        checksums.append(f"{sha256(path)}  {path.relative_to(destination)}")
    (destination / "SHA256SUMS").write_text("\n".join(checksums) + "\n")
    print(destination)


def generate(args, cfg):
    name = names(args, cfg)[0]
    entry = cfg["datasets"][name]
    if not entry["raw"]:
        raise SystemExit(f"raw source for {name} is not resolved in sources.json")
    for executable in ("c++", "zstd"):
        if not shutil.which(executable):
            raise SystemExit(f"required executable not found: {executable}")
    raw_paths = [DATA_ROOT / path for path in entry["raw"]]
    missing = [str(path) for path in raw_paths if not path.exists()]
    if missing:
        raise SystemExit("missing raw input(s): " + ", ".join(missing))
    work = DATA_ROOT / "cache" / "kmc" / name
    work.mkdir(parents=True, exist_ok=True)
    processed = DATA_ROOT / "processed"
    processed.mkdir(parents=True, exist_ok=True)
    binary = DATA_ROOT / "cache" / "bin" / "count_kmers"
    binary.parent.mkdir(parents=True, exist_ok=True)
    source = Path(__file__).with_name("count_kmers.cpp")
    if not binary.exists() or source.stat().st_mtime > binary.stat().st_mtime:
        subprocess.run(["c++", "-O3", "-std=c++17", str(source), "-lz", "-o", str(binary)], check=True)
    stats_path = work / "stats.tsv"
    output = processed / f"{name}_k15.tsv.zst"
    with subprocess.Popen(
        [str(binary), str(stats_path), entry["format"], *map(str, raw_paths)], stdout=subprocess.PIPE
    ) as counter:
        assert counter.stdout is not None
        compressor = subprocess.run(
            ["zstd", "-T0", "-19", "-f", "-o", str(output)],
            stdin=counter.stdout,
        )
        counter.stdout.close()
        counter_status = counter.wait()
    if counter_status or compressor.returncode:
        raise RuntimeError("k-mer counting or compression failed")
    stats = {}
    for line in stats_path.read_text().splitlines():
        key, value = line.split("\t")
        stats[key] = int(value)
    manifests = DATA_ROOT / "manifests"
    manifests.mkdir(parents=True, exist_ok=True)
    manifest = {
        "name": name,
        "accession": entry["accession"],
        "k": 15,
        "canonical": False,
        "case_policy": "normalize-ASCII-ACGT-to-uppercase",
        "ambiguous_base_policy": "omit-window",
        "records": stats["records"],
        "distinct_values": stats["distinct_values"],
        "modal_value": stats["modal_value"],
        "alpha": stats["modal_records"] / stats["records"],
        "sha256": sha256(output),
        "content_sha256": decompressed_sha256(output),
        "raw": [
            {"path": str(path.relative_to(DATA_ROOT)), "sha256": sha256(path)}
            for path in raw_paths
        ],
        "generator": "datasets/count_kmers.cpp (dense 2-bit forward-strand k=15 counter)",
    }
    manifest_path = manifests / f"{name}_k15.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(output)
    print(manifest_path)


def iter_rows(path):
    process = subprocess.Popen(["zstd", "-dc", str(path)], stdout=subprocess.PIPE, text=True)
    assert process.stdout is not None
    for line_number, line in enumerate(process.stdout, 1):
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 2:
            process.kill()
            raise ValueError(f"{path}:{line_number}: expected two tab-separated fields")
        kmer, value = fields
        if len(kmer) != 15 or any(base not in "ACGT" for base in kmer):
            process.kill()
            raise ValueError(f"{path}:{line_number}: invalid 15-mer {kmer!r}")
        yield kmer, int(value)
    if process.wait() != 0:
        raise RuntimeError(f"zstd failed while reading {path}")


def validate_one(name, entry, manifest_only=False):
    path = DATA_ROOT / "processed" / f"{name}_k15.tsv.zst"
    manifest_path = DATA_ROOT / "manifests" / f"{name}_k15.json"
    if not path.exists() or not manifest_path.exists():
        raise FileNotFoundError(f"missing processed table or manifest for {name}")
    with manifest_path.open() as handle:
        manifest = json.load(handle)
    actual_hash = sha256(path)
    if manifest.get("sha256") != actual_hash:
        raise ValueError(f"{name}: SHA-256 mismatch")
    if manifest_only:
        print(f"{name}: checksum ok")
        return
    if manifest.get("content_sha256"):
        if manifest["content_sha256"] != decompressed_sha256(path):
            raise ValueError(f"{name}: decompressed content SHA-256 mismatch")
    histogram = collections.Counter()
    records = 0
    previous = None
    for kmer, value in iter_rows(path):
        if previous is not None and kmer <= previous:
            raise ValueError(f"{name}: table is not strictly sorted at {kmer}")
        previous = kmer
        histogram[value] += 1
        records += 1
    modal_value, modal_records = min(histogram.items(), key=lambda item: (-item[1], item[0]))
    observed = {
        "records": records,
        "distinct_values": len(histogram),
        "modal_value": modal_value,
        "alpha": modal_records / records,
    }
    for key in ("records", "distinct_values", "modal_value"):
        if key in manifest and manifest[key] != observed[key]:
            raise ValueError(f"{name}: manifest {key} does not match table")
    expected = entry.get("expected", {})
    for key in ("records", "distinct_values", "modal_value"):
        if key in expected and expected[key] != observed[key]:
            raise ValueError(f"{name}: expected {key}={expected[key]}, observed {observed[key]}")
    if "alpha" in expected and abs(expected["alpha"] - observed["alpha"]) > 0.001:
        raise ValueError(f"{name}: expected alpha≈{expected['alpha']}, observed {observed['alpha']}")
    print(f"{name}: {records:,} records, alpha={observed['alpha']:.6f}, {len(histogram)} values")


def profile_one(name):
    path = DATA_ROOT / "processed" / f"{name}_k15.tsv.zst"
    manifest_path = DATA_ROOT / "manifests" / f"{name}_k15.json"
    histogram = collections.Counter(value for _, value in iter_rows(path))
    manifest = json.loads(manifest_path.read_text())
    manifest["histogram"] = {str(value): frequency for value, frequency in sorted(histogram.items())}
    manifest["content_sha256"] = decompressed_sha256(path)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"{name}: recorded {len(histogram)} histogram entries")


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    list_parser = sub.add_parser("list")
    list_parser.add_argument("--dataset")
    download_parser = sub.add_parser("download")
    download_parser.add_argument("--dataset")
    download_parser.add_argument("--include-raw", action="store_true")
    fetch_parser = sub.add_parser("fetch")
    fetch_parser.add_argument("--dataset")
    fetch_parser.add_argument("--force", action="store_true")
    generate_parser = sub.add_parser("generate")
    generate_parser.add_argument("--dataset", required=True)
    validate_parser = sub.add_parser("validate")
    validate_parser.add_argument("--dataset")
    validate_parser.add_argument("--manifest-only", action="store_true")
    stage_parser = sub.add_parser("stage")
    stage_parser.add_argument("--output", type=Path, default=ROOT / "results" / "huggingface" / "autocsf-genomics")
    profile_parser = sub.add_parser("profile")
    profile_parser.add_argument("--dataset")
    args = parser.parse_args()
    cfg = config()
    if args.command == "list":
        for name in names(args, cfg):
            print(name, cfg["datasets"][name]["accession"])
    elif args.command == "download":
        download(args, cfg)
    elif args.command == "fetch":
        fetch(args, cfg)
    elif args.command == "generate":
        generate(args, cfg)
    elif args.command == "validate":
        for name in names(args, cfg):
            validate_one(name, cfg["datasets"][name], args.manifest_only)
    elif args.command == "stage":
        stage(args, cfg)
    elif args.command == "profile":
        for name in names(args, cfg):
            profile_one(name)


if __name__ == "__main__":
    main()
