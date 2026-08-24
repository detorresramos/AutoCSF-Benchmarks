#!/usr/bin/env python3
"""Compatibility entry point for legacy Codex stop-hook configurations.

Artifact verification moved to ``scripts/verify.py``.  The old stop hook only
needs a valid JSON response; explicit verification should use ``make verify``.
"""

import argparse


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hook", action="store_true")
    args = parser.parse_args()
    if args.hook:
        print("{}")
        return
    parser.error("completion_check.py was replaced by scripts/verify.py")


if __name__ == "__main__":
    main()
