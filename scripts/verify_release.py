"""Verify the SHA-256 integrity manifest for this release."""

from __future__ import annotations

import csv
import hashlib
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "SHA256SUMS.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tree_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    for child in sorted(item for item in path.rglob("*") if item.is_file()):
        relative_path = child.relative_to(path).as_posix()
        digest.update(f"{relative_path}\t{sha256(child)}\n".encode("utf-8"))
    return digest.hexdigest()


def main() -> int:
    with MANIFEST.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    failures = []
    for row in rows:
        path = ROOT / row["relative_path"]
        if row["kind"] == "file":
            observed = sha256(path) if path.is_file() else "MISSING"
        elif row["kind"] == "tree":
            observed = tree_sha256(path) if path.is_dir() else "MISSING"
        else:
            observed = "UNKNOWN_KIND"
        if observed != row["sha256"]:
            failures.append((row["relative_path"], observed, row["sha256"]))
    if failures:
        for relative_path, observed, expected in failures:
            print(f"FAIL {relative_path}: observed={observed} expected={expected}")
        return 1
    print(f"Verified {len(rows)} release files.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
