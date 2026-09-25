#!/usr/bin/env python3
"""Download predeclared discovery inputs without inspecting their contents."""

from __future__ import annotations

import argparse
import csv
import hashlib
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", type=Path, default=Path("inputs/source-plan.tsv"))
    parser.add_argument("--destination", type=Path, default=Path("inputs/discovery-data"))
    parser.add_argument("--manifest", type=Path, default=Path("inputs/discovery-manifest.tsv"))
    parser.add_argument("--unit-map", type=Path, default=Path("inputs/discovery-unit-map.tsv"))
    args = parser.parse_args()
    args.destination.mkdir(parents=True, exist_ok=True)
    with args.plan.open(encoding="utf-8", newline="") as handle:
        planned = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["role"].startswith("discovery")
        ]
    with args.unit_map.open(encoding="utf-8", newline="") as handle:
        unit_map = {row["accession"]: row for row in csv.DictReader(handle, delimiter="\t")}
    output: list[dict[str, str | int]] = []
    for row in planned:
        target = args.destination / row["file"]
        subprocess.run(
            ["curl", "-L", "--fail", "--retry", "5", "--continue-at", "-", row["url"], "-o", str(target)],
            check=True,
        )
        output.append(
            {
                **row,
                **{key: value for key, value in unit_map[row["accession"]].items() if key != "accession"},
                "retrieved_utc": datetime.now(timezone.utc).isoformat(),
                "bytes": target.stat().st_size,
                "sha256": digest(target),
                "content_opened": "no",
            }
        )
    fields = list(planned[0]) + [
        "biological_unit",
        "condition_map",
        "inclusion_decision",
        "access_terms",
        "retrieved_utc",
        "bytes",
        "sha256",
        "content_opened",
    ]
    with args.manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output)


if __name__ == "__main__":
    main()
