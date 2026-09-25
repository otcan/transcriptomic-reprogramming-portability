#!/usr/bin/env python3
"""Parse GEO sample SOFT exports into a deterministic metadata table."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_soft(path: Path) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    characteristics: list[str] = []
    descriptions: list[str] = []
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if raw.startswith("^SAMPLE = "):
            if current is not None:
                current["characteristics"] = " | ".join(characteristics)
                current["description"] = " | ".join(descriptions)
                records.append(current)
            current = {"accession": raw.split(" = ", 1)[1]}
            characteristics = []
            descriptions = []
        elif current is None or not raw.startswith("!Sample_"):
            continue
        else:
            key, _, value = raw.partition(" = ")
            field = key.removeprefix("!Sample_")
            if field == "characteristics_ch1":
                characteristics.append(value)
            elif field == "description":
                descriptions.append(value)
            elif field in {"title", "source_name_ch1", "organism_ch1"}:
                current[field] = value
    if current is not None:
        current["characteristics"] = " | ".join(characteristics)
        current["description"] = " | ".join(descriptions)
        records.append(current)
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    fields = [
        "series",
        "accession",
        "title",
        "source_name_ch1",
        "organism_ch1",
        "characteristics",
        "description",
    ]
    rows: list[dict[str, str]] = []
    for path in sorted(args.inputs):
        series = path.name.split(".", 1)[0]
        for record in parse_soft(path):
            rows.append({"series": series, **record})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()

