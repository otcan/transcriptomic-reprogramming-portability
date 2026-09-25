#!/usr/bin/env python3
"""Download the frozen primary temporal-validation inputs and record provenance."""

from __future__ import annotations

import concurrent.futures
import hashlib
import os
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import urlopen

import pandas as pd


SERIES_FILES = {
    "GSE297984": "GSE297984_logCPM_2c_7c.csv.gz",
    "GSE297233": "GSE297233_raw_counts_matrix.csv.gz",
    "GSE300625": "GSE300625_raw_genecounts.csv.gz",
    "GSE304042": "GSE304042_5_ARPE_single_triple_OSK_30-1011800743.csv.gz",
    "GSE304043": "GSE304043_6_mRPE_Old-GFP_Old-GSTA4_merged_featureCounts.csv.gz",
}

SINGLE_CELL_FILES = {
    "GSM8986586": "GSM8986586_GM00731_D0_filtered_feature_bc_matrix.h5",
    "GSM8986587": "GSM8986587_GM00731_D3_filtered_feature_bc_matrix.h5",
    "GSM8986588": "GSM8986588_GM00731_D7_filtered_feature_bc_matrix.h5",
    "GSM8986589": "GSM8986589_GM00731_D10_filtered_feature_bc_matrix.h5",
    "GSM8986590": "GSM8986590_GM23815_D0_filtered_feature_bc_matrix.h5",
    "GSM8986591": "GSM8986591_GM23815_D3_filtered_feature_bc_matrix.h5",
    "GSM8986592": "GSM8986592_GM23815_D7_filtered_feature_bc_matrix.h5",
    "GSM8986593": "GSM8986593_GM23815_D10_filtered_feature_bc_matrix.h5",
}


def series_prefix(accession: str) -> str:
    return accession[:-3] + "nnn"


def download(item: tuple[str, str, str], output: Path) -> dict[str, object]:
    accession, filename, url = item
    target = output / filename
    temporary = output / f".{filename}.partial"
    digest = hashlib.sha256()
    size = 0
    with urlopen(url) as response, temporary.open("wb") as handle:
        while chunk := response.read(1024 * 1024):
            handle.write(chunk)
            digest.update(chunk)
            size += len(chunk)
    os.replace(temporary, target)
    return {
        "accession": accession,
        "filename": filename,
        "url": url,
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "bytes": size,
        "sha256": digest.hexdigest(),
    }


def main() -> None:
    root = Path(".").resolve()
    output = root / "inputs" / "validation-data"
    output.mkdir(parents=True, exist_ok=False)
    items: list[tuple[str, str, str]] = []
    for accession, filename in SERIES_FILES.items():
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{series_prefix(accession)}/{accession}/suppl/{filename}"
        items.append((accession, filename, url))
    for accession in ["GSE297984", "GSE297233", "GSE297234", "GSE300625", "GSE304042", "GSE304043"]:
        filename = f"{accession}_family.soft.gz"
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{series_prefix(accession)}/{accession}/soft/{filename}"
        items.append((accession, filename, url))
    for accession, filename in SINGLE_CELL_FILES.items():
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8986nnn/{accession}/suppl/{filename}"
        items.append((accession, filename, url))

    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
        records = list(executor.map(lambda item: download(item, output), items))
    manifest = pd.DataFrame(records).sort_values(["accession", "filename"])
    manifest.to_csv(root / "inputs" / "validation-manifest.tsv", sep="\t", index=False)
    print(manifest[["accession", "filename", "bytes", "sha256"]].to_string(index=False))


if __name__ == "__main__":
    main()
