#!/usr/bin/env python3
"""Create publication tables that preserve results without source-unit identifiers."""

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
source = ROOT / "results" / "benchmark-b1" / "contrast-matrix.csv"
target = ROOT / "manuscript" / "tables" / "supplementary-table-5-contrast-matrix.csv"

table = pd.read_csv(source)
if len(table) != 184 or table["contrast_id"].duplicated().any():
    raise RuntimeError("Expected 184 uniquely identified frozen contrasts")

method_columns = [
    "ridge_rank_a1",
    "elastic_net_rank",
    "pca50_ridge",
    "meta_effect_projection",
    "gse113957_only_ridge",
    "gse226189_only_ridge",
]
release = table[["family", "dataset", "category", "favourable", *method_columns]].copy()
release.insert(0, "contrast_code", [f"C{i:04d}" for i in range(1, len(release) + 1)])
if release.isna().any().any():
    raise RuntimeError("Release contrast table unexpectedly contains missing values")
target.parent.mkdir(parents=True, exist_ok=True)
release.to_csv(target, index=False)
print(f"wrote {target.relative_to(ROOT)} with {len(release)} rows")
