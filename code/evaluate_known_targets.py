#!/usr/bin/env python3
"""Evaluate the frozen five-target benchmark against an abundance-matched null."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from core import read_gtf_gene_map
from run_discovery import load_gse113957, load_gse226189, load_metadata


TARGETS = ["GSTA4", "DDX21", "TOMM70A", "SERBP1", "PHGDH"]


def assign_expression_deciles(abundance: pd.Series) -> pd.Series:
    ordered = abundance.rename_axis("gene").reset_index(name="abundance").sort_values(
        ["abundance", "gene"], kind="mergesort"
    )
    ordered["decile"] = pd.qcut(np.arange(len(ordered)), 10, labels=False)
    return ordered.set_index("gene")["decile"].sort_index().astype(int)


def matched_null(
    percentiles: pd.Series,
    deciles: pd.Series,
    targets: list[str],
    n_resamples: int = 100_000,
    seed: int = 271828,
) -> tuple[float, np.ndarray]:
    missing = [target for target in targets if target not in percentiles.index]
    if missing:
        raise ValueError(f"Targets absent from frozen rank: {missing}")
    target_set = set(targets)
    pool_genes = {
        decile: np.array(sorted(gene for gene in deciles.index[deciles == decile] if gene not in target_set))
        for decile in sorted(deciles.unique())
    }
    pools = {decile: percentiles.loc[genes].to_numpy(dtype=float) for decile, genes in pool_genes.items()}
    rng = np.random.default_rng(seed)
    null = np.empty(n_resamples, dtype=float)
    target_deciles = deciles.loc[targets].tolist()
    decile_counts = pd.Series(target_deciles).value_counts().sort_index().to_dict()
    for iteration in range(n_resamples):
        selected = np.concatenate(
            [
            rng.choice(pools[decile], size=int(count), replace=False)
            for decile, count in decile_counts.items()
            ]
        )
        null[iteration] = selected.mean()
    observed = float(percentiles.loc[targets].mean())
    return observed, null


def main() -> None:
    root = Path(".").resolve()
    metadata = load_metadata(root / "receipts" / "sample-metadata.tsv")
    human_map = read_gtf_gene_map(root / "inputs" / "annotations" / "Homo_sapiens.GRCh38.111.gtf.gz")
    reference_a, _, _ = load_gse113957(root / "inputs" / "discovery-data" / "GSE113957_fpkm.txt.gz", metadata)
    reference_b, _, _ = load_gse226189(root / "inputs" / "discovery-data" / "GSE226189_RAW.tar", metadata, human_map)
    rank = pd.read_csv(root / "results" / "discovery" / "candidate-rank.csv", index_col="gene")

    cohort_a = reference_a.rank(axis=0, method="average", pct=True).median(axis=1)
    cohort_b = reference_b.rank(axis=0, method="average", pct=True).median(axis=1)
    abundance = pd.concat([cohort_a.rename("a"), cohort_b.rename("b")], axis=1).mean(axis=1)
    abundance = abundance.reindex(rank.index)
    if abundance.isna().any():
        raise ValueError("Candidate abundance is incomplete")
    deciles = assign_expression_deciles(abundance)
    evaluable = [target for target in TARGETS if target in rank.index]
    missing = [target for target in TARGETS if target not in rank.index]
    observed, null = matched_null(rank["percentile"], deciles, evaluable)
    p_value = float((1 + np.sum(null >= observed)) / (len(null) + 1))

    target_table = rank.reindex(TARGETS)[["rank", "percentile", "T", "A", "R", "C"]].copy()
    target_table.insert(0, "expression_decile", deciles.reindex(TARGETS))
    target_table.insert(0, "evaluable", target_table["rank"].notna())
    output = root / "results" / "target-evaluation"
    output.mkdir(parents=True, exist_ok=True)
    target_table.to_csv(output / "known-target-ranks.csv", index_label="gene")
    pd.DataFrame({"mean_percentile": null}).to_csv(output / "matched-null.csv", index=False)
    summary = {
        "targets": TARGETS,
        "evaluable_targets": evaluable,
        "missing_targets": missing,
        "complete_five_target_endpoint": len(missing) == 0,
        "candidate_universe": int(len(rank)),
        "observed_mean_percentile_evaluable_only": observed,
        "null_mean": float(null.mean()),
        "null_interval_95": [float(np.quantile(null, 0.025)), float(np.quantile(null, 0.975))],
        "one_sided_empirical_p": p_value,
        "resamples": len(null),
        "seed": 271828,
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
