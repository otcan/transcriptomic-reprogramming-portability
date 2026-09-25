#!/usr/bin/env python3
"""Fit and export Adaptive Amendment A1 using only the two age references."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from core import assign_stratified_folds, load_metadata, read_gtf_gene_map
from ridge_y import fit_ridge_rank_model, save_model
from run_discovery import load_gse113957, load_gse226189


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    root = Path(".").resolve()
    output = root / "results" / "adaptive-y-a1"
    output.mkdir(parents=True, exist_ok=True)
    metadata = load_metadata(root / "receipts" / "sample-metadata.tsv")
    human_map = read_gtf_gene_map(root / "inputs" / "annotations" / "Homo_sapiens.GRCh38.111.gtf.gz")
    expression_a, _, meta_a = load_gse113957(root / "inputs" / "discovery-data" / "GSE113957_fpkm.txt.gz", metadata)
    expression_b, _, meta_b = load_gse226189(root / "inputs" / "discovery-data" / "GSE226189_RAW.tar", metadata, human_map)
    meta_b = meta_b.copy()
    meta_b["age"] = meta_b["age (years)"]
    protein_coding = set(human_map.loc[human_map["biotype"] == "protein_coding", "symbol"])
    features = sorted(set(expression_a.index).intersection(expression_b.index).intersection(protein_coding))
    alphas = np.logspace(-2, 6, 17)
    fold_a = assign_stratified_folds(meta_a, seed=1729)
    fold_b = assign_stratified_folds(meta_b, seed=1729)
    predictions: list[pd.DataFrame] = []
    fold_alphas: list[dict[str, float | int]] = []
    for fold in range(5):
        train_a, train_b = fold_a != fold, fold_b != fold
        model, _ = fit_ridge_rank_model(
            expression_a.loc[:, train_a],
            expression_b.loc[:, train_b],
            pd.to_numeric(meta_a.loc[train_a, "age"]),
            pd.to_numeric(meta_b.loc[train_b, "age"]),
            features,
            alphas,
        )
        fold_alphas.append({"fold": fold, "alpha": float(model.alpha_)})
        for study, expression, meta, folds in [
            ("GSE113957", expression_a, meta_a, fold_a),
            ("GSE226189", expression_b, meta_b, fold_b),
        ]:
            held = folds[folds == fold].index
            x = expression.loc[features, held].rank(axis=0, pct=True).T
            older_prediction = model.predict(x)
            predictions.append(
                pd.DataFrame(
                    {
                        "sample": held,
                        "study": study,
                        "age": pd.to_numeric(meta.loc[held, "age"]).to_numpy(),
                        "Y_A1": -older_prediction,
                        "fold": fold,
                    }
                )
            )
    prediction_table = pd.concat(predictions, ignore_index=True)
    prediction_table.to_csv(output / "outer-predictions.csv", index=False)
    pd.DataFrame(fold_alphas).to_csv(output / "outer-alphas.csv", index=False)
    correlations = {
        "pooled": float(stats.spearmanr(prediction_table["Y_A1"], prediction_table["age"]).statistic)
    }
    for study, group in prediction_table.groupby("study"):
        correlations[study] = float(stats.spearmanr(group["Y_A1"], group["age"]).statistic)

    final_model, feature_means = fit_ridge_rank_model(
        expression_a,
        expression_b,
        pd.to_numeric(meta_a["age"]),
        pd.to_numeric(meta_b["age"]),
        features,
        alphas,
    )
    save_model(
        output / "model.csv",
        features,
        feature_means,
        final_model.coef_,
        final_model.intercept_,
        final_model.alpha_,
    )
    summary = {
        "adaptive": True,
        "protocol_v1_y_status": "FAILED_ZERO_GENES",
        "n_features": len(features),
        "final_alpha": float(final_model.alpha_),
        "outer_spearman_y_vs_age": correlations,
        "inputs": {
            "discovery_manifest_sha256": sha256(root / "inputs" / "discovery-manifest.tsv"),
            "adaptive_config_sha256": sha256(root / "protocol" / "adaptive-config-a1.json"),
        },
        "validation_unsealed": False,
        "intervention_state_outcomes_inspected_before_lock": False,
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

