#!/usr/bin/env python3
"""Repair B1 outer CV target standardisation without changing final models or contrasts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from core import assign_stratified_folds
from run_benchmark_b1 import METHODS, fit_method, predict_y, reference_data


EVALUABLE_METHODS = [method for method in METHODS if method != "elastic_net_rank"]


def standardise_training_age(age: pd.Series, study: pd.Series) -> pd.Series:
    result = pd.Series(index=age.index, dtype=float)
    for cohort, members in study.groupby(study).groups.items():
        values = age.loc[members].astype(float)
        result.loc[members] = (values - values.mean()) / values.std(ddof=1)
    return result


def folds_for_references(root: Path, study: pd.Series, metadata: pd.DataFrame) -> pd.Series:
    meta_a = metadata.loc[study.index[study == "GSE113957"]]
    meta_b = metadata.loc[study.index[study == "GSE226189"]].copy()
    meta_b["age"] = meta_b["age (years)"]
    return pd.concat(
        [assign_stratified_folds(meta_a, seed=1729), assign_stratified_folds(meta_b, seed=1729)]
    ).loc[study.index]


def run_fold(root: Path, fold: int) -> None:
    output = root / "results" / "benchmark-b1" / "cv-repair"
    output.mkdir(parents=True, exist_ok=True)
    x, _, study, _, age_frame, metadata = reference_data(root)
    age = age_frame["age"]
    folds = folds_for_references(root, study, metadata)
    train = folds.index[folds != fold]
    held = folds.index[folds == fold]
    y_train = standardise_training_age(age.loc[train], study.loc[train])
    rows = []
    for method in EVALUABLE_METHODS:
        model = fit_method(method, x.loc[train], y_train, study.loc[train])
        prediction = predict_y(model, x.loc[held])
        for sample in held:
            rows.append(
                {
                    "method": method,
                    "sample": sample,
                    "study": study.loc[sample],
                    "age": float(age.loc[sample]),
                    "fold": fold,
                    "Y": float(prediction.loc[sample]),
                }
            )
    pd.DataFrame(rows).to_csv(output / f"fold-{fold}.csv", index=False)


def combine(root: Path) -> None:
    base = root / "results" / "benchmark-b1"
    paths = [base / "cv-repair" / f"fold-{fold}.csv" for fold in range(5)]
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError(missing)
    predictions = pd.concat([pd.read_csv(path) for path in paths], ignore_index=True)
    predictions.to_csv(base / "outer-cv-predictions-r2.csv", index=False)
    correlations = []
    for (method, cohort), group in predictions.groupby(["method", "study"]):
        correlations.append(
            {"method": method, "scope": cohort, "spearman_y_age": stats.spearmanr(group["Y"], group["age"]).statistic}
        )
    for method, group in predictions.groupby("method"):
        correlations.append(
            {"method": method, "scope": "pooled", "spearman_y_age": stats.spearmanr(group["Y"], group["age"]).statistic}
        )
    correlation_table = pd.DataFrame(correlations)
    correlation_table.to_csv(base / "outer-cv-correlations-r2.csv", index=False)

    sign = pd.read_csv(base / "method-sign-agreement.csv", index_col=0).loc[EVALUABLE_METHODS, EVALUABLE_METHODS]
    upper = sign.to_numpy()[np.triu_indices(len(EVALUABLE_METHODS), k=1)]
    age_controls = pd.read_csv(base / "age-control-means.csv", index_col=0)[EVALUABLE_METHODS]
    fractions = pd.read_csv(base / "method-favourable-fractions.csv", index_col=0).loc[EVALUABLE_METHODS, "positive_family_fraction"]
    summary = json.loads((base / "summary.json").read_text())
    summary.update(
        {
            "revision": "r2_training_only_cv_and_constant_elastic_exclusion",
            "evaluable_methods": EVALUABLE_METHODS,
            "excluded_method": "elastic_net_rank",
            "excluded_method_reason": "Final alpha 10.0/l1_ratio 0.1 model has zero non-zero coefficients and constant intervention outputs.",
            "outer_cv_pooled_spearman": correlation_table[correlation_table["scope"] == "pooled"].set_index("method")["spearman_y_age"].to_dict(),
            "methods_orienting_both_age_controls": int(age_controls.gt(0).all(axis=0).sum()),
            "median_pairwise_method_sign_agreement": float(np.median(upper)),
            "transportability_gate_pass": bool(age_controls.gt(0).all(axis=0).sum() >= 4 and np.median(upper) >= 0.80),
            "method_favourable_family_fractions": fractions.to_dict(),
            "methods_at_or_above_80_percent_families": int((fractions >= 0.80).sum()),
            "universal_direction_gate_pass": int((fractions >= 0.80).sum()) >= 4,
        }
    )
    (base / "summary-r2.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fold", type=int, choices=range(5))
    parser.add_argument("--combine", action="store_true")
    args = parser.parse_args()
    root = Path(".").resolve()
    if args.fold is not None:
        run_fold(root, args.fold)
    elif args.combine:
        combine(root)
    else:
        parser.error("use --fold N or --combine")


if __name__ == "__main__":
    main()
