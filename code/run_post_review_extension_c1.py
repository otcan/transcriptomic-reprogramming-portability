#!/usr/bin/env python3
"""Run the locked post-review C1 sensitivity and contemporary-clock extension."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import rdata
from scipy import stats

from core import read_gtf_gene_map
from run_benchmark_b1 import (
    FAVOURABLE_FAMILIES,
    fit_method,
    prepare_expressions,
    reference_data,
    score_expression,
)


PASTA_FILES = {
    "v_genes_model.rda": "42289278a0a0d5171af23b108b086deea7b413d287559a9ec39cf7c76269cb5f",
    "beta_Pasta.rda": "b4a1891b7ed0737fbb784a89f6e00a9dfa07a15fa53359975f78f42fe92722e8",
    "cvfit_Pasta.rda": "158a9ec9d9e53786c85bfc8af3593bfa99ae4bc949c4ad3f2afb767c14befa56",
}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def pair_names(methods: list[str]) -> list[tuple[str, str]]:
    return list(itertools.combinations(methods, 2))


def pairwise_sign_agreement(table: pd.DataFrame, methods: list[str]) -> pd.DataFrame:
    rows = []
    for first, second in pair_names(methods):
        valid = table[[first, second]].dropna()
        agreement = ((valid[first] > 0) == (valid[second] > 0)).mean() if len(valid) else np.nan
        rows.append({"method_1": first, "method_2": second, "n": len(valid), "agreement": agreement})
    return pd.DataFrame(rows)


def family_pairwise_agreement(table: pd.DataFrame, methods: list[str]) -> pd.DataFrame:
    rows = []
    for family, group in table.groupby("family", sort=True):
        for row in pairwise_sign_agreement(group, methods).to_dict("records"):
            rows.append({"family": family, **row})
    return pd.DataFrame(rows)


def equal_family_statistic(family_pairs: pd.DataFrame) -> tuple[float, pd.DataFrame]:
    per_pair = (
        family_pairs.groupby(["method_1", "method_2"], as_index=False)
        .agg(agreement=("agreement", "mean"), families=("family", "nunique"))
    )
    return float(per_pair["agreement"].median()), per_pair


def cluster_bootstrap(
    family_pairs: pd.DataFrame, draws: int, seed: int
) -> tuple[pd.DataFrame, dict[str, float]]:
    matrix = family_pairs.pivot(index="family", columns=["method_1", "method_2"], values="agreement")
    rng = np.random.default_rng(seed)
    values = matrix.to_numpy(dtype=float)
    results = np.empty(draws, dtype=float)
    for draw in range(draws):
        selected = rng.integers(0, len(matrix), size=len(matrix))
        results[draw] = np.nanmedian(np.nanmean(values[selected], axis=0))
    frame = pd.DataFrame({"draw": np.arange(draws), "equal_family_median_agreement": results})
    summary = {
        "draws": draws,
        "seed": seed,
        "lower_95": float(np.quantile(results, 0.025)),
        "median": float(np.quantile(results, 0.5)),
        "upper_95": float(np.quantile(results, 0.975)),
    }
    return frame, summary


def leave_one_family_out(family_pairs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for held in sorted(family_pairs["family"].unique()):
        statistic, _ = equal_family_statistic(family_pairs[family_pairs["family"] != held])
        rows.append({"held_family": held, "equal_family_median_agreement": statistic})
    return pd.DataFrame(rows)


def agreement_threshold_table(
    contrast_statistic: float, family_statistic: float, thresholds: list[float]
) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "threshold": thresholds,
            "contrast_weighted_pass": [contrast_statistic >= value for value in thresholds],
            "equal_family_pass": [family_statistic >= value for value in thresholds],
        }
    )


def favourable_threshold_table(family_means: pd.DataFrame, methods: list[str]) -> pd.DataFrame:
    positive = family_means[methods].gt(0).sum(axis=0)
    rows = []
    for required in range(4, len(family_means) + 1):
        rows.append(
            {
                "required_positive_families": required,
                "threshold_fraction": required / len(family_means),
                "methods_meeting_threshold": int((positive >= required).sum()),
            }
        )
    return pd.DataFrame(rows)


def stratified_summary(
    family_means: pd.DataFrame, methods: list[str], strata: dict[str, dict[str, str]]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    family_rows = []
    summary_rows = []
    for dimension in ["species", "intervention", "context"]:
        values = sorted({mapping[dimension] for mapping in strata.values()})
        for value in values:
            families = sorted(family for family, mapping in strata.items() if mapping[dimension] == value)
            if len(families) < 2:
                continue
            subset = family_means.loc[families, methods]
            agreements = pairwise_sign_agreement(subset.reset_index(drop=True), methods)
            summary_rows.append(
                {
                    "dimension": dimension,
                    "stratum": value,
                    "family_count": len(families),
                    "median_pairwise_sign_agreement": float(agreements["agreement"].median()),
                    "minimum_pairwise_sign_agreement": float(agreements["agreement"].min()),
                    "maximum_pairwise_sign_agreement": float(agreements["agreement"].max()),
                }
            )
            for family in families:
                for method in methods:
                    family_rows.append(
                        {
                            "dimension": dimension,
                            "stratum": value,
                            "family": family,
                            "method": method,
                            "effect": family_means.loc[family, method],
                            "positive": bool(family_means.loc[family, method] > 0),
                        }
                    )
    return pd.DataFrame(summary_rows), pd.DataFrame(family_rows)


def spearman_vectors(first: pd.Series, second: pd.Series, minimum: int = 1000) -> tuple[float, int]:
    joined = pd.concat([first.rename("first"), second.rename("second")], axis=1).dropna()
    if len(joined) < minimum:
        return np.nan, len(joined)
    return float(stats.spearmanr(joined["first"], joined["second"]).statistic), len(joined)


def spearman_brown(value: float) -> float:
    if not np.isfinite(value) or value <= -1:
        return np.nan
    return float(2 * value / (1 + value))


def unique_balanced_splits(count: int, maximum: int, seed: int) -> list[tuple[tuple[int, ...], tuple[int, ...]]]:
    half = count // 2
    combinations = list(itertools.combinations(range(count), half))
    # Mirrored halves represent the same reliability comparison when sizes are equal.
    if count % 2 == 0:
        combinations = [combo for combo in combinations if 0 in combo]
    if len(combinations) > maximum:
        rng = np.random.default_rng(seed)
        chosen = np.sort(rng.choice(len(combinations), size=maximum, replace=False))
        combinations = [combinations[index] for index in chosen]
    all_indices = set(range(count))
    return [(combo, tuple(sorted(all_indices.difference(combo)))) for combo in combinations]


def reliability_analyses(
    expressions: dict[str, pd.DataFrame], contrasts: list, favourable: list[str], config: dict
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    contrast_vectors: dict[str, list[tuple[str, pd.Series]]] = {family: [] for family in favourable}
    contrast_lookup = {}
    for contrast in contrasts:
        expression = expressions[contrast.expression_key]
        delta = expression[contrast.treated].mean(axis=1) - expression[contrast.control].mean(axis=1)
        contrast_lookup[contrast.contrast_id] = (contrast, expression, delta)
        if contrast.family in contrast_vectors:
            contrast_vectors[contrast.family].append((contrast.contrast_id, delta))

    family_rows = []
    for family in favourable:
        vectors = contrast_vectors[family]
        if len(vectors) < 4:
            family_rows.append(
                {
                    "family": family,
                    "contrast_units": len(vectors),
                    "split_count": 0,
                    "median_spearman": np.nan,
                    "median_spearman_brown": np.nan,
                    "status": "INELIGIBLE_LT4_CONTRAST_UNITS",
                }
            )
            continue
        splits = unique_balanced_splits(
            len(vectors), config["family_split_draws"], config["family_split_seed"]
        )
        correlations = []
        for first_indices, second_indices in splits:
            first = pd.concat([vectors[index][1] for index in first_indices], axis=1).mean(axis=1)
            second = pd.concat([vectors[index][1] for index in second_indices], axis=1).mean(axis=1)
            correlation, genes = spearman_vectors(first, second)
            correlations.append((correlation, genes))
        raw = np.array([value[0] for value in correlations], dtype=float)
        family_rows.append(
            {
                "family": family,
                "contrast_units": len(vectors),
                "split_count": len(splits),
                "median_spearman": float(np.nanmedian(raw)),
                "median_spearman_brown": float(np.nanmedian([spearman_brown(value) for value in raw])),
                "minimum_common_genes": int(min(value[1] for value in correlations)),
                "status": "ESTIMATED",
            }
        )

    arm_rows = []
    draw_rows = []
    for contrast_id, (contrast, expression, _) in contrast_lookup.items():
        if contrast.family not in favourable or len(contrast.treated) < 2 or len(contrast.control) < 2:
            continue
        treated_splits = unique_balanced_splits(
            len(contrast.treated), config["arm_split_max_draws"], config["arm_split_seed"]
        )
        control_splits = unique_balanced_splits(
            len(contrast.control), config["arm_split_max_draws"], config["arm_split_seed"] + 1
        )
        pairs = list(itertools.product(treated_splits, control_splits))
        if len(pairs) > config["arm_split_max_draws"]:
            rng = np.random.default_rng(config["arm_split_seed"])
            selected = np.sort(rng.choice(len(pairs), config["arm_split_max_draws"], replace=False))
            pairs = [pairs[index] for index in selected]
        values = []
        for draw, ((ta, tb), (ca, cb)) in enumerate(pairs):
            delta_a = expression[[contrast.treated[index] for index in ta]].mean(axis=1) - expression[
                [contrast.control[index] for index in ca]
            ].mean(axis=1)
            delta_b = expression[[contrast.treated[index] for index in tb]].mean(axis=1) - expression[
                [contrast.control[index] for index in cb]
            ].mean(axis=1)
            correlation, genes = spearman_vectors(delta_a, delta_b)
            values.append(correlation)
            draw_rows.append(
                {
                    "contrast_id": contrast_id,
                    "family": contrast.family,
                    "draw": draw,
                    "spearman": correlation,
                    "spearman_brown": spearman_brown(correlation),
                    "common_genes": genes,
                }
            )
        arm_rows.append(
            {
                "contrast_id": contrast_id,
                "family": contrast.family,
                "treated_units": len(contrast.treated),
                "control_units": len(contrast.control),
                "split_count": len(values),
                "median_spearman": float(np.nanmedian(values)),
                "median_spearman_brown": float(np.nanmedian([spearman_brown(value) for value in values])),
            }
        )
    return pd.DataFrame(family_rows), pd.DataFrame(arm_rows), pd.DataFrame(draw_rows)


def feature_coverage(
    expressions: dict[str, pd.DataFrame], features: list[str], contrasts: list
) -> tuple[pd.DataFrame, dict[str, str]]:
    contrast_keys = {contrast.contrast_id: contrast.expression_key for contrast in contrasts}
    rows = []
    feature_set = set(features)
    for key, expression in sorted(expressions.items()):
        observed = len(feature_set.intersection(expression.index.astype(str)))
        rows.append(
            {
                "expression_key": key,
                "observed_features": observed,
                "total_features": len(features),
                "coverage": observed / len(features),
                "sample_count": expression.shape[1],
            }
        )
    return pd.DataFrame(rows), contrast_keys


def coverage_sensitivity(
    contrast_table: pd.DataFrame,
    contrast_keys: dict[str, str],
    coverage: pd.DataFrame,
    methods: list[str],
    thresholds: list[float],
) -> pd.DataFrame:
    coverage_lookup = coverage.set_index("expression_key")["coverage"].to_dict()
    work = contrast_table.copy()
    work["expression_key"] = work["contrast_id"].map(contrast_keys)
    work["coverage"] = work["expression_key"].map(coverage_lookup)
    rows = []
    for threshold in thresholds:
        subset = work[work["coverage"] >= threshold]
        contrast_pairs = pairwise_sign_agreement(subset, methods)
        family_pairs = family_pairwise_agreement(subset, methods)
        family_statistic, _ = equal_family_statistic(family_pairs)
        favourable = subset[subset["family"].isin(FAVOURABLE_FAMILIES)]
        rows.append(
            {
                "coverage_threshold": threshold,
                "expression_keys": subset["expression_key"].nunique(),
                "contrasts": len(subset),
                "families": subset["family"].nunique(),
                "favourable_families": favourable["family"].nunique(),
                "contrast_weighted_median_agreement": float(contrast_pairs["agreement"].median()),
                "equal_family_median_agreement": family_statistic,
            }
        )
    return pd.DataFrame(rows)


def score_models_on_contrasts(expressions: dict, contrasts: list, models: dict, features: list[str]) -> pd.DataFrame:
    lookup = {}
    for key, expression in expressions.items():
        for method, model in models.items():
            lookup[(key, method)] = score_expression(expression, features, model)
    rows = []
    for contrast in contrasts:
        row = {
            "family": contrast.family,
            "contrast_id": contrast.contrast_id,
            "dataset": contrast.dataset,
            "category": contrast.category,
            "favourable": contrast.favourable,
        }
        for method in models:
            values = lookup[(contrast.expression_key, method)]
            row[method] = float(values.loc[contrast.treated].mean() - values.loc[contrast.control].mean())
        rows.append(row)
    return pd.DataFrame(rows)


def common_feature_sensitivity(root: Path, expressions: dict, contrasts: list, methods: list[str], minimum: int) -> tuple[pd.DataFrame, dict]:
    x, y, study, original_features, _, _ = reference_data(root)
    common = set(original_features)
    for expression in expressions.values():
        common.intersection_update(expression.index.astype(str))
    features = sorted(common)
    summary = {
        "original_features": len(original_features),
        "common_features": len(features),
        "minimum_required": minimum,
        "status": "RUN" if len(features) >= minimum else "INFEASIBLE_BELOW_MINIMUM",
    }
    if len(features) < minimum:
        return pd.DataFrame(), summary
    models = {method: fit_method(method, x[features], y, study) for method in methods}
    table = score_models_on_contrasts(expressions, contrasts, models, features)
    contrast_stat = float(pairwise_sign_agreement(table, methods)["agreement"].median())
    family_stat, _ = equal_family_statistic(family_pairwise_agreement(table, methods))
    summary.update(
        {
            "contrast_weighted_median_agreement": contrast_stat,
            "equal_family_median_agreement": family_stat,
        }
    )
    return table, summary


def load_rdata(path: Path, name: str):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return rdata.conversion.convert(rdata.parser.parse_file(path))[name]


def load_pasta_model(pasta_dir: Path) -> dict:
    data = pasta_dir / "data"
    for filename, expected in PASTA_FILES.items():
        actual = digest(data / filename)
        if actual != expected:
            raise RuntimeError(f"Pasta model hash mismatch for {filename}: {actual}")
    genes = load_rdata(data / "v_genes_model.rda", "v_genes_model").astype(str)
    beta = float(load_rdata(data / "beta_Pasta.rda", "beta_Pasta")[0])
    cvfit = load_rdata(data / "cvfit_Pasta.rda", "cvfit_Pasta")
    fit = cvfit["glmnet.fit"]
    column = int(cvfit["index"].sel(dim_0="min").values[0]) - 1
    sparse = fit["beta"]
    coefficients = np.zeros(len(genes), dtype=float)
    start, end = int(sparse.p[column]), int(sparse.p[column + 1])
    coefficients[sparse.i[start:end]] = sparse.x[start:end]
    return {
        "genes": genes,
        "beta": beta,
        "intercept": float(fit["a0"][column]),
        "coefficients": coefficients,
        "lambda_min": float(fit["lambda"][column]),
    }


def pasta_input(expression: pd.DataFrame, genes: np.ndarray, gene_to_symbol: dict[str, str]) -> tuple[pd.DataFrame, int]:
    index = expression.index.astype(str).str.split(".").str[0]
    direct = pd.DataFrame(expression.to_numpy(), index=index, columns=expression.columns)
    direct = direct.groupby(level=0, sort=False).median()
    rows = []
    observed = 0
    for gene in genes:
        if gene in direct.index:
            rows.append(direct.loc[gene].rename(gene))
            observed += 1
            continue
        symbol = gene_to_symbol.get(gene)
        if symbol is not None and symbol in expression.index:
            rows.append(expression.loc[symbol].rename(gene))
            observed += 1
        else:
            rows.append(pd.Series(np.nan, index=expression.columns, name=gene))
    matrix = pd.DataFrame(rows)
    median = float(np.nanmedian(matrix.to_numpy(dtype=float)))
    matrix = matrix.fillna(median).rank(axis=0, method="average")
    return matrix, observed


def pasta_predict(matrix: pd.DataFrame, model: dict) -> pd.Series:
    age_shift = (model["intercept"] + matrix.T.to_numpy(dtype=float) @ model["coefficients"]) * model["beta"]
    return pd.Series(-age_shift, index=matrix.columns, name="pasta_youth_direction")


def validate_pasta_example(pasta_dir: Path, model: dict) -> dict:
    example = load_rdata(pasta_dir / "data" / "ES_GSE103938.rda", "ES_GSE103938")
    array = example.assayData["exprs"]
    expression = pd.DataFrame(
        array.values,
        index=array.coords["dim_0"].values.astype(str),
        columns=array.coords["dim_1"].values.astype(str),
    )
    matrix, observed = pasta_input(expression, model["genes"], {})
    predicted_age_shift = -pasta_predict(matrix, model)
    expected = {"Prolif_1": 13.5843466008, "OSKM_1": -8.3980706048, "OSKM1nM_1": -15.6617237438}
    maximum_error = max(abs(float(predicted_age_shift.loc[sample]) - value) for sample, value in expected.items())
    return {
        "status": "PASS" if maximum_error < 1e-8 else "FAIL",
        "observed_model_genes": observed,
        "maximum_absolute_error": maximum_error,
        "expected": expected,
        "observed": {sample: float(predicted_age_shift.loc[sample]) for sample in expected},
    }


def run_pasta(
    root: Path, pasta_dir: Path, expressions: dict, contrasts: list
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    model = load_pasta_model(pasta_dir)
    validation = validate_pasta_example(pasta_dir, model)
    if validation["status"] != "PASS":
        raise RuntimeError(f"Pasta implementation validation failed: {validation}")
    human_map = read_gtf_gene_map(root / "inputs" / "annotations" / "Homo_sapiens.GRCh38.111.gtf.gz")
    gene_to_symbol = human_map.drop_duplicates("gene_id").set_index("gene_id")["symbol"].to_dict()
    score_lookup = {}
    coverage_rows = []
    score_rows = []
    for key, expression in sorted(expressions.items()):
        matrix, observed = pasta_input(expression, model["genes"], gene_to_symbol)
        scores = pasta_predict(matrix, model)
        score_lookup[key] = scores
        coverage_rows.append(
            {
                "expression_key": key,
                "observed_features": observed,
                "total_features": len(model["genes"]),
                "coverage": observed / len(model["genes"]),
                "sample_count": expression.shape[1],
            }
        )
        for sample, value in scores.items():
            score_rows.append({"expression_key": key, "sample": sample, "pasta_youth_direction": value})
    contrast_rows = []
    for contrast in contrasts:
        scores = score_lookup[contrast.expression_key]
        effect = float(scores.loc[contrast.treated].mean() - scores.loc[contrast.control].mean())
        contrast_rows.append(
            {
                "family": contrast.family,
                "contrast_id": contrast.contrast_id,
                "dataset": contrast.dataset,
                "expression_key": contrast.expression_key,
                "category": contrast.category,
                "favourable": contrast.favourable,
                "pasta_youth_direction": effect,
            }
        )
    contrast_table = pd.DataFrame(contrast_rows)
    family = (
        contrast_table[contrast_table["family"].isin(FAVOURABLE_FAMILIES)]
        .groupby("family")["pasta_youth_direction"]
        .mean()
        .to_frame()
    )
    age = contrast_table[contrast_table["category"] == "age_control"].groupby("family")[
        "pasta_youth_direction"
    ].mean()
    summary = {
        "implementation_validation": validation,
        "lambda_min": model["lambda_min"],
        "model_genes": len(model["genes"]),
        "age_control_effects": age.to_dict(),
        "both_age_controls_positive": bool((age > 0).all()),
        "positive_favourable_families": int((family["pasta_youth_direction"] > 0).sum()),
        "favourable_family_count": len(family),
        "positive_favourable_fraction": float((family["pasta_youth_direction"] > 0).mean()),
    }
    return pd.DataFrame(score_rows), pd.DataFrame(coverage_rows), contrast_table, family, summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--pasta-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    output = root / "results" / "post-review-c1"
    output.mkdir(parents=True, exist_ok=True)
    config = json.loads((root / "protocol" / "post-review-extension-c1.json").read_text())
    methods = config["evaluable_methods"]

    contrast_table = pd.read_csv(root / "results" / "benchmark-b1" / "contrast-matrix.csv")
    contrast_pairs = pairwise_sign_agreement(contrast_table, methods)
    contrast_statistic = float(contrast_pairs["agreement"].median())
    family_pairs = family_pairwise_agreement(contrast_table, methods)
    family_statistic, equal_family_pairs = equal_family_statistic(family_pairs)
    bootstrap, bootstrap_summary = cluster_bootstrap(
        family_pairs, config["cluster_bootstrap_draws"], config["cluster_bootstrap_seed"]
    )
    loo = leave_one_family_out(family_pairs)
    thresholds = agreement_threshold_table(
        contrast_statistic, family_statistic, config["agreement_thresholds"]
    )
    favourable_means = (
        contrast_table[contrast_table["family"].isin(FAVOURABLE_FAMILIES)]
        .groupby("family")[methods]
        .mean()
    )
    favourable_thresholds = favourable_threshold_table(favourable_means, methods)
    strata_summary, strata_signs = stratified_summary(
        favourable_means, methods, config["favourable_family_strata"]
    )

    expressions, contrasts = prepare_expressions(root)
    _, _, _, features, _, _ = reference_data(root)
    coverage, contrast_keys = feature_coverage(expressions, features, contrasts)
    coverage_results = coverage_sensitivity(
        contrast_table, contrast_keys, coverage, methods, config["coverage_thresholds"]
    )
    family_reliability, arm_reliability, arm_draws = reliability_analyses(
        expressions, contrasts, FAVOURABLE_FAMILIES, config
    )
    common_table, common_summary = common_feature_sensitivity(
        root, expressions, contrasts, methods, config["common_feature_minimum"]
    )
    pasta_scores, pasta_coverage, pasta_contrasts, pasta_family, pasta_summary = run_pasta(
        root, args.pasta_dir.resolve(), expressions, contrasts
    )

    outputs = {
        "contrast-pair-agreement.csv": contrast_pairs,
        "family-pair-agreement.csv": family_pairs,
        "equal-family-pair-agreement.csv": equal_family_pairs,
        "family-cluster-bootstrap.csv": bootstrap,
        "leave-one-family-out-agreement.csv": loo,
        "agreement-threshold-sensitivity.csv": thresholds,
        "favourable-threshold-sensitivity.csv": favourable_thresholds,
        "stratified-summary.csv": strata_summary,
        "stratified-family-signs.csv": strata_signs,
        "feature-coverage.csv": coverage,
        "coverage-threshold-sensitivity.csv": coverage_results,
        "family-split-reliability.csv": family_reliability,
        "arm-split-reliability.csv": arm_reliability,
        "arm-split-reliability-draws.csv": arm_draws,
        "pasta-sample-scores.csv": pasta_scores,
        "pasta-feature-coverage.csv": pasta_coverage,
        "pasta-contrast-matrix.csv": pasta_contrasts,
        "pasta-favourable-family-means.csv": pasta_family.reset_index(),
    }
    if len(common_table):
        outputs["common-feature-contrast-matrix.csv"] = common_table
    for filename, frame in outputs.items():
        frame.to_csv(output / filename, index=False)

    summary = {
        "amendment": "C1",
        "post_review_outcome_aware": True,
        "methods": methods,
        "contrast_count": len(contrast_table),
        "benchmark_family_count": contrast_table["family"].nunique(),
        "contrast_weighted_median_pairwise_agreement": contrast_statistic,
        "equal_family_median_pairwise_agreement": family_statistic,
        "cluster_bootstrap": bootstrap_summary,
        "leave_one_family_out_minimum": float(loo["equal_family_median_agreement"].min()),
        "leave_one_family_out_maximum": float(loo["equal_family_median_agreement"].max()),
        "favourable_positive_counts": favourable_means.gt(0).sum(axis=0).astype(int).to_dict(),
        "family_split_reliability_estimable": int((family_reliability["status"] == "ESTIMATED").sum()),
        "arm_split_reliability_contrasts": len(arm_reliability),
        "feature_coverage_minimum": float(coverage["coverage"].min()),
        "feature_coverage_maximum": float(coverage["coverage"].max()),
        "common_feature_sensitivity": common_summary,
        "pasta": pasta_summary,
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
