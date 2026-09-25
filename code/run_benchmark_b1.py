#!/usr/bin/env python3
"""Run the locked adaptive B1 transportability benchmark."""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.linear_model import ElasticNetCV, RidgeCV

from core import assign_stratified_folds, load_homologene_one_to_one, load_metadata, read_gtf_gene_map
from run_discovery import (
    load_geo_expression,
    load_gse113957,
    load_gse226189,
    load_gse246954,
    parse_final_state,
)
from score_adaptive_discovery import log2_cpm, pseudobulk_h5ad
from score_validation import (
    filter_and_log_counts,
    map_human_ensembl_counts,
    map_mouse_counts,
    read_10x_pseudobulk,
    read_soft_records,
    read_table,
)


METHODS = [
    "ridge_rank_a1",
    "elastic_net_rank",
    "pca50_ridge",
    "meta_effect_projection",
    "gse113957_only_ridge",
    "gse226189_only_ridge",
]
FAVOURABLE_FAMILIES = [
    "GSE165177_MPTR",
    "GSE246954_chemical",
    "GSE176206_SOKM",
    "GSE297984_chemical",
    "GSE297233_OSK",
    "GSE297234_early_OSKM",
    "GSE304042_ARPE_OSK",
    "GSE304043_GSTA4",
]


@dataclass
class Contrast:
    family: str
    contrast_id: str
    dataset: str
    expression_key: str
    treated: list[str]
    control: list[str]
    category: str
    favourable: bool = False


def reference_data(root: Path) -> tuple[pd.DataFrame, pd.Series, pd.Series, list[str], pd.DataFrame, pd.DataFrame]:
    metadata = load_metadata(root / "receipts" / "sample-metadata.tsv")
    human_map = read_gtf_gene_map(root / "inputs" / "annotations" / "Homo_sapiens.GRCh38.111.gtf.gz")
    expression_a, _, meta_a = load_gse113957(
        root / "inputs" / "discovery-data" / "GSE113957_fpkm.txt.gz", metadata
    )
    expression_b, _, meta_b = load_gse226189(
        root / "inputs" / "discovery-data" / "GSE226189_RAW.tar", metadata, human_map
    )
    meta_b = meta_b.copy()
    meta_b["age"] = meta_b["age (years)"]
    protein_coding = set(human_map.loc[human_map["biotype"] == "protein_coding", "symbol"])
    features = sorted(set(expression_a.index).intersection(expression_b.index).intersection(protein_coding))
    x_a = expression_a.loc[features].rank(axis=0, pct=True).T
    x_b = expression_b.loc[features].rank(axis=0, pct=True).T
    x = pd.concat([x_a, x_b])
    y_a = (pd.to_numeric(meta_a["age"]) - pd.to_numeric(meta_a["age"]).mean()) / pd.to_numeric(meta_a["age"]).std(ddof=1)
    y_b = (pd.to_numeric(meta_b["age"]) - pd.to_numeric(meta_b["age"]).mean()) / pd.to_numeric(meta_b["age"]).std(ddof=1)
    y = pd.concat([y_a, y_b])
    study = pd.Series(["GSE113957"] * len(x_a) + ["GSE226189"] * len(x_b), index=x.index)
    age = pd.concat([pd.to_numeric(meta_a["age"]), pd.to_numeric(meta_b["age"])])
    return x, y, study, features, age.to_frame("age"), pd.concat([meta_a, meta_b])


def effect_weights(x: pd.DataFrame, y: pd.Series, study: pd.Series) -> pd.Series:
    effects = []
    for cohort in ["GSE113957", "GSE226189"]:
        subset = study.index[study == cohort]
        values = x.loc[subset].to_numpy(dtype=float)
        target = y.loc[subset].to_numpy(dtype=float)
        values -= values.mean(axis=0)
        target -= target.mean()
        denominator = np.sqrt((values**2).sum(axis=0) * (target**2).sum())
        correlation = np.divide(values.T @ target, denominator, out=np.zeros(values.shape[1]), where=denominator > 0)
        correlation = np.clip(correlation, -0.999999, 0.999999)
        z_value = correlation * np.sqrt((len(subset) - 2) / (1 - correlation**2))
        effects.append(z_value)
    a, b = effects
    weight = (a + b) / 2
    weight[np.sign(a) != np.sign(b)] = 0
    return pd.Series(weight, index=x.columns)


def fit_method(method: str, x: pd.DataFrame, y: pd.Series, study: pd.Series) -> dict[str, object]:
    means = x.mean(axis=0)
    alphas = np.logspace(-2, 6, 17)
    if method in {"ridge_rank_a1", "gse113957_only_ridge", "gse226189_only_ridge"}:
        train = x.index
        if method == "gse113957_only_ridge":
            train = study.index[study == "GSE113957"]
        elif method == "gse226189_only_ridge":
            train = study.index[study == "GSE226189"]
        model = RidgeCV(alphas=alphas, fit_intercept=True).fit(x.loc[train], y.loc[train])
        return {
            "kind": "linear",
            "means": means,
            "coefficients": pd.Series(model.coef_, index=x.columns),
            "intercept": float(model.intercept_),
            "hyperparameters": {"alpha": float(model.alpha_)},
        }
    if method == "elastic_net_rank":
        model = ElasticNetCV(
            l1_ratio=[0.1, 0.5, 0.9],
            alphas=np.logspace(-4, 1, 25),
            cv=5,
            max_iter=50_000,
            random_state=1729,
            n_jobs=1,
        ).fit(x, y)
        return {
            "kind": "linear",
            "means": means,
            "coefficients": pd.Series(model.coef_, index=x.columns),
            "intercept": float(model.intercept_),
            "hyperparameters": {"alpha": float(model.alpha_), "l1_ratio": float(model.l1_ratio_)},
        }
    if method == "pca50_ridge":
        n_components = min(50, len(x) - 1)
        pca = PCA(n_components=n_components, random_state=1729).fit(x)
        transformed = pca.transform(x)
        model = RidgeCV(alphas=alphas, fit_intercept=True).fit(transformed, y)
        return {
            "kind": "pca",
            "means": means,
            "pca": pca,
            "ridge": model,
            "hyperparameters": {"alpha": float(model.alpha_), "n_components": n_components},
        }
    if method == "meta_effect_projection":
        weights = effect_weights(x, y, study)
        return {
            "kind": "projection",
            "means": means,
            "weights": weights,
            "hyperparameters": {"nonzero_weights": int((weights != 0).sum())},
        }
    raise ValueError(method)


def predict_y(model: dict[str, object], x: pd.DataFrame) -> pd.Series:
    kind = model["kind"]
    if kind == "linear":
        older = float(model["intercept"]) + x.dot(model["coefficients"])
    elif kind == "pca":
        older = pd.Series(model["ridge"].predict(model["pca"].transform(x)), index=x.index)
    elif kind == "projection":
        weights = model["weights"]
        older = (x - model["means"]).dot(weights) / weights.abs().sum()
    else:
        raise ValueError(str(kind))
    return -pd.Series(older, index=x.index)


def score_expression(expression: pd.DataFrame, features: list[str], model: dict[str, object]) -> pd.Series:
    observed = [gene for gene in features if gene in expression.index]
    if len(observed) / len(features) < 0.60:
        return pd.Series(np.nan, index=expression.columns, dtype=float)
    ranks = expression.rank(axis=0, pct=True)
    x = pd.DataFrame(index=expression.columns, columns=features, dtype=float)
    x.loc[:, observed] = ranks.loc[observed].T
    x = x.fillna(model["means"])
    return predict_y(model, x)


def fit_and_cross_validate(root: Path, output: Path) -> tuple[dict[str, dict[str, object]], list[str]]:
    x, y, study, features, age, metadata = reference_data(root)
    meta_a = metadata.loc[study.index[study == "GSE113957"]]
    meta_b = metadata.loc[study.index[study == "GSE226189"]].copy()
    meta_b["age"] = meta_b["age (years)"]
    fold_a = assign_stratified_folds(meta_a, seed=1729)
    fold_b = assign_stratified_folds(meta_b, seed=1729)
    folds = pd.concat([fold_a, fold_b]).loc[x.index]
    prediction_rows = []
    final_models: dict[str, dict[str, object]] = {}
    hyperparameters = []
    for method in METHODS:
        for fold in range(5):
            train = folds.index[folds != fold]
            held = folds.index[folds == fold]
            model = fit_method(method, x.loc[train], y.loc[train], study.loc[train])
            prediction = predict_y(model, x.loc[held])
            for sample in held:
                prediction_rows.append(
                    {
                        "method": method,
                        "sample": sample,
                        "study": study.loc[sample],
                        "age": float(age.loc[sample, "age"]),
                        "fold": fold,
                        "Y": float(prediction.loc[sample]),
                    }
                )
        final = fit_method(method, x, y, study)
        final_models[method] = final
        hyperparameters.append({"method": method, **final["hyperparameters"]})
    predictions = pd.DataFrame(prediction_rows)
    predictions.to_csv(output / "outer-cv-predictions.csv", index=False)
    correlations = []
    for (method, cohort), group in predictions.groupby(["method", "study"]):
        correlations.append(
            {"method": method, "scope": cohort, "spearman_y_age": stats.spearmanr(group["Y"], group["age"]).statistic}
        )
    for method, group in predictions.groupby("method"):
        correlations.append(
            {"method": method, "scope": "pooled", "spearman_y_age": stats.spearmanr(group["Y"], group["age"]).statistic}
        )
    pd.DataFrame(correlations).to_csv(output / "outer-cv-correlations.csv", index=False)
    pd.DataFrame(hyperparameters).to_csv(output / "final-hyperparameters.csv", index=False)
    joblib.dump(final_models, output / "final-models.joblib", compress=3)
    return final_models, features


def prepare_expressions(root: Path) -> tuple[dict[str, pd.DataFrame], list[Contrast]]:
    discovery = root / "inputs" / "discovery-data"
    validation = root / "inputs" / "validation-data"
    mouse_map = read_gtf_gene_map(root / "inputs" / "annotations" / "Mus_musculus.GRCm39.111.gtf.gz")
    human_map = read_gtf_gene_map(root / "inputs" / "annotations" / "Homo_sapiens.GRCh38.111.gtf.gz")
    orthology = load_homologene_one_to_one(root / "inputs" / "annotations" / "homologene_build68.data")
    ortho_map = orthology.set_index("mouse_symbol")["human_symbol"]
    expressions: dict[str, pd.DataFrame] = {}
    contrasts: list[Contrast] = []

    part1 = load_geo_expression(discovery / "GSE165177_Log2_RPM_Transient_reprogramming.txt.gz")
    part2 = load_geo_expression(discovery / "GSE165177_Log2_RPM_Transient_reprogramming_part2_170621.txt.gz")
    expressions["gse165177"] = pd.concat([part1, part2], axis=1)
    for state, family, favourable, category in [
        ("transiently_reprogrammed", "GSE165177_MPTR", True, "favourable"),
        ("failed_to_transiently_reprogram", "GSE165177_failed", False, "negative"),
    ]:
        for sample in expressions["gse165177"].columns:
            parsed = parse_final_state(sample, state)
            if parsed is None:
                continue
            donor, day, experiment = parsed
            control = f"{donor}_negative_control_{day}days_{experiment}"
            if control in expressions["gse165177"].columns:
                contrasts.append(Contrast(family, f"{family}|{donor}|{day}|{experiment}", "GSE165177", "gse165177", [sample], [control], category, favourable))

    chemical, _ = load_gse246954(root / "data" / "working" / "GSE246954_raw_genecounts.xls", mouse_map, orthology)
    expressions["gse246954"] = chemical
    for replicate in range(1, 5):
        for treatment in ["2c", "7c"]:
            contrasts.append(Contrast("GSE246954_chemical", f"GSE246954|A{replicate}|{treatment}", "GSE246954", "gse246954", [f"A{replicate}_{treatment}"], [f"A{replicate}control"], "favourable", True))
        contrasts.append(Contrast("GSE246954_age", f"GSE246954|age|{replicate}", "GSE246954", "gse246954", [f"B{replicate}control"], [f"A{replicate}control"], "age_control"))

    adipo_counts, _ = pseudobulk_h5ad(root / "data" / "working" / "GSE176206_adipo_screen.h5ad", ["age", "experiment", "combination_short"], ortho_map)
    msc_counts, _ = pseudobulk_h5ad(root / "data" / "working" / "GSE176206_msc_screen.h5ad", ["age", "batch", "combination_short"], ortho_map)
    expressions["gse176206_adipo"] = log2_cpm(adipo_counts)
    expressions["gse176206_msc"] = log2_cpm(msc_counts)
    for key in ["gse176206_adipo", "gse176206_msc"]:
        for sample in expressions[key].columns:
            combo = sample.split("|")[-1]
            if combo == "NT":
                continue
            control = "|".join([*sample.split("|")[:-1], "NT"])
            favourable = combo == "SOKM"
            family = "GSE176206_SOKM" if favourable else "GSE176206_other_factors"
            contrasts.append(Contrast(family, f"{key}|{sample}", "GSE176206", key, [sample], [control], "favourable" if favourable else "exploratory", favourable))

    human_chemical = read_table(validation / "GSE297984_logCPM_2c_7c.csv.gz")
    expressions["gse297984"] = human_chemical
    titles = {re.match(r"(RNA_LS_\d+)", r["title"]).group(1): r["title"] for r in read_soft_records(validation / "GSE297984_family.soft.gz")}
    groups: dict[tuple[int, int, str], list[str]] = {}
    for sample, title in titles.items():
        match = re.match(r"RNA_LS_\d+_D(6|14)_(\d+)y_(DMSO|2c|7c)", title)
        if match:
            day, age, treatment = match.groups()
            if int(age) in [56, 83]:
                groups.setdefault((int(age), int(day), treatment), []).append(sample)
    for age in [56, 83]:
        for day in [6, 14]:
            for treatment in ["2c", "7c"]:
                contrasts.append(Contrast("GSE297984_chemical", f"GSE297984|{age}|{day}|{treatment}", "GSE297984", "gse297984", groups[(age, day, treatment)], groups[(age, day, "DMSO")], "favourable", True))

    osk_counts = read_table(validation / "GSE297233_raw_counts_matrix.csv.gz", index_col="GeneId")
    expressions["gse297233"] = filter_and_log_counts(map_human_ensembl_counts(osk_counts, human_map))
    for treatment in ["OSK", "O4YRSK"]:
        contrasts.append(Contrast("GSE297233_OSK", f"GSE297233|{treatment}", "GSE297233", "gse297233", [f"{treatment}_D4_1", f"{treatment}_D4_2"], [f"{treatment}_D0_1", f"{treatment}_D0_2"], "favourable", True))

    trajectory_series = {}
    for path in sorted(validation.glob("GSM89865*_filtered_feature_bc_matrix.h5")):
        match = re.match(r"GSM\d+_(GM\d+)_D(\d+)_", path.name)
        line, day = match.groups()
        trajectory_series[f"{line}|{day}"] = read_10x_pseudobulk(path)
    expressions["gse297234"] = filter_and_log_counts(pd.DataFrame(trajectory_series).fillna(0))
    for day in [3, 7, 10]:
        favourable = day in [3, 7]
        contrasts.append(Contrast("GSE297234_early_OSKM" if favourable else "GSE297234_day10", f"GSE297234|old|{day}", "GSE297234", "gse297234", [f"GM00731|{day}"], ["GM00731|0"], "favourable" if favourable else "trajectory", favourable))
    contrasts.append(Contrast("GSE297234_age", "GSE297234|age|day0", "GSE297234", "gse297234", ["GM23815|0"], ["GM00731|0"], "age_control"))

    adverse_raw = read_table(validation / "GSE300625_raw_genecounts.csv.gz")
    expressions["gse300625"] = filter_and_log_counts(map_mouse_counts(adverse_raw, mouse_map, orthology))
    adverse_meta = {}
    for record in read_soft_records(validation / "GSE300625_family.soft.gz"):
        library = record["description"].removeprefix("Library name: ")
        title = record["title"].lower()
        adverse_meta[library] = ("liver" if "liver" in title else "kidney", "7c" if title.startswith("7c") else "vehicle")
    for tissue in ["liver", "kidney"]:
        treated = [sample for sample, value in adverse_meta.items() if value == (tissue, "7c")]
        control = [sample for sample, value in adverse_meta.items() if value == (tissue, "vehicle")]
        contrasts.append(Contrast("GSE300625_adverse", f"GSE300625|{tissue}", "GSE300625", "gse300625", treated, control, "adverse"))

    expressions["gse304042"] = read_table(validation / "GSE304042_5_ARPE_single_triple_OSK_30-1011800743.csv.gz", index_col="Feature.ID")
    contrasts.append(Contrast("GSE304042_ARPE_OSK", "GSE304042|OSK", "GSE304042", "gse304042", ["OSK-mC-a", "OSK-mC-b"], ["GFP-a", "GFP-b", "GFP-c"], "favourable", True))
    gsta_raw = read_table(validation / "GSE304043_6_mRPE_Old-GFP_Old-GSTA4_merged_featureCounts.csv.gz", index_col="Gene name")
    expressions["gse304043"] = filter_and_log_counts(map_mouse_counts(gsta_raw, mouse_map, orthology))
    contrasts.append(Contrast("GSE304043_GSTA4", "GSE304043|GSTA4", "GSE304043", "gse304043", [c for c in expressions["gse304043"].columns if "GSTA4" in c], [c for c in expressions["gse304043"].columns if "GFP" in c], "favourable", True))
    return expressions, contrasts


def evaluate(root: Path, output: Path, models: dict[str, dict[str, object]], features: list[str]) -> dict[str, object]:
    expressions, contrasts = prepare_expressions(root)
    score_records = []
    score_lookup: dict[tuple[str, str], pd.Series] = {}
    for key, expression in expressions.items():
        for method, model in models.items():
            values = score_expression(expression, features, model)
            score_lookup[(key, method)] = values
            for sample, value in values.items():
                score_records.append({"expression_key": key, "sample": sample, "method": method, "Y": value})
    pd.DataFrame(score_records).to_csv(output / "all-sample-scores.csv", index=False)

    contrast_rows = []
    gene_deltas: dict[str, list[pd.Series]] = {}
    for contrast in contrasts:
        expression = expressions[contrast.expression_key]
        gene_delta = expression[contrast.treated].mean(axis=1) - expression[contrast.control].mean(axis=1)
        gene_deltas.setdefault(contrast.family, []).append(gene_delta)
        row = {
            "family": contrast.family,
            "contrast_id": contrast.contrast_id,
            "dataset": contrast.dataset,
            "category": contrast.category,
            "favourable": contrast.favourable,
        }
        for method in METHODS:
            values = score_lookup[(contrast.expression_key, method)]
            row[method] = float(values.loc[contrast.treated].mean() - values.loc[contrast.control].mean())
        contrast_rows.append(row)
    contrast_table = pd.DataFrame(contrast_rows)
    contrast_table.to_csv(output / "contrast-matrix.csv", index=False)

    method_values = contrast_table[METHODS]
    sign_agreement = pd.DataFrame(index=METHODS, columns=METHODS, dtype=float)
    correlations = pd.DataFrame(index=METHODS, columns=METHODS, dtype=float)
    for first in METHODS:
        for second in METHODS:
            sign_agreement.loc[first, second] = ((method_values[first] > 0) == (method_values[second] > 0)).mean()
            correlations.loc[first, second] = stats.spearmanr(method_values[first], method_values[second]).statistic
    sign_agreement.to_csv(output / "method-sign-agreement.csv", index_label="method")
    correlations.to_csv(output / "method-spearman.csv", index_label="method")

    family = contrast_table[contrast_table["family"].isin(FAVOURABLE_FAMILIES)].groupby("family")[METHODS].mean()
    family.to_csv(output / "favourable-family-means.csv", index_label="family")
    family_fraction = family.gt(0).mean(axis=0)
    family_fraction.rename("positive_family_fraction").to_csv(output / "method-favourable-fractions.csv", header=True)
    age_controls = contrast_table[contrast_table["category"] == "age_control"].groupby("family")[METHODS].mean()
    age_controls.to_csv(output / "age-control-means.csv", index_label="family")

    vector_series = {}
    for family_name in FAVOURABLE_FAMILIES:
        vectors = gene_deltas[family_name]
        family_vector = pd.concat(vectors, axis=1).mean(axis=1)
        vector_series[family_name] = (family_vector - family_vector.mean()) / family_vector.std(ddof=1)
    vectors = pd.DataFrame(vector_series)
    gene_correlations = vectors.corr(method="spearman", min_periods=1000)
    gene_correlations.to_csv(output / "family-gene-effect-spearman.csv", index_label="family")
    complete_vectors = vectors.dropna()
    loo = []
    for held in FAVOURABLE_FAMILIES:
        consensus = complete_vectors.drop(columns=held).mean(axis=1)
        loo.append({"held_family": held, "spearman": stats.spearmanr(consensus, complete_vectors[held]).statistic, "common_genes": len(complete_vectors)})
    pd.DataFrame(loo).to_csv(output / "leave-one-family-out-gene-consensus.csv", index=False)

    upper = sign_agreement.to_numpy()[np.triu_indices(len(METHODS), k=1)]
    cv = pd.read_csv(output / "outer-cv-correlations.csv")
    age_control_correct = age_controls.gt(0).all(axis=0)
    methods_with_both_age_controls = int(age_control_correct.sum())
    methods_at_80 = int((family_fraction >= 0.80).sum())
    summary = {
        "adaptive_post_primary_failure": True,
        "methods": METHODS,
        "contrast_count": int(len(contrast_table)),
        "favourable_family_count": len(FAVOURABLE_FAMILIES),
        "outer_cv_pooled_spearman": cv[cv["scope"] == "pooled"].set_index("method")["spearman_y_age"].to_dict(),
        "methods_orienting_both_age_controls": methods_with_both_age_controls,
        "median_pairwise_method_sign_agreement": float(np.median(upper)),
        "transportability_gate_pass": bool(methods_with_both_age_controls >= 4 and np.median(upper) >= 0.80),
        "method_favourable_family_fractions": family_fraction.to_dict(),
        "methods_at_or_above_80_percent_families": methods_at_80,
        "universal_direction_gate_pass": methods_at_80 >= 4,
        "median_pairwise_family_gene_effect_spearman": float(np.nanmedian(gene_correlations.to_numpy()[np.triu_indices(len(FAVOURABLE_FAMILIES), k=1)])),
        "median_leave_one_family_out_gene_consensus_spearman": float(pd.DataFrame(loo)["spearman"].median()),
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    root = Path(".").resolve()
    output = root / "results" / "benchmark-b1"
    output.mkdir(parents=True, exist_ok=True)
    models, features = fit_and_cross_validate(root, output)
    summary = evaluate(root, output, models, features)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
