#!/usr/bin/env python3
"""Apply frozen Adaptive Y A1 to discovery interventions and factor screens."""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats

from core import load_homologene_one_to_one, percentile_interval, rank_score, read_gtf_gene_map
from ridge_y import score_ridge_rank_model
from run_discovery import (
    load_geo_expression,
    load_gse246954,
    matched_mptr_differences,
    paired_chemical_differences,
)


def load_programmes(root: Path) -> tuple[list[str], list[str], dict[str, list[str]]]:
    base = root / "results" / "discovery"
    identity = pd.read_csv(base / "identity-programme.csv")["gene"].tolist()
    pluripotency = pd.read_csv(base / "pluripotency-programme.csv")["gene"].tolist()
    damage_table = pd.read_csv(base / "damage-programmes.csv")
    damage = {go_id: group["gene"].tolist() for go_id, group in damage_table.groupby("go_id")}
    return identity, pluripotency, damage


def load_model(root: Path) -> tuple[list[str], pd.Series, pd.Series, float]:
    table = pd.read_csv(root / "results" / "adaptive-y-a1" / "model.csv")
    metadata = json.loads((root / "results" / "adaptive-y-a1" / "model.json").read_text())
    features = table["feature"].tolist()
    means = table.set_index("feature")["training_mean"]
    coefficients = table.set_index("feature")["coefficient"]
    return features, means, coefficients, float(metadata["intercept"])


def score_non_y(expression: pd.DataFrame, identity: list[str], pluripotency: list[str], damage: dict[str, list[str]]) -> pd.DataFrame:
    scores = pd.DataFrame(index=expression.columns)
    scores["I"] = rank_score(expression, identity)
    scores["P"] = rank_score(expression, pluripotency)
    d_columns: list[str] = []
    for go_id, genes in damage.items():
        column = f"D_{go_id.replace(':', '_')}"
        scores[column] = rank_score(expression, genes)
        d_columns.append(column)
    scores["D"] = scores[d_columns].mean(axis=1, skipna=False)
    return scores


def pseudobulk_h5ad(path: Path, grouping: list[str], ortho_map: pd.Series) -> tuple[pd.DataFrame, pd.Series]:
    """Return pseudobulk counts and the minimum NT-cell detection fraction per gene.

    The detection fraction implements the frozen identity-programme eligibility
    rule.  It is calculated separately in every available age/pool NT group and
    then reduced to the minimum, so a gene must be detected in at least 25% of
    cells in every control pool to be eligible.
    """
    data = ad.read_h5ad(path, backed="r")
    missing = [column for column in grouping if column not in data.obs]
    if missing:
        raise ValueError(f"Missing grouping fields in {path.name}: {missing}")
    matrix = data.layers["counts"]
    records: list[np.ndarray] = []
    names: list[str] = []
    nt_detection: list[np.ndarray] = []
    for key, indices in data.obs.groupby(grouping, observed=True).indices.items():
        key_tuple = key if isinstance(key, tuple) else (key,)
        indices = np.sort(indices)
        chunk = matrix[indices]
        if not sparse.issparse(chunk):
            chunk = sparse.csr_matrix(chunk)
        records.append(np.asarray(chunk.sum(axis=0)).ravel())
        names.append("|".join(map(str, key_tuple)))
        if str(key_tuple[-1]) == "NT":
            nt_detection.append(np.asarray((chunk > 0).mean(axis=0)).ravel())
    if not nt_detection:
        data.file.close()
        raise ValueError(f"No NT groups found in {path.name}")
    counts = pd.DataFrame(np.column_stack(records), index=data.var_names.astype(str), columns=names)
    detection = pd.Series(np.vstack(nt_detection).min(axis=0), index=data.var_names.astype(str))
    data.file.close()
    human = counts.index.to_series().map(ortho_map).fillna("")
    counts.index = human.to_numpy()
    counts = counts.loc[counts.index != ""].groupby(level=0).sum()
    detection.index = human.to_numpy()
    detection = detection.loc[detection.index != ""].groupby(level=0).max()
    return counts, detection


def log2_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    return np.log2(counts.divide(counts.sum(axis=0), axis=1) * 1_000_000 + 0.5)


def factor_deltas(scores: pd.DataFrame) -> pd.DataFrame:
    rows: list[pd.Series] = []
    for name in scores.index:
        parts = name.split("|")
        if parts[-1] in {"NT", "no_transgene_detected"}:
            continue
        control = "|".join([*parts[:-1], "NT"])
        if control not in scores.index:
            continue
        delta = scores.loc[name] - scores.loc[control]
        delta.name = name
        rows.append(delta)
    return pd.DataFrame(rows)


def main() -> None:
    root = Path(".").resolve()
    output = root / "results" / "adaptive-discovery-a1"
    output.mkdir(parents=True, exist_ok=True)
    features, means, coefficients, intercept = load_model(root)
    identity, pluripotency, damage = load_programmes(root)

    part1 = load_geo_expression(root / "inputs" / "discovery-data" / "GSE165177_Log2_RPM_Transient_reprogramming.txt.gz")
    part2 = load_geo_expression(root / "inputs" / "discovery-data" / "GSE165177_Log2_RPM_Transient_reprogramming_part2_170621.txt.gz")
    mptr = pd.concat([part1, part2], axis=1)
    mptr_scores = score_non_y(mptr, identity, pluripotency, damage)
    mptr_scores.insert(0, "Y_A1", score_ridge_rank_model(mptr, features, means, coefficients, intercept))
    mptr_delta, _ = matched_mptr_differences(mptr, mptr_scores, "transiently_reprogrammed")
    failed_delta, _ = matched_mptr_differences(mptr, mptr_scores, "failed_to_transiently_reprogram")
    donor = mptr_delta.assign(donor=mptr_delta.index.str.split("|").str[0]).groupby("donor").mean(numeric_only=True)
    mptr_interval = percentile_interval(donor, seed=1729)
    mptr_scores.to_csv(output / "gse165177-state-scores.csv", index_label="sample")
    mptr_delta.to_csv(output / "gse165177-successful-paired-deltas.csv", index_label="pair")
    failed_delta.to_csv(output / "gse165177-failed-paired-deltas.csv", index_label="pair")
    donor.to_csv(output / "gse165177-donor-deltas.csv", index_label="donor")
    mptr_interval.to_csv(output / "gse165177-donor-bootstrap.csv", index_label="axis")

    mouse_map = read_gtf_gene_map(root / "inputs" / "annotations" / "Mus_musculus.GRCm39.111.gtf.gz")
    orthology = load_homologene_one_to_one(root / "inputs" / "annotations" / "homologene_build68.data")
    chemical, _ = load_gse246954(root / "data" / "working" / "GSE246954_raw_genecounts.xls", mouse_map, orthology)
    chemical_scores = score_non_y(chemical, identity, pluripotency, damage)
    chemical_scores.insert(0, "Y_A1", score_ridge_rank_model(chemical, features, means, coefficients, intercept))
    chemical_delta, _ = paired_chemical_differences(chemical, chemical_scores)
    age_control = pd.DataFrame(
        {
            f"A{replicate}|B{replicate}": chemical_scores.loc[f"B{replicate}control"]
            - chemical_scores.loc[f"A{replicate}control"]
            for replicate in range(1, 5)
        }
    ).T
    chemical_scores.to_csv(output / "gse246954-state-scores.csv", index_label="sample")
    chemical_delta.to_csv(output / "gse246954-old-paired-deltas.csv", index_label="pair")
    age_control.to_csv(output / "gse246954-young-minus-old-controls.csv", index_label="pair")
    age_control_interval = percentile_interval(age_control, seed=1729)
    age_control_interval.to_csv(output / "gse246954-young-minus-old-bootstrap.csv", index_label="axis")
    intervals: list[pd.DataFrame] = []
    for treatment in ["2c", "7c"]:
        interval = percentile_interval(chemical_delta[chemical_delta.index.str.endswith(treatment)], seed=1729).reset_index(names="axis")
        interval.insert(0, "treatment", treatment)
        intervals.append(interval)
    chemical_interval = pd.concat(intervals, ignore_index=True)
    chemical_interval.to_csv(output / "gse246954-paired-bootstrap.csv", index=False)

    ortho_map = orthology.set_index("mouse_symbol")["human_symbol"]
    adipo_counts, adipo_detection = pseudobulk_h5ad(root / "data" / "working" / "GSE176206_adipo_screen.h5ad", ["age", "experiment", "combination_short"], ortho_map)
    msc_counts, msc_detection = pseudobulk_h5ad(root / "data" / "working" / "GSE176206_msc_screen.h5ad", ["age", "batch", "combination_short"], ortho_map)
    adipo, msc = log2_cpm(adipo_counts), log2_cpm(msc_counts)

    adipo_nt = adipo[[column for column in adipo if column.endswith("|NT")]].mean(axis=1)
    msc_nt = msc[[column for column in msc if column.endswith("|NT")]].mean(axis=1)
    common = adipo_nt.index.intersection(msc_nt.index)
    difference = adipo_nt.loc[common] - msc_nt.loc[common]
    adipo_eligible = difference.index[(difference >= 1) & (adipo_detection.reindex(difference.index).fillna(0) >= 0.25)]
    msc_eligible = difference.index[(difference <= -1) & (msc_detection.reindex(difference.index).fillna(0) >= 0.25)]
    adipo_identity = difference.loc[adipo_eligible].nlargest(100).index.tolist()
    msc_identity = (-difference.loc[msc_eligible]).nlargest(100).index.tolist()
    pd.DataFrame({"gene": adipo_identity}).to_csv(output / "gse176206-adipogenic-identity-programme.csv", index=False)
    pd.DataFrame({"gene": msc_identity}).to_csv(output / "gse176206-msc-identity-programme.csv", index=False)

    factor_outputs: list[pd.DataFrame] = []
    for cell_type, expression, identity_programme in [
        ("adipogenic", adipo, adipo_identity),
        ("MSC", msc, msc_identity),
    ]:
        scores = score_non_y(expression, identity_programme, pluripotency, damage)
        scores.insert(0, "Y_A1", score_ridge_rank_model(expression, features, means, coefficients, intercept))
        scores.insert(0, "cell_type", cell_type)
        scores.to_csv(output / f"gse176206-{cell_type.lower()}-state-scores.csv", index_label="pseudobulk")
        numeric = scores.drop(columns="cell_type")
        deltas = factor_deltas(numeric)
        deltas.insert(0, "cell_type", cell_type)
        deltas.to_csv(output / f"gse176206-{cell_type.lower()}-factor-deltas.csv", index_label="pseudobulk")
        factor_outputs.append(deltas)
    all_factors = pd.concat(factor_outputs)
    factor_correlation = float(stats.spearmanr(all_factors["Y_A1"], all_factors["I"], nan_policy="omit").statistic)

    separation_pass = bool(
        mptr_interval.loc["Y_A1", "ci_low"] > 0
        and mptr_interval.loc["I", "ci_low"] > -0.05
        and mptr_interval.loc["P", "ci_high"] < 0.05
    )
    guarded_pass = bool(separation_pass and mptr_interval.loc["D", "ci_high"] < 0)
    summary = {
        "adaptive_model": "Y_A1",
        "mptr_matched_pairs": len(mptr_delta),
        "mptr_failed_matched_pairs": len(failed_delta),
        "mptr_donors": len(donor),
        "mptr_separation_pass": separation_pass,
        "mptr_guarded_pass": guarded_pass,
        "gse176206_adipo_pseudobulks": adipo.shape[1],
        "gse176206_msc_pseudobulks": msc.shape[1],
        "gse176206_factor_delta_rows": len(all_factors),
        "gse176206_spearman_delta_y_delta_i": factor_correlation,
        "gse246954_young_minus_old_y_estimate": float(age_control_interval.loc["Y_A1", "estimate"]),
        "gse246954_young_minus_old_y_ci": [
            float(age_control_interval.loc["Y_A1", "ci_low"]),
            float(age_control_interval.loc["Y_A1", "ci_high"]),
        ],
        "gse246954_2c_y_estimate": float(chemical_interval.query("treatment == '2c' and axis == 'Y_A1'")["estimate"].iloc[0]),
        "gse246954_7c_y_estimate": float(chemical_interval.query("treatment == '7c' and axis == 'Y_A1'")["estimate"].iloc[0]),
        "validation_unsealed": False,
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
