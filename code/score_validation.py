#!/usr/bin/env python3
"""Apply the frozen Paper 2 model once to post-2024 temporal validation data."""

from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from core import collapse_expression, load_homologene_one_to_one, log2_cpm, read_gtf_gene_map
from ridge_y import score_ridge_rank_model
from score_adaptive_discovery import load_model, load_programmes, score_non_y


def read_soft_records(path: Path) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    records.append(current)
                current = {"accession": line.split(" = ", 1)[1]}
            elif current is not None and line.startswith("!Sample_"):
                key, _, value = line.partition(" = ")
                field = key.removeprefix("!Sample_")
                if field in {"title", "description"}:
                    current[field] = value
    if current is not None:
        records.append(current)
    return records


def read_table(path: Path, index_col: int | str = 0) -> pd.DataFrame:
    table = pd.read_csv(path, index_col=index_col)
    table.index = table.index.astype(str).str.removeprefix("\ufeff")
    return table.apply(pd.to_numeric, errors="raise")


def filter_and_log_counts(counts: pd.DataFrame) -> pd.DataFrame:
    cpm = counts.divide(counts.sum(axis=0).replace(0, np.nan), axis=1) * 1_000_000
    keep = (cpm >= 1).mean(axis=1) >= 0.20
    return log2_cpm(counts.loc[keep])


def map_human_ensembl_counts(counts: pd.DataFrame, gene_map: pd.DataFrame) -> pd.DataFrame:
    symbols = counts.index.to_series().str.split(".").str[0].map(gene_map.set_index("gene_id")["symbol"]).fillna("")
    return collapse_expression(counts, symbols, "sum")


def map_mouse_counts(counts: pd.DataFrame, gene_map: pd.DataFrame, orthology: pd.DataFrame) -> pd.DataFrame:
    identifiers = counts.index.to_series().str.split(".").str[0]
    ensembl_lookup = gene_map.set_index("gene_id")["symbol"]
    mouse_symbols = identifiers.map(ensembl_lookup)
    mouse_symbols = mouse_symbols.where(mouse_symbols.notna(), counts.index.to_series())
    mouse = collapse_expression(counts, mouse_symbols.fillna(""), "sum")
    human_symbols = mouse.index.to_series().map(orthology.set_index("mouse_symbol")["human_symbol"]).fillna("")
    return collapse_expression(mouse, human_symbols, "sum")


def read_10x_pseudobulk(path: Path) -> pd.Series:
    with h5py.File(path, "r") as handle:
        matrix = handle["matrix"]
        shape = matrix["shape"][:]
        genes = np.char.decode(matrix["features"]["name"][:], "utf-8")
        totals = np.bincount(
            matrix["indices"][:], weights=matrix["data"][:].astype(float), minlength=int(shape[0])
        )
    return pd.Series(totals, index=genes).groupby(level=0).sum()


def independent_bootstrap_delta(
    treated: pd.DataFrame,
    control: pd.DataFrame,
    n_resamples: int = 10_000,
    seed: int = 1729,
) -> pd.DataFrame:
    treated_array = treated.to_numpy(dtype=float)
    control_array = control.to_numpy(dtype=float)
    rng = np.random.default_rng(seed)
    draws = np.empty((n_resamples, treated.shape[1]), dtype=float)
    for idx in range(n_resamples):
        treated_draw = rng.integers(0, len(treated_array), len(treated_array))
        control_draw = rng.integers(0, len(control_array), len(control_array))
        draws[idx] = treated_array[treated_draw].mean(axis=0) - control_array[control_draw].mean(axis=0)
    return pd.DataFrame(
        {
            "estimate": treated_array.mean(axis=0) - control_array.mean(axis=0),
            "ci_low": np.percentile(draws, 2.5, axis=0),
            "ci_high": np.percentile(draws, 97.5, axis=0),
        },
        index=treated.columns,
    )


def score(expression: pd.DataFrame, model: tuple, programmes: tuple, include_identity: bool = True) -> pd.DataFrame:
    features, means, coefficients, intercept = model
    identity, pluripotency, damage = programmes
    scores = score_non_y(expression, identity, pluripotency, damage)
    scores.insert(0, "Y_A1", score_ridge_rank_model(expression, features, means, coefficients, intercept))
    if not include_identity:
        scores = scores.drop(columns="I")
    return scores


def main() -> None:
    root = Path(".").resolve()
    inputs = root / "inputs" / "validation-data"
    output = root / "results" / "validation"
    output.mkdir(parents=True, exist_ok=True)
    model = load_model(root)
    programmes = load_programmes(root)
    human_map = read_gtf_gene_map(root / "inputs" / "annotations" / "Homo_sapiens.GRCh38.111.gtf.gz")
    mouse_map = read_gtf_gene_map(root / "inputs" / "annotations" / "Mus_musculus.GRCm39.111.gtf.gz")
    orthology = load_homologene_one_to_one(root / "inputs" / "annotations" / "homologene_build68.data")

    # Positive human chemical challenge: nested cultures are collapsed within donor/day/treatment.
    chemical = read_table(inputs / "GSE297984_logCPM_2c_7c.csv.gz")
    chemical_scores = score(chemical, model, programmes)
    chemical_records = read_soft_records(inputs / "GSE297984_family.soft.gz")
    title_by_column = {
        re.match(r"(RNA_LS_\d+)", record["title"]).group(1): record["title"]
        for record in chemical_records
    }
    chemical_meta = []
    for column in chemical_scores.index:
        match = re.match(r"RNA_LS_\d+_D(6|14)_(\d+)y_(DMSO|2c|7c)", title_by_column[column])
        if not match:
            raise ValueError(f"Unparsed GSE297984 title: {title_by_column[column]}")
        day, age, treatment = match.groups()
        chemical_meta.append({"sample": column, "day": int(day), "age": int(age), "treatment": treatment})
    chemical_meta = pd.DataFrame(chemical_meta).set_index("sample")
    chemical_scores.to_csv(output / "gse297984-state-scores.csv", index_label="sample")
    chemical_group = chemical_scores.join(chemical_meta).groupby(["age", "day", "treatment"]).mean(numeric_only=True)
    chemical_deltas = []
    for age in [56, 83]:
        for day in [6, 14]:
            for treatment in ["2c", "7c"]:
                delta = chemical_group.loc[(age, day, treatment)] - chemical_group.loc[(age, day, "DMSO")]
                chemical_deltas.append({"age": age, "day": day, "treatment": treatment, **delta.to_dict()})
    chemical_deltas = pd.DataFrame(chemical_deltas)
    chemical_deltas.to_csv(output / "gse297984-donor-day-deltas.csv", index=False)
    positive_by_cocktail = {}
    for treatment in ["2c", "7c"]:
        subset = chemical_deltas[chemical_deltas["treatment"] == treatment]
        donor_positive = subset.groupby("age")["Y_A1"].mean().gt(0).all()
        day_positive = subset.groupby("day")["Y_A1"].mean().gt(0).all()
        positive_by_cocktail[treatment] = bool(donor_positive and day_positive)
    positive_challenge_pass = any(positive_by_cocktail.values())

    # Bulk OSK/O4YRSK two-replicate descriptive contrast.
    osk_counts = read_table(inputs / "GSE297233_raw_counts_matrix.csv.gz", index_col="GeneId")
    osk = filter_and_log_counts(map_human_ensembl_counts(osk_counts, human_map))
    osk_scores = score(osk, model, programmes)
    osk_scores.to_csv(output / "gse297233-state-scores.csv", index_label="sample")
    osk_deltas = pd.DataFrame(
        {
            treatment: osk_scores.loc[[f"{treatment}_D4_1", f"{treatment}_D4_2"]].mean()
            - osk_scores.loc[[f"{treatment}_D0_1", f"{treatment}_D0_2"]].mean()
            for treatment in ["OSK", "O4YRSK"]
        }
    ).T
    osk_deltas.to_csv(output / "gse297233-deltas.csv", index_label="treatment")

    # Single-cell trajectory is pseudobulked by the only available donor/time library.
    trajectories: dict[str, pd.Series] = {}
    for path in sorted(inputs.glob("GSM89865*_filtered_feature_bc_matrix.h5")):
        match = re.match(r"GSM\d+_(GM\d+)_D(\d+)_", path.name)
        if not match:
            raise ValueError(f"Unparsed 10x filename: {path.name}")
        line, day = match.groups()
        trajectories[f"{line}|{day}"] = read_10x_pseudobulk(path)
    trajectory_counts = pd.DataFrame(trajectories).fillna(0)
    trajectory = filter_and_log_counts(trajectory_counts)
    trajectory_scores = score(trajectory, model, programmes)
    trajectory_scores.to_csv(output / "gse297234-state-scores.csv", index_label="sample")
    # GEO: GM00731 is the 96-year donor; GM23815 is the 22-year donor.
    old_line = "GM00731"
    old_baseline = trajectory_scores.loc[f"{old_line}|0"]
    old_deltas = pd.DataFrame(
        {int(day): trajectory_scores.loc[f"{old_line}|{day}"] - old_baseline for day in [3, 7, 10]}
    ).T
    old_deltas.index.name = "day"
    old_deltas.to_csv(output / "gse297234-old-trajectory-deltas.csv")
    youth_days = old_deltas.index[old_deltas["Y_A1"] > 0].tolist()
    risk_days = old_deltas.index[(old_deltas["P"] > 0.05) | (old_deltas["I"] < -0.05)].tolist()
    trajectory_order_pass = bool(youth_days and (not risk_days or min(youth_days) < min(risk_days)))

    # In-vivo adverse challenge; tissue identity is unavailable by design.
    adverse_counts_raw = read_table(inputs / "GSE300625_raw_genecounts.csv.gz")
    adverse = filter_and_log_counts(map_mouse_counts(adverse_counts_raw, mouse_map, orthology))
    adverse_scores = score(adverse, model, programmes, include_identity=False)
    adverse_scores.to_csv(output / "gse300625-state-scores.csv", index_label="sample")
    adverse_records = read_soft_records(inputs / "GSE300625_family.soft.gz")
    adverse_meta = {}
    for record in adverse_records:
        library = record["description"].removeprefix("Library name: ")
        title = record["title"].lower()
        adverse_meta[library] = {
            "tissue": "liver" if "liver" in title else "kidney",
            "treatment": "7c" if title.startswith("7c") else "vehicle",
        }
    adverse_meta = pd.DataFrame(adverse_meta).T
    adverse_intervals = []
    for tissue in ["liver", "kidney"]:
        treated = adverse_scores.loc[adverse_meta.index[(adverse_meta.tissue == tissue) & (adverse_meta.treatment == "7c")]]
        control = adverse_scores.loc[adverse_meta.index[(adverse_meta.tissue == tissue) & (adverse_meta.treatment == "vehicle")]]
        interval = independent_bootstrap_delta(treated, control).reset_index(names="axis")
        interval.insert(0, "tissue", tissue)
        adverse_intervals.append(interval)
    adverse_intervals = pd.concat(adverse_intervals, ignore_index=True)
    adverse_intervals.to_csv(output / "gse300625-bootstrap.csv", index=False)
    adverse_guarded_calls = {"liver": False, "kidney": False}  # I is unavailable, hence never pass.

    # Human ARPE OSK sensitivity; RPE identity is unavailable.
    arpe = read_table(inputs / "GSE304042_5_ARPE_single_triple_OSK_30-1011800743.csv.gz", index_col="Feature.ID")
    arpe_scores = score(arpe, model, programmes, include_identity=False)
    arpe_scores.to_csv(output / "gse304042-state-scores.csv", index_label="sample")
    arpe_delta = arpe_scores.loc[["OSK-mC-a", "OSK-mC-b"]].mean() - arpe_scores.loc[["GFP-a", "GFP-b", "GFP-c"]].mean()
    arpe_delta.rename("OSK_minus_GFP").to_csv(output / "gse304042-osk-delta.csv", header=True)

    # Mouse RPE GSTA4 effect; independent biological animals.
    gsta_counts_raw = read_table(inputs / "GSE304043_6_mRPE_Old-GFP_Old-GSTA4_merged_featureCounts.csv.gz", index_col="Gene name")
    gsta = filter_and_log_counts(map_mouse_counts(gsta_counts_raw, mouse_map, orthology))
    gsta_scores = score(gsta, model, programmes, include_identity=False)
    gsta_scores.to_csv(output / "gse304043-state-scores.csv", index_label="sample")
    gsta_interval = independent_bootstrap_delta(
        gsta_scores.loc[[column for column in gsta_scores.index if "GSTA4" in column]],
        gsta_scores.loc[[column for column in gsta_scores.index if "GFP" in column]],
    )
    gsta_interval.to_csv(output / "gse304043-bootstrap.csv", index_label="axis")

    summary = {
        "model": "Y_A1 frozen at adaptive-y-a1",
        "gse297984_positive_by_cocktail": positive_by_cocktail,
        "gse297984_positive_challenge_pass": positive_challenge_pass,
        "gse297234_first_youthward_day": min(youth_days) if youth_days else None,
        "gse297234_first_identity_or_pluripotency_risk_day": min(risk_days) if risk_days else None,
        "gse297234_trajectory_order_pass": trajectory_order_pass,
        "gse297234_young_minus_old_baseline_y": float(
            trajectory_scores.loc["GM23815|0", "Y_A1"] - trajectory_scores.loc["GM00731|0", "Y_A1"]
        ),
        "gse300625_guarded_calls": adverse_guarded_calls,
        "gse300625_adverse_specificity_pass": not any(adverse_guarded_calls.values()),
        "gse300625_endpoint_informative": False,
        "gse300625_endpoint_limit": "Tissue identity axis unavailable; guarded status is non-evaluable rather than an empirical negative call.",
        "validation_refit": False,
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
