#!/usr/bin/env python3
"""Run the locked bulk-discovery analysis and freeze the candidate rank."""

from __future__ import annotations

import argparse
import gzip
import json
import re
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from core import (
    assign_stratified_folds,
    collapse_expression,
    genes_for_go_terms,
    linear_age_effect,
    load_homologene_one_to_one,
    load_metadata,
    log2_cpm,
    parse_go_descendants,
    percentile_interval,
    rank_score,
    read_gtf_gene_map,
    select_y_signature,
)


ROOTS = ["GO:0006979", "GO:0006974", "GO:0034976", "GO:0090398"]


def load_gse113957(path: Path, metadata: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
    symbols = raw["Annotation/Divergence"].astype(str).str.split("|").str[0]
    sample_meta = metadata[(metadata["series"] == "GSE113957") & (metadata["disease"].str.lower() == "normal")].copy()
    sample_meta = sample_meta.set_index("title")
    samples = [name for name in raw.columns if name in sample_meta.index]
    fpkm = collapse_expression(raw[samples].astype(float), symbols, "median")
    keep = (fpkm >= 1).mean(axis=1) >= 0.20
    return np.log2(fpkm.loc[keep] + 0.5), fpkm.loc[keep], sample_meta.loc[samples]


def load_gse226189(path: Path, metadata: pd.DataFrame, human_map: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sample_meta = metadata[metadata["series"] == "GSE226189"].copy().set_index("accession")
    series: list[pd.Series] = []
    with tarfile.open(path) as archive:
        members = sorted(
            (member for member in archive.getmembers() if member.name.endswith("_geneCOUNT.txt.gz")),
            key=lambda member: member.name,
        )
        for member in members:
            accession = member.name.split("_", 1)[0]
            extracted = archive.extractfile(member)
            if extracted is None:
                raise RuntimeError(f"Cannot extract {member.name}")
            with gzip.open(extracted, "rt") as handle:
                table = pd.read_csv(handle, sep="\t")
            values = pd.Series(table.iloc[:, 1].to_numpy(dtype=float), index=table.iloc[:, 0].str.split(".").str[0])
            values.name = accession
            series.append(values)
    counts = pd.concat(series, axis=1).fillna(0)
    symbol_map = human_map.set_index("gene_id")["symbol"]
    symbols = counts.index.to_series().map(symbol_map).fillna("")
    counts = collapse_expression(counts, symbols, "sum")
    cpm = counts.divide(counts.sum(axis=0), axis=1) * 1_000_000
    keep = (cpm >= 1).mean(axis=1) >= 0.20
    counts = counts.loc[keep]
    ordered_meta = sample_meta.loc[counts.columns].copy()
    ordered_meta.index = counts.columns
    return log2_cpm(counts), counts, ordered_meta


def load_geo_expression(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
    annotation_count = 12
    return collapse_expression(raw.iloc[:, annotation_count:].astype(float), raw["Feature"], "median")


def construct_ip(sendai: pd.DataFrame, damage_genes: set[str]) -> tuple[list[str], list[str], pd.DataFrame]:
    lines = ["N2", "N3", "Y1", "Y2", "O1", "O2"]
    changes: dict[str, pd.Series] = {}
    baselines: dict[str, pd.Series] = {}
    late_values: dict[str, pd.Series] = {}
    for line in lines:
        baseline = f"{line}_Fib_Sendai_Exp2"
        late = [f"{line}_d47_SSEA4_Sendai_Exp2", f"{line}_d54_SSEA4_Sendai_Exp2"]
        missing = [sample for sample in [baseline, *late] if sample not in sendai.columns]
        if missing:
            raise ValueError(f"Missing I/P anchors: {missing}")
        late_median = sendai[late].median(axis=1)
        changes[line] = late_median - sendai[baseline]
        baselines[line] = sendai[baseline]
        late_values[line] = late_median
    delta = pd.DataFrame(changes)
    baseline_expr = pd.DataFrame(baselines)
    late_expr = pd.DataFrame(late_values)
    p_mask = (delta.median(axis=1) >= 1) & ((delta >= 1).sum(axis=1) >= 5) & ((late_expr >= 1).sum(axis=1) >= 5)
    pluripotency = sorted(set(delta.index[p_mask]) - {"POU5F1", "SOX2", "KLF4", "MYC"})
    i_mask = (delta.median(axis=1) <= -1) & ((delta <= -1).sum(axis=1) >= 5) & ((baseline_expr >= 1).sum(axis=1) >= 5)
    identity = sorted(set(delta.index[i_mask]) - set(pluripotency) - damage_genes)
    return identity, pluripotency, delta


def score_axes(
    expression: pd.DataFrame,
    y_signature: dict[str, list[str]],
    identity: list[str],
    pluripotency: list[str],
    damage: dict[str, list[str]],
) -> pd.DataFrame:
    scores = pd.DataFrame(index=expression.columns)
    scores["Y"] = rank_score(expression, y_signature["young_up"], y_signature["old_up"])
    scores["I"] = rank_score(expression, identity)
    scores["P"] = rank_score(expression, pluripotency)
    d_columns: list[str] = []
    for go_id, genes in damage.items():
        column = f"D_{go_id.replace(':', '_')}"
        scores[column] = rank_score(expression, genes)
        d_columns.append(column)
    scores["D"] = scores[d_columns].mean(axis=1) if not scores[d_columns].isna().any(axis=None) else scores[d_columns].mean(axis=1, skipna=False)
    return scores


def parse_final_state(name: str, state: str) -> tuple[str, int, str] | None:
    pattern = rf"^(O[123])_{state}_(\d+)days_(exp\d+)$"
    match = re.match(pattern, name)
    if not match:
        return None
    return match.group(1), int(match.group(2)), match.group(3)


def matched_mptr_differences(expression: pd.DataFrame, scores: pd.DataFrame, state: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    score_rows: list[pd.Series] = []
    gene_rows: list[pd.Series] = []
    for sample in expression.columns:
        parsed = parse_final_state(sample, state)
        if parsed is None:
            continue
        donor, day, experiment = parsed
        control = f"{donor}_negative_control_{day}days_{experiment}"
        if control not in expression.columns:
            continue
        score_delta = scores.loc[sample] - scores.loc[control]
        score_delta.name = f"{donor}|{day}|{experiment}"
        score_rows.append(score_delta)
        gene_delta = expression[sample] - expression[control]
        gene_delta.name = score_delta.name
        gene_rows.append(gene_delta)
    return pd.DataFrame(score_rows), pd.concat(gene_rows, axis=1)


def load_gse246954(path: Path, mouse_map: pd.DataFrame, orthology: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_excel(path, sheet_name="as.data.frame(count_matrix)")
    counts = raw.set_index("geneID")
    symbol_map = mouse_map.set_index("gene_id")["symbol"]
    mouse_symbols = counts.index.to_series().str.split(".").str[0].map(symbol_map).fillna("")
    mouse_counts = collapse_expression(counts, mouse_symbols, "sum")
    ortho_map = orthology.set_index("mouse_symbol")["human_symbol"]
    human_symbols = mouse_counts.index.to_series().map(ortho_map).fillna("")
    human_counts = collapse_expression(mouse_counts, human_symbols, "sum")
    cpm = human_counts.divide(human_counts.sum(axis=0), axis=1) * 1_000_000
    keep = (cpm >= 1).mean(axis=1) >= 0.20
    return log2_cpm(human_counts.loc[keep]), human_counts.loc[keep]


def paired_chemical_differences(expression: pd.DataFrame, scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    score_rows: list[pd.Series] = []
    gene_rows: list[pd.Series] = []
    for replicate in range(1, 5):
        for treatment in ["2c", "7c"]:
            treated, control = f"A{replicate}_{treatment}", f"A{replicate}control"
            delta = scores.loc[treated] - scores.loc[control]
            delta.name = f"A{replicate}|{treatment}"
            score_rows.append(delta)
            gene = expression[treated] - expression[control]
            gene.name = delta.name
            gene_rows.append(gene)
    return pd.DataFrame(score_rows), pd.concat(gene_rows, axis=1)


def zscore_columns(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.apply(lambda column: (column - column.mean()) / column.std(ddof=1), axis=0)


def candidate_rank(
    effect_a: pd.DataFrame,
    effect_b: pd.DataFrame,
    mptr_positive: pd.DataFrame,
    chemical_positive: pd.DataFrame,
    mptr_negative: pd.DataFrame,
) -> pd.DataFrame:
    universe = effect_a.index.intersection(effect_b.index).intersection(mptr_positive.index).intersection(chemical_positive.index)
    a, b = effect_a.loc[universe], effect_b.loc[universe]
    concordant = np.sign(a["beta_age"]) == np.sign(b["beta_age"])
    meta_z = (a["z_age"] + b["z_age"]) / np.sqrt(2)
    age_direction = np.sign((a["beta_age"] + b["beta_age"]) / 2)
    age_raw = meta_z.abs().where(concordant, 0)

    mp = zscore_columns(mptr_positive.loc[universe]).median(axis=1)
    chem = zscore_columns(chemical_positive.loc[universe]).median(axis=1)
    mp_signed = -age_direction * mp
    chem_signed = -age_direction * chem
    reversal_raw = pd.concat([mp_signed, chem_signed], axis=1).median(axis=1).clip(lower=0)

    eligible = pd.concat([zscore_columns(mptr_positive.loc[universe]), zscore_columns(chemical_positive.loc[universe])], axis=1)
    eligible_reversal = eligible.mul(-age_direction, axis=0).gt(0).mean(axis=1)
    if mptr_negative.empty:
        negative_reversal = pd.Series(0.0, index=universe)
    else:
        negative_reversal = zscore_columns(mptr_negative.loc[universe]).mul(-age_direction, axis=0).gt(0).mean(axis=1)
    consistency_raw = (eligible_reversal - negative_reversal).clip(0, 1)

    rank = pd.DataFrame(index=universe)
    rank["age_beta_meta"] = (a["beta_age"] + b["beta_age"]) / 2
    rank["A"] = age_raw.rank(method="average", pct=True)
    rank["R"] = reversal_raw.rank(method="average", pct=True)
    rank["C"] = consistency_raw.rank(method="average", pct=True)
    rank["T"] = 0.40 * rank["A"] + 0.40 * rank["R"] + 0.20 * rank["C"]
    rank["rank"] = rank["T"].rank(method="max", ascending=False).astype(int)
    rank["percentile"] = rank["T"].rank(method="average", pct=True)
    return rank.sort_values(["rank", "T"], ascending=[True, False])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    args = parser.parse_args()
    root = args.root.resolve()
    inputs = root / "inputs"
    results = root / "results" / "discovery"
    results.mkdir(parents=True, exist_ok=True)

    metadata = load_metadata(root / "receipts" / "sample-metadata.tsv")
    human_map = read_gtf_gene_map(inputs / "annotations" / "Homo_sapiens.GRCh38.111.gtf.gz")
    mouse_map = read_gtf_gene_map(inputs / "annotations" / "Mus_musculus.GRCm39.111.gtf.gz")
    orthology = load_homologene_one_to_one(inputs / "annotations" / "homologene_build68.data")

    ref1_expr, _, ref1_meta = load_gse113957(inputs / "discovery-data" / "GSE113957_fpkm.txt.gz", metadata)
    ref2_expr, _, ref2_meta = load_gse226189(inputs / "discovery-data" / "GSE226189_RAW.tar", metadata, human_map)
    ref1_meta = ref1_meta.copy()
    ref2_meta = ref2_meta.copy()
    ref2_meta["age"] = ref2_meta["age (years)"]

    fold1 = assign_stratified_folds(ref1_meta, seed=1729)
    fold2 = assign_stratified_folds(ref2_meta, seed=1729)
    predictions: list[pd.DataFrame] = []
    fold_gene_counts: list[dict[str, int]] = []
    for fold in range(5):
        effect1 = linear_age_effect(ref1_expr.loc[:, fold1 != fold], ref1_meta.loc[fold1 != fold])
        effect2 = linear_age_effect(ref2_expr.loc[:, fold2 != fold], ref2_meta.loc[fold2 != fold])
        signature = select_y_signature(effect1, effect2)
        fold_gene_counts.append({"fold": fold, "young_up": len(signature["young_up"]), "old_up": len(signature["old_up"])})
        if min(len(signature["young_up"]), len(signature["old_up"])) < 25:
            continue
        for study, expression, meta, folds in [
            ("GSE113957", ref1_expr, ref1_meta, fold1),
            ("GSE226189", ref2_expr, ref2_meta, fold2),
        ]:
            held = folds[folds == fold].index
            score = rank_score(expression[held], signature["young_up"], signature["old_up"])
            predictions.append(pd.DataFrame({"sample": held, "study": study, "age": pd.to_numeric(meta.loc[held, "age"]).values, "Y": score.loc[held].values, "fold": fold}))
    prediction_table = pd.concat(predictions, ignore_index=True) if predictions else pd.DataFrame(columns=["sample", "study", "age", "Y", "fold"])
    prediction_table.to_csv(results / "y-cross-validation-predictions.csv", index=False)
    pd.DataFrame(fold_gene_counts).to_csv(results / "y-cross-validation-feature-counts.csv", index=False)

    correlations: dict[str, float] = {}
    if not prediction_table.empty:
        correlations["pooled"] = float(stats.spearmanr(prediction_table["Y"], prediction_table["age"]).statistic)
        for study, group in prediction_table.groupby("study"):
            correlations[study] = float(stats.spearmanr(group["Y"], group["age"]).statistic)

    final_effect1 = linear_age_effect(ref1_expr, ref1_meta)
    final_effect2 = linear_age_effect(ref2_expr, ref2_meta)
    final_signature = select_y_signature(final_effect1, final_effect2)
    final_effect1.to_csv(results / "age-effects-gse113957.csv", index_label="gene")
    final_effect2.to_csv(results / "age-effects-gse226189.csv", index_label="gene")
    pd.DataFrame(
        [(gene, "young_up") for gene in final_signature["young_up"]] + [(gene, "old_up") for gene in final_signature["old_up"]],
        columns=["gene", "direction"],
    ).to_csv(results / "y-signature.csv", index=False)

    descendants = parse_go_descendants(inputs / "annotations" / "go-basic-2024-01-17.obo", ROOTS)
    d_human = genes_for_go_terms(inputs / "annotations" / "goa_human-2024-01-17.gaf.gz", descendants)
    damage = {root_id: sorted(genes) for root_id, genes in d_human.items()}
    damage_union = set().union(*d_human.values())
    pd.DataFrame([(root_id, gene) for root_id, genes in damage.items() for gene in genes], columns=["go_id", "gene"]).to_csv(results / "damage-programmes.csv", index=False)

    sendai = load_geo_expression(inputs / "discovery-data" / "GSE165176_Log2_RPM_Sendai_reprogramming.txt.gz")
    identity, pluripotency, ip_delta = construct_ip(sendai, damage_union)
    pd.DataFrame({"gene": identity}).to_csv(results / "identity-programme.csv", index=False)
    pd.DataFrame({"gene": pluripotency}).to_csv(results / "pluripotency-programme.csv", index=False)
    ip_delta.to_csv(results / "ip-anchor-donor-deltas.csv", index_label="gene")

    part1 = load_geo_expression(inputs / "discovery-data" / "GSE165177_Log2_RPM_Transient_reprogramming.txt.gz")
    part2 = load_geo_expression(inputs / "discovery-data" / "GSE165177_Log2_RPM_Transient_reprogramming_part2_170621.txt.gz")
    mptr = pd.concat([part1, part2], axis=1)
    if mptr.columns.duplicated().any():
        raise ValueError("Duplicated GSE165177 sample columns")
    mptr_scores = score_axes(mptr, final_signature, identity, pluripotency, damage)
    mptr_scores.to_csv(results / "gse165177-state-scores.csv", index_label="sample")
    mptr_score_delta, mptr_gene_delta = matched_mptr_differences(mptr, mptr_scores, "transiently_reprogrammed")
    failed_score_delta, failed_gene_delta = matched_mptr_differences(mptr, mptr_scores, "failed_to_transiently_reprogram")
    mptr_donor = mptr_score_delta.assign(donor=mptr_score_delta.index.str.split("|").str[0]).groupby("donor").mean(numeric_only=True)
    mptr_interval = percentile_interval(mptr_donor, seed=1729)
    mptr_score_delta.to_csv(results / "gse165177-successful-paired-deltas.csv", index_label="pair")
    failed_score_delta.to_csv(results / "gse165177-failed-paired-deltas.csv", index_label="pair")
    mptr_interval.to_csv(results / "gse165177-donor-bootstrap.csv", index_label="axis")

    chemical, _ = load_gse246954(root / "data" / "working" / "GSE246954_raw_genecounts.xls", mouse_map, orthology)
    chemical_scores = score_axes(chemical, final_signature, identity, pluripotency, damage)
    chemical_scores.to_csv(results / "gse246954-state-scores.csv", index_label="sample")
    chemical_score_delta, chemical_gene_delta = paired_chemical_differences(chemical, chemical_scores)
    chemical_score_delta.to_csv(results / "gse246954-old-paired-deltas.csv", index_label="pair")
    chemical_intervals: list[pd.DataFrame] = []
    for treatment in ["2c", "7c"]:
        subset = chemical_score_delta[chemical_score_delta.index.str.endswith(treatment)]
        interval = percentile_interval(subset, seed=1729).reset_index(names="axis")
        interval.insert(0, "treatment", treatment)
        chemical_intervals.append(interval)
    pd.concat(chemical_intervals).to_csv(results / "gse246954-paired-bootstrap.csv", index=False)

    ranking = candidate_rank(final_effect1, final_effect2, mptr_gene_delta, chemical_gene_delta, failed_gene_delta)
    ranking.to_csv(results / "candidate-rank.csv", index_label="gene")

    y_pass = (
        correlations.get("pooled", 1) <= -0.35
        and correlations.get("GSE113957", 1) <= -0.15
        and correlations.get("GSE226189", 1) <= -0.15
        and all(min(row["young_up"], row["old_up"]) >= 25 for row in fold_gene_counts)
    )
    separation_pass = bool(
        mptr_interval.loc["Y", "ci_low"] > 0
        and mptr_interval.loc["I", "ci_low"] > -0.05
        and mptr_interval.loc["P", "ci_high"] < 0.05
    )
    summary = {
        "y_cross_validation_correlations": correlations,
        "y_cross_validation_pass": y_pass,
        "y_final_genes": {key: len(value) for key, value in final_signature.items()},
        "identity_genes": len(identity),
        "pluripotency_genes": len(pluripotency),
        "damage_genes": {key: len(value) for key, value in damage.items()},
        "mptr_matched_pairs": len(mptr_score_delta),
        "mptr_donors": len(mptr_donor),
        "discovery_separation_pass": separation_pass,
        "candidate_universe": len(ranking),
        "validation_unsealed": False,
    }
    (results / "discovery-summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
