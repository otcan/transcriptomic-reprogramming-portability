"""Deterministic analysis primitives for Paper 2."""

from __future__ import annotations

import gzip
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def bh_fdr(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.clip(q, 0, 1)
    return out


def parse_characteristics(value: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in str(value).split(" | "):
        key, sep, val = item.partition(": ")
        if sep:
            result[key.strip().lower()] = val.strip()
    return result


def load_metadata(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    parsed = table["characteristics"].map(parse_characteristics)
    for key in sorted({key for record in parsed for key in record}):
        table[key] = parsed.map(lambda record: record.get(key, ""))
    return table


def read_gtf_gene_map(path: Path) -> pd.DataFrame:
    rows: list[tuple[str, str, str]] = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]+)"', fields[8]))
            gene_id = attrs.get("gene_id", "").split(".")[0]
            symbol = attrs.get("gene_name", "")
            biotype = attrs.get("gene_biotype", attrs.get("gene_type", ""))
            if gene_id and symbol:
                rows.append((gene_id, symbol, biotype))
    return pd.DataFrame(rows, columns=["gene_id", "symbol", "biotype"]).drop_duplicates("gene_id")


def load_homologene_one_to_one(path: Path) -> pd.DataFrame:
    columns = ["group", "tax_id", "gene_id", "symbol", "protein_gi", "protein_accession"]
    table = pd.read_csv(path, sep="\t", names=columns, dtype=str)
    subset = table[table["tax_id"].isin(["9606", "10090"])]
    counts = subset.groupby(["group", "tax_id"]).size().unstack(fill_value=0)
    valid = counts[(counts.get("9606", 0) == 1) & (counts.get("10090", 0) == 1)].index
    pivot = subset[subset["group"].isin(valid)].pivot(index="group", columns="tax_id", values="symbol")
    return pivot.rename(columns={"9606": "human_symbol", "10090": "mouse_symbol"}).reset_index()


def collapse_expression(matrix: pd.DataFrame, symbols: pd.Series, method: str = "median") -> pd.DataFrame:
    work = matrix.copy()
    work.index = symbols.astype(str).values
    work = work.loc[(work.index != "") & (work.index != "nan")]
    if method == "sum":
        return work.groupby(level=0, sort=True).sum()
    return work.groupby(level=0, sort=True).median()


def log2_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    libraries = counts.sum(axis=0).replace(0, np.nan)
    return np.log2(counts.divide(libraries, axis=1) * 1_000_000 + 0.5)


def linear_age_effect(expression: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    common = [column for column in expression.columns if column in metadata.index]
    expr = expression[common].T
    meta = metadata.loc[common]
    age = pd.to_numeric(meta["age"], errors="raise").to_numpy(dtype=float)
    age_z = (age - age.mean()) / age.std(ddof=1)
    sex = meta["sex"].str.lower().replace({"m": "male", "f": "female"})
    sex_design = pd.get_dummies(sex, prefix="sex", drop_first=True, dtype=float)
    x = np.column_stack([np.ones(len(meta)), age_z, sex_design.to_numpy()])
    y = expr.to_numpy(dtype=float)
    mean = np.nanmean(y, axis=0)
    sd = np.nanstd(y, axis=0, ddof=1)
    valid = np.isfinite(sd) & (sd > 0) & np.isfinite(y).all(axis=0)
    y = (y[:, valid] - mean[valid]) / sd[valid]
    xtx_inv = np.linalg.pinv(x.T @ x)
    beta = xtx_inv @ x.T @ y
    residual = y - x @ beta
    df = len(meta) - x.shape[1]
    sigma2 = np.sum(residual**2, axis=0) / df
    se = np.sqrt(sigma2 * xtx_inv[1, 1])
    t_value = beta[1] / se
    p_value = 2 * stats.t.sf(np.abs(t_value), df)
    genes = expr.columns[valid]
    result = pd.DataFrame(
        {"beta_age": beta[1], "se_age": se, "z_age": t_value, "p_age": p_value}, index=genes
    )
    result["q_age"] = bh_fdr(result["p_age"].to_numpy())
    return result.sort_index()


def select_y_signature(effect_a: pd.DataFrame, effect_b: pd.DataFrame, q_max: float = 0.10) -> dict[str, list[str]]:
    common = effect_a.index.intersection(effect_b.index)
    a = effect_a.loc[common]
    b = effect_b.loc[common]
    eligible = (a["q_age"] <= q_max) & (b["q_age"] <= q_max) & (
        np.sign(a["beta_age"]) == np.sign(b["beta_age"])
    )
    young = common[eligible & (a["beta_age"] < 0)].tolist()
    old = common[eligible & (a["beta_age"] > 0)].tolist()
    return {"young_up": young, "old_up": old}


def rank_score(expression: pd.DataFrame, positive: list[str], negative: list[str] | None = None) -> pd.Series:
    ranks = expression.rank(axis=0, method="average", pct=True)
    pos = [gene for gene in positive if gene in ranks.index]
    if len(pos) < 15 or len(pos) < 0.5 * len(positive):
        return pd.Series(np.nan, index=expression.columns, dtype=float)
    if negative is None:
        return ranks.loc[pos].mean(axis=0) - 0.5
    neg = [gene for gene in negative if gene in ranks.index]
    if len(neg) < 15 or len(neg) < 0.5 * len(negative):
        return pd.Series(np.nan, index=expression.columns, dtype=float)
    return 0.5 * (ranks.loc[pos].mean(axis=0) - ranks.loc[neg].mean(axis=0))


def assign_stratified_folds(metadata: pd.DataFrame, n_folds: int = 5, seed: int = 1729) -> pd.Series:
    rng = np.random.default_rng(seed)
    folds = pd.Series(index=metadata.index, dtype=int)
    age_bin = pd.qcut(pd.to_numeric(metadata["age"]), q=n_folds, labels=False, duplicates="drop")
    strata = metadata["sex"].str.lower().astype(str) + "|" + age_bin.astype(str)
    for _, members in strata.groupby(strata).groups.items():
        values = np.asarray(list(members), dtype=object)
        rng.shuffle(values)
        for offset, sample in enumerate(values):
            folds.loc[sample] = offset % n_folds
    return folds.astype(int)


def percentile_interval(values: pd.DataFrame, n_resamples: int = 10000, seed: int = 1729) -> pd.DataFrame:
    array = values.to_numpy(dtype=float)
    rng = np.random.default_rng(seed)
    means = np.empty((n_resamples, array.shape[1]), dtype=float)
    for idx in range(n_resamples):
        draw = rng.integers(0, len(array), len(array))
        means[idx] = np.nanmean(array[draw], axis=0)
    return pd.DataFrame(
        {
            "estimate": np.nanmean(array, axis=0),
            "ci_low": np.nanpercentile(means, 2.5, axis=0),
            "ci_high": np.nanpercentile(means, 97.5, axis=0),
        },
        index=values.columns,
    )


def parse_go_descendants(obo_path: Path, roots: list[str]) -> dict[str, set[str]]:
    parents: dict[str, set[str]] = defaultdict(set)
    current = ""
    obsolete = False
    with obo_path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line == "[Term]":
                current = ""
                obsolete = False
            elif line.startswith("id: GO:"):
                current = line.split("id: ", 1)[1]
            elif line == "is_obsolete: true":
                obsolete = True
            elif current and line.startswith("is_a: GO:"):
                parents[current].add(line.split()[1])
            elif current and line.startswith("relationship: part_of GO:"):
                parents[current].add(line.split()[2])
            if obsolete and current:
                parents.pop(current, None)

    memo: dict[tuple[str, str], bool] = {}

    def descends(term: str, root: str, trail: frozenset[str] = frozenset()) -> bool:
        key = (term, root)
        if key in memo:
            return memo[key]
        if term == root:
            memo[key] = True
        elif term in trail:
            memo[key] = False
        else:
            memo[key] = any(descends(parent, root, trail | {term}) for parent in parents.get(term, set()))
        return memo[key]

    terms = set(parents) | {parent for values in parents.values() for parent in values}
    return {root: {term for term in terms if descends(term, root)} | {root} for root in roots}


def genes_for_go_terms(gaf_path: Path, descendants: dict[str, set[str]]) -> dict[str, set[str]]:
    term_to_roots = {
        term: [root for root, terms in descendants.items() if term in terms]
        for term in set().union(*descendants.values())
    }
    result = {root: set() for root in descendants}
    opener = gzip.open if str(gaf_path).endswith(".gz") else open
    with opener(gaf_path, "rt") as handle:
        for line in handle:
            if line.startswith("!"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7 or "NOT" in fields[3].split("|") or fields[6] == "ND":
                continue
            symbol, go_id = fields[2], fields[4]
            for root in term_to_roots.get(go_id, []):
                result[root].add(symbol)
    return result

