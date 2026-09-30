#!/usr/bin/env python3
"""Create deterministic manuscript figures from committed aggregate results."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "manuscript" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

METHODS = [
    "ridge_rank_a1",
    "pca50_ridge",
    "meta_effect_projection",
    "gse113957_only_ridge",
    "gse226189_only_ridge",
]
METHOD_LABELS = {
    "ridge_rank_a1": "Ridge-rank",
    "pca50_ridge": "PCA50-ridge",
    "meta_effect_projection": "Meta-effect",
    "gse113957_only_ridge": "Cohort A ridge",
    "gse226189_only_ridge": "Cohort B ridge",
}
FAMILY_LABELS = {
    "GSE165177_MPTR": "MPTR",
    "GSE176206_SOKM": "SOKM screen",
    "GSE246954_chemical": "Mouse 2c/7c",
    "GSE297233_OSK": "Human OSK/O4YRSK",
    "GSE297234_early_OSKM": "Human OSKM trajectory",
    "GSE297984_chemical": "Human 2c/7c",
    "GSE304042_ARPE_OSK": "RPE OSK",
    "GSE304043_GSTA4": "RPE GSTA4",
    "GSE165177_failed": "Failed MPTR",
    "GSE176206_other_factors": "Other factor screens",
    "GSE246954_age": "Mouse age control",
    "GSE297234_age": "Human age control",
    "GSE297234_day10": "OSKM day 10",
    "GSE300625_adverse": "Adverse 7c",
}

mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "axes.titlesize": 9.5,
        "axes.labelsize": 8.5,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "figure.dpi": 150,
        "savefig.dpi": 600,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)

# Okabe-Ito-derived palette. Direction is also encoded by signs, labels or position.
BLUE = "#0072B2"
ORANGE = "#D55E00"
SKY = "#56B4E9"
PURPLE = "#CC79A7"
RED = ORANGE
GREEN = BLUE
GREY = "#667085"
PALE = "#F2F4F7"


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.10, 1.10, label, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", bbox_inches="tight")
    plt.close(fig)


def box(ax: plt.Axes, xy: tuple[float, float], width: float, height: float, text: str, color: str) -> None:
    patch = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.012,rounding_size=0.015",
        linewidth=1,
        edgecolor=color,
        facecolor=mpl.colors.to_rgba(color, 0.10),
    )
    ax.add_patch(patch)
    ax.text(xy[0] + width / 2, xy[1] + height / 2, text, ha="center", va="center", fontsize=7)


def figure1_design() -> None:
    fig = plt.figure(figsize=(7.2, 5.6))
    grid = fig.add_gridspec(3, 1, height_ratios=[1.15, 0.85, 1.0], hspace=0.55)
    ax = fig.add_subplot(grid[0])
    panel_label(ax, "a")
    ax.set_title("Temporally separated construction, challenge and robustness analysis", loc="left", pad=8)
    x = np.arange(5)
    ax.plot(x, np.zeros(5), color="#98A2B3", lw=1.4, zorder=1)
    colors = [BLUE, GREY, ORANGE, BLUE, PURPLE]
    ax.scatter(x, np.zeros(5), s=95, color=colors, edgecolor="white", linewidth=1.2, zorder=2)
    top = ["Reference data\nthrough 2024", "Protocol and\ndiscovery freeze", "2025–2026\ntemporal challenge",
           "Adaptive B1\nbenchmark", "Post-review C1\nrobustness"]
    bottom = ["2 age cohorts", "No outcome refitting", "Primary claim failed",
              "5 evaluable methods", "Family weighting + Pasta"]
    for i, (heading, detail) in enumerate(zip(top, bottom)):
        ax.text(i, 0.17, heading, ha="center", va="bottom", fontsize=7.5, fontweight="semibold")
        ax.text(i, -0.17, detail, ha="center", va="top", fontsize=6.8, color=GREY)
    ax.set_xlim(-0.35, 4.35)
    ax.set_ylim(-0.42, 0.48)
    ax.set_axis_off()
    ax.text(0, -0.38, "Publication conclusions were known; temporal separation was neither blinded nor prospective.",
            fontsize=6.8, color=GREY)

    ax = fig.add_subplot(grid[1])
    panel_label(ax, "b")
    ax.set_title("Four axes kept distinct biological questions separate", loc="left", pad=8)
    labels = [
        ("Y", "Chronological-age\ndirection", BLUE),
        ("I", "Somatic\nidentity", SKY),
        ("P", "Pluripotency /\ndedifferentiation", PURPLE),
        ("D", "Measured stress /\ndamage programmes", ORANGE),
    ]
    for index, (symbol, desc, color) in enumerate(labels):
        xpos = index
        ax.text(xpos, 0.61, symbol, ha="center", va="center", fontsize=19, fontweight="bold", color=color)
        ax.text(xpos, 0.20, desc, ha="center", va="center", fontsize=7.4)
        if index < 3:
            ax.axvline(xpos + 0.5, color="#D0D5DD", lw=0.8, ymin=0.15, ymax=0.85)
    ax.set_xlim(-0.5, 3.5)
    ax.set_ylim(-0.15, 1.0)
    ax.set_axis_off()
    ax.text(-0.45, -0.08, "Y is dimensionless; it is not biological age or years rejuvenated.", fontsize=6.8, color=GREY)

    ax = fig.add_subplot(grid[2])
    panel_label(ax, "c")
    ax.set_title("Qualification proceeds from prediction to outcome relevance", loc="left", pad=8)
    steps = [
        ("1", "Age association", "Held-out donors"),
        ("2", "Control orientation", "Young–old contrasts"),
        ("3", "Intervention portability", "Methods and families"),
        ("4", "Outcome triangulation", "State, function, safety"),
    ]
    ax.plot(np.arange(4), np.repeat(0.56, 4), color="#98A2B3", lw=1.4, zorder=1)
    for i, (num, title, desc) in enumerate(steps):
        color = [BLUE, SKY, PURPLE, ORANGE][i]
        ax.scatter(i, 0.56, s=135, color=color, edgecolor="white", linewidth=1.2, zorder=2)
        ax.text(i, 0.56, num, ha="center", va="center", color="white", fontsize=8, fontweight="bold", zorder=3)
        ax.text(i, 0.88, title, ha="center", fontsize=7.5, fontweight="semibold")
        ax.text(i, 0.18, desc, ha="center", fontsize=6.9, color=GREY)
    ax.set_xlim(-0.35, 3.35)
    ax.set_ylim(-0.05, 1.10)
    ax.set_axis_off()
    ax.text(-0.32, -0.02, "Success at one level does not establish validity at the next.", fontsize=6.8, color=GREY)
    save(fig, "figure-1-study-design")


def figure2_primary() -> None:
    pred = pd.read_csv(ROOT / "results/adaptive-y-a1/outer-predictions.csv")
    chem = pd.read_csv(ROOT / "results/validation/gse297984-donor-day-deltas.csv")
    traj = pd.read_csv(ROOT / "results/validation/gse297234-old-trajectory-deltas.csv")
    fig, axes = plt.subplots(
        2,
        2,
        figsize=(7.2, 7.0),
        gridspec_kw={"height_ratios": [1.0, 1.05]},
        constrained_layout=True,
    )

    ax = axes[0, 0]
    panel_label(ax, "a")
    for study, group in pred.groupby("study"):
        ax.scatter(group["age"], group["Y_A1"], s=11, alpha=0.65, label=study)
    ax.set_xlabel("Chronological age (years)")
    ax.set_ylabel("Held-out Y_A1")
    ax.legend(frameon=False, fontsize=6, loc="upper right")
    ax.set_title("Age-reference outer cross-validation")
    ax.text(0.03, 0.03, "Pooled Spearman = −0.730", transform=ax.transAxes, fontsize=7,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8})

    ax = axes[0, 1]
    panel_label(ax, "b")
    endpoint = pd.DataFrame(
        {
            "Endpoint": [
                "Protocol-v1 Y genes",
                "A1 age association",
                "MPTR discovery separation",
                "Human chemical validation",
                "Youth-before-risk trajectory",
                "Adverse 7c specificity",
                "Five-target recovery",
            ],
            "status": [0, 1, 1, -1, -1, 0, -1],
            "label": ["Failed", "Pass", "Pass", "Failed", "Failed", "Not informative", "Failed / incomplete"],
        }
    )
    colors = endpoint.status.map({1: GREEN, 0: GREY, -1: RED})
    ax.barh(np.arange(len(endpoint)), np.ones(len(endpoint)), color=colors, height=0.65)
    ax.set_yticks(np.arange(len(endpoint)), endpoint.Endpoint)
    ax.set_xlim(0, 1)
    ax.set_xticks([])
    ax.invert_yaxis()
    for i, label in enumerate(endpoint.label):
        ax.text(0.5, i, label, color="white", ha="center", va="center", fontsize=7, fontweight="bold")
    ax.set_title("Frozen primary scorecard")

    ax = axes[1, 0]
    panel_label(ax, "c")
    chem["group"] = chem.apply(lambda r: f"{int(r.age)}y d{int(r.day)} {r.treatment}", axis=1)
    colors = [GREEN if value > 0 else RED for value in chem.Y_A1]
    ax.barh(np.arange(len(chem)), chem.Y_A1, color=colors)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks(np.arange(len(chem)), chem.group)
    ax.invert_yaxis()
    ax.set_xlabel("ΔY_A1, treatment − matched control")
    ax.set_title("Post-2024 human chemical challenge\n(6 of 8 effects negative)")

    ax = axes[1, 1]
    panel_label(ax, "d")
    ax.plot(traj.day, traj.Y_A1, marker="o", color=BLUE, label="Y_A1")
    ax.plot(traj.day, traj.I, marker="s", color=SKY, label="Identity (I)")
    ax.plot(traj.day, traj.P, marker="o", color=ORANGE, label="Pluripotency (P)")
    ax.axhline(0, color="black", lw=0.7)
    ax.axhline(0.05, color=GREY, lw=0.6, ls="--")
    ax.axhline(-0.05, color=GREY, lw=0.6, ls="--")
    ax.set_xticks(traj.day)
    ax.set_xlabel("OSKM day versus day 0 (96-year donor)")
    ax.set_ylabel("Within-donor score change")
    ax.legend(frameon=False, fontsize=6, ncol=3, loc="lower center")
    ax.set_title("Youthward movement did not precede risk axes")
    save(fig, "figure-2-frozen-framework")


def figure3_benchmark() -> None:
    corr = pd.read_csv(ROOT / "results/benchmark-b1/outer-cv-correlations-r2.csv")
    age = pd.read_csv(ROOT / "results/benchmark-b1/age-control-means.csv", index_col=0)[METHODS]
    family = pd.read_csv(ROOT / "results/benchmark-b1/favourable-family-means.csv", index_col=0)[METHODS]
    fractions = pd.read_csv(ROOT / "results/benchmark-b1/method-favourable-fractions.csv", index_col=0).loc[METHODS]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 7.0), constrained_layout=True)

    ax = axes[0, 0]
    panel_label(ax, "a")
    pooled = corr[corr.scope == "pooled"].set_index("method").loc[METHODS]
    ax.barh(range(len(METHODS)), pooled.spearman_y_age, color=BLUE)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks(range(len(METHODS)), [METHOD_LABELS[m] for m in METHODS])
    ax.invert_yaxis()
    ax.set_xlabel("Outer-CV Spearman(score, age)")
    ax.set_title("All evaluable models encode chronological age")

    ax = axes[0, 1]
    panel_label(ax, "b")
    sign = np.sign(age.T)
    sns.heatmap(sign, cmap=LinearSegmentedColormap.from_list("sign", [ORANGE, "white", BLUE]), vmin=-1, vmax=1,
                annot=age.T, fmt=".3f", cbar=False, linewidths=0.5, ax=ax)
    ax.set_yticklabels([METHOD_LABELS[m] for m in METHODS], rotation=0)
    ax.set_xticklabels(["Mouse fibroblasts", "Human fibroblasts"], rotation=20, ha="right")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Both young-minus-old controls orient correctly")

    ax = axes[1, 0]
    panel_label(ax, "c")
    ordered = [name for name in FAMILY_LABELS if name in family.index]
    raw = family.loc[ordered, METHODS]
    sign = np.sign(raw)
    annotation = raw.map(lambda x: "+" if x > 0 else ("−" if x < 0 else "0"))
    sns.heatmap(sign, cmap=LinearSegmentedColormap.from_list("sign2", [ORANGE, "white", BLUE]), vmin=-1, vmax=1,
                annot=annotation, fmt="", cbar=False, linewidths=0.5, ax=ax)
    ax.set_yticklabels([FAMILY_LABELS[f] for f in ordered], rotation=0)
    ax.set_xticklabels([METHOD_LABELS[m] for m in METHODS], rotation=35, ha="right")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Source-labelled favourable families disagree")

    ax = axes[1, 1]
    panel_label(ax, "d")
    values = fractions.positive_family_fraction
    ax.barh(range(len(METHODS)), values, color=[GREEN if value >= 0.8 else ORANGE for value in values])
    ax.axvline(0.8, color=RED, ls="--", lw=1, label="Locked gate")
    ax.set_yticks(range(len(METHODS)), [METHOD_LABELS[m] for m in METHODS])
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction of 8 families with positive mean")
    ax.legend(frameon=False, fontsize=6)
    ax.set_title("No model reaches the universal-direction gate")
    save(fig, "figure-3-model-benchmark")


def figure4_concordance() -> None:
    sign = pd.read_csv(ROOT / "results/benchmark-b1/method-sign-agreement.csv", index_col=0).loc[METHODS, METHODS]
    gene = pd.read_csv(ROOT / "results/benchmark-b1/family-gene-effect-spearman.csv", index_col=0)
    loo = pd.read_csv(ROOT / "results/benchmark-b1/leave-one-family-out-gene-consensus.csv")
    family_rel = pd.read_csv(ROOT / "results/post-review-c1/family-split-reliability.csv")
    arm_rel = pd.read_csv(ROOT / "results/post-review-c1/arm-split-reliability.csv")
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 7.2), constrained_layout=True)

    ax = axes[0, 0]
    panel_label(ax, "a")
    sns.heatmap(sign, cmap="Blues", vmin=0.5, vmax=1, annot=True, fmt=".2f", square=True, cbar_kws={"label": "Sign agreement"}, ax=ax)
    ax.set_xticklabels([METHOD_LABELS[m] for m in METHODS], rotation=35, ha="right")
    ax.set_yticklabels([METHOD_LABELS[m] for m in METHODS], rotation=0)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Method sign agreement (184 contrasts)")

    ax = axes[0, 1]
    panel_label(ax, "b")
    sns.heatmap(gene, cmap="vlag", center=0, vmin=-0.25, vmax=0.55, annot=False, square=True,
                cbar_kws={"label": "Spearman"}, ax=ax)
    ax.set_xticklabels([FAMILY_LABELS.get(x, x) for x in gene.columns], rotation=45, ha="right")
    ax.set_yticklabels([FAMILY_LABELS.get(x, x) for x in gene.index], rotation=0)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Gene-effect similarity across families")

    ax = axes[1, 0]
    panel_label(ax, "c")
    labels = [FAMILY_LABELS.get(x, x) for x in loo.held_family]
    ax.barh(range(len(loo)), loo.spearman, color=BLUE)
    ax.axvline(0, color="black", lw=0.7)
    ax.set_yticks(range(len(loo)), labels)
    ax.invert_yaxis()
    ax.set_xlabel("LOFO consensus-to-held-family Spearman")
    ax.set_title("Cross-family consensus transfer")
    ax.text(0.98, 0.04, "Median = 0.197", transform=ax.transAxes, ha="right", fontsize=7)

    ax = axes[1, 1]
    panel_label(ax, "d")
    family_values = family_rel.loc[family_rel.status == "ESTIMATED", "median_spearman"].to_numpy()
    arm_values = arm_rel["median_spearman"].to_numpy()
    rng = np.random.default_rng(44)
    ax.scatter(np.repeat(0, len(family_values)) + rng.uniform(-0.08, 0.08, len(family_values)),
               family_values, s=42, color=BLUE, edgecolor="white", linewidth=0.6, zorder=3)
    ax.scatter(np.repeat(1, len(arm_values)) + rng.uniform(-0.08, 0.08, len(arm_values)),
               arm_values, s=42, marker="s", color=PURPLE, edgecolor="white", linewidth=0.6, zorder=3)
    ax.scatter(2, 0.0538123679329506, s=58, marker="D", color=ORANGE, edgecolor="white", linewidth=0.6, zorder=3)
    for xpos, values in [(0, family_values), (1, arm_values)]:
        ax.hlines(np.median(values), xpos - 0.18, xpos + 0.18, color="black", lw=1.4)
    ax.axhline(0, color="black", lw=0.7)
    ax.set_xticks([0, 1, 2], ["Family split\n(n=4)", "Arm split\n(n=8)", "Cross-family\nmedian"])
    ax.set_ylim(-0.12, 1.0)
    ax.set_ylabel("Raw Spearman correlation")
    ax.set_title("Within-context reliability versus cross-family similarity")
    save(fig, "figure-4-concordance-and-interpretation")


def figure5_post_review() -> None:
    c1 = ROOT / "results" / "post-review-c1"
    summary = json.loads((c1 / "summary.json").read_text())
    loo = pd.read_csv(c1 / "leave-one-family-out-agreement.csv")
    coverage = pd.read_csv(c1 / "coverage-threshold-sensitivity.csv")
    pasta = pd.read_csv(c1 / "pasta-favourable-family-means.csv").set_index("family")
    historical = pd.read_csv(ROOT / "results/benchmark-b1/favourable-family-means.csv", index_col=0)[METHODS]

    fig, axes = plt.subplots(2, 2, figsize=(8.4, 7.2), constrained_layout=True)

    ax = axes[0, 0]
    panel_label(ax, "a")
    bootstrap = summary["cluster_bootstrap"]
    values = [summary["contrast_weighted_median_pairwise_agreement"], summary["equal_family_median_pairwise_agreement"]]
    ax.scatter([0, 1], values, s=70, color=[SKY, BLUE], edgecolor="white", linewidth=0.8, zorder=3,
               label="Observed")
    ax.errorbar(1, bootstrap["median"],
                yerr=[[bootstrap["median"] - bootstrap["lower_95"]],
                      [bootstrap["upper_95"] - bootstrap["median"]]],
                color="black", capsize=5, lw=1.2, marker="_", ms=10, label="Family bootstrap 95% interval")
    ax.axhline(0.8, color=ORANGE, ls="--", lw=1, label="Locked 0.80 gate")
    ax.set_xticks([0, 1], ["Contrast-\nweighted", "Equal-family"])
    ax.set_ylim(0.45, 0.86)
    ax.set_ylabel("Median pairwise sign agreement")
    ax.legend(frameon=False, fontsize=6.2, loc="lower left")
    ax.set_title("Substantial agreement remains below the locked gate")

    ax = axes[0, 1]
    panel_label(ax, "b")
    loo_labels = [FAMILY_LABELS.get(x, x.replace("GSE", "")) for x in loo.held_family]
    y = np.arange(len(loo))
    ax.hlines(y, 0.45, loo.equal_family_median_agreement, color="#D0D5DD", lw=1)
    ax.scatter(loo.equal_family_median_agreement, y, color=BLUE, s=34, zorder=3)
    ax.axvline(0.8, color=ORANGE, ls="--", lw=1)
    ax.set_xlim(0.45, 0.86)
    ax.set_yticks(y, loo_labels)
    ax.invert_yaxis()
    ax.set_xlabel("Equal-family median agreement")
    ax.set_title("No single family explains the result")

    ax = axes[1, 0]
    panel_label(ax, "c")
    categories = ["60%", "70%", "80%", "90%", "Common\n8,427 genes"]
    x = np.arange(len(categories))
    ax.plot(x[:4], coverage.contrast_weighted_median_agreement,
            marker="o", color=SKY, label="Contrast-weighted")
    ax.plot(x[:4], coverage.equal_family_median_agreement,
            marker="s", color=BLUE, label="Equal-family")
    common = summary["common_feature_sensitivity"]
    ax.scatter([x[4]], [common["contrast_weighted_median_agreement"]], marker="o", color=SKY, s=44)
    ax.scatter([x[4]], [common["equal_family_median_agreement"]], marker="s", color=BLUE, s=44)
    ax.axvline(3.5, color="#D0D5DD", lw=0.8)
    ax.axhline(0.8, color=ORANGE, ls="--", lw=1)
    ax.set_xlim(-0.35, 4.35)
    ax.set_ylim(0.45, 0.86)
    ax.set_xticks(x, categories)
    ax.set_xlabel("Eligibility threshold or separate refit")
    ax.set_ylabel("Median agreement")
    ax.legend(frameon=False, fontsize=6.2, loc="lower left")
    ax.set_title("Coverage choices do not restore portability")

    ax = axes[1, 1]
    panel_label(ax, "d")
    combined = historical.copy()
    combined["pasta"] = pasta.loc[combined.index, "pasta_youth_direction"]
    annotation = combined.map(lambda x: "+" if x > 0 else ("−" if x < 0 else "0"))
    sns.heatmap(np.sign(combined), cmap=LinearSegmentedColormap.from_list("sign3", [ORANGE, "white", BLUE]),
                vmin=-1, vmax=1, annot=annotation, fmt="", cbar=False, linewidths=0.5, ax=ax)
    ax.set_yticklabels([FAMILY_LABELS.get(x, x) for x in combined.index], rotation=0)
    ax.set_xticklabels([*[METHOD_LABELS[m] for m in METHODS], "Pasta\n(exploratory)"], rotation=40, ha="right")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Contemporary Pasta remains split (4/8)")
    save(fig, "figure-5-post-review-robustness")


def supplementary_figures() -> None:
    matrix = pd.read_csv(ROOT / "results/benchmark-b1/contrast-matrix.csv")
    values = matrix[METHODS]
    scale = values.abs().median().replace(0, 1)
    plot = values.div(scale).clip(-4, 4)
    fig, ax = plt.subplots(figsize=(7.2, 9.0))
    sns.heatmap(plot, cmap="vlag", center=0, vmin=-4, vmax=4, yticklabels=False,
                cbar_kws={"label": "Effect / method median |effect| (clipped)"}, ax=ax)
    ax.set_xticklabels([METHOD_LABELS[m] for m in METHODS], rotation=35, ha="right")
    ax.set_ylabel("184 prespecified intervention contrasts")
    ax.set_title("Supplementary Fig. 1 | Complete intervention matrix")
    save(fig, "supplementary-figure-1-full-contrast-matrix")

    targets = pd.read_csv(ROOT / "results/target-evaluation/known-target-ranks.csv")
    fig, ax = plt.subplots(figsize=(5.2, 3.2))
    present = targets[targets["rank"].notna()].copy()
    missing = targets[targets["rank"].isna()].copy()
    ax.bar(present.gene, present.percentile, color=BLUE)
    for gene in missing.gene:
        ax.text(gene, 0.04, "not in\nuniverse", ha="center", va="bottom", color=RED, fontsize=7)
    ax.axhline(0.5, color=GREY, ls="--", lw=0.8)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Frozen candidate-rank percentile")
    ax.set_title("Supplementary Fig. 2 | Known later target benchmark")
    save(fig, "supplementary-figure-2-target-recovery")

    g = pd.read_csv(ROOT / "results/adaptive-discovery-a1/gse165177-donor-deltas.csv")
    cols = ["Y_A1", "I", "P", "D"]
    long = g.melt(id_vars=[column for column in g.columns if column not in cols], value_vars=cols,
                  var_name="axis", value_name="delta")
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    sns.stripplot(data=long, x="axis", y="delta", hue="axis", palette=[BLUE, SKY, PURPLE, ORANGE],
                  jitter=0.08, size=7, legend=False, ax=ax)
    ax.axhline(0, color="black", lw=0.8)
    ax.axhline(0.05, color=GREY, lw=0.6, ls="--")
    ax.axhline(-0.05, color=GREY, lw=0.6, ls="--")
    ax.set_xlabel("")
    ax.set_ylabel("Donor-mean MPTR change")
    ax.set_title("Supplementary Fig. 3 | Discovery-stage state separation")
    save(fig, "supplementary-figure-3-discovery-state-axes")


def main() -> None:
    sns.set_theme(
        style="whitegrid",
        context="paper",
        font_scale=0.78,
        rc={"grid.linewidth": 0.35, "grid.color": "#D0D5DD"},
    )
    figure1_design()
    figure2_primary()
    figure3_benchmark()
    figure4_concordance()
    figure5_post_review()
    supplementary_figures()
    manifest = {}
    for path in sorted(OUT.iterdir()):
        if path.is_file():
            import hashlib

            manifest[path.name] = {"bytes": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
