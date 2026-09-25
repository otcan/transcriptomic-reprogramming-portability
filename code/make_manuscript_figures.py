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
}

mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "figure.dpi": 150,
        "savefig.dpi": 600,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)

BLUE = "#2F6B8A"
ORANGE = "#D17A22"
RED = "#B54749"
GREEN = "#3D8061"
GREY = "#667085"
PALE = "#F2F4F7"


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.08, 1.05, label, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")


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
    fig = plt.figure(figsize=(7.2, 6.2))
    grid = fig.add_gridspec(3, 1, height_ratios=[1.2, 0.9, 1.0], hspace=0.44)
    ax = fig.add_subplot(grid[0])
    ax.set_axis_off()
    panel_label(ax, "a")
    ax.set_title("Temporally locked construction, external challenge and adaptive benchmark", loc="left", pad=8)
    box(ax, (0.01, 0.52), 0.27, 0.30, "Discovery through 2024\n2 human age references\n3 reprogramming families", BLUE)
    box(ax, (0.365, 0.52), 0.22, 0.30, "Freeze\nprotocol-v1.0\nthen discovery-a1-v1", GREY)
    box(ax, (0.68, 0.52), 0.30, 0.30, "2025–2026 validation\nchemical, OSK/OSKM, adverse 7c,\nGSTA4", ORANGE)
    ax.annotate("", xy=(0.36, 0.67), xytext=(0.285, 0.67), arrowprops=dict(arrowstyle="->", lw=1.2))
    ax.annotate("", xy=(0.675, 0.67), xytext=(0.59, 0.67), arrowprops=dict(arrowstyle="->", lw=1.2))
    box(ax, (0.365, 0.08), 0.30, 0.24, "Primary claim failed\nNo refitting or endpoint substitution", RED)
    box(ax, (0.70, 0.08), 0.28, 0.24, "Post-primary B1\n6 attempted methods; 5 evaluable\n184 contrasts; 8 families", GREEN)
    ax.annotate("", xy=(0.51, 0.34), xytext=(0.79, 0.51), arrowprops=dict(arrowstyle="->", lw=1.0))
    ax.annotate("", xy=(0.695, 0.20), xytext=(0.67, 0.20), arrowprops=dict(arrowstyle="->", lw=1.0))
    ax.text(0.01, 0.02, "Known publication conclusions were not blinded; expression matrices were held until the frozen discovery release.", fontsize=6.8, color=GREY)

    ax = fig.add_subplot(grid[1])
    ax.set_axis_off()
    panel_label(ax, "b")
    ax.set_title("The frozen state vector kept distinct biological questions separate", loc="left", pad=8)
    labels = [
        ("Y", "youth-associated\nchronological-age direction", BLUE),
        ("I", "retained somatic\nidentity", GREEN),
        ("P", "endogenous pluripotency /\ndedifferentiation", ORANGE),
        ("D", "four measured stress /\ndamage programmes", RED),
    ]
    for index, (symbol, desc, color) in enumerate(labels):
        x = 0.01 + index * 0.247
        box(ax, (x, 0.23), 0.22, 0.48, f"{symbol}\n{desc}", color)
    ax.text(0.01, 0.05, "Y is dimensionless and is not an estimate of biological age or years rejuvenated.", fontsize=7, color=GREY)

    ax = fig.add_subplot(grid[2])
    ax.set_axis_off()
    panel_label(ax, "c")
    ax.set_title("Qualification logic used in the final benchmark paper", loc="left", pad=8)
    steps = [
        ("1", "Age association", "Held-out donors\nand both cohorts"),
        ("2", "Control orientation", "Independent young–old\ncontrasts"),
        ("3", "Intervention portability", "Direction across\nmethods and families"),
        ("4", "Biological triangulation", "Identity, pluripotency,\nfunction and safety"),
    ]
    for i, (num, title, desc) in enumerate(steps):
        x = 0.01 + i * 0.247
        box(ax, (x, 0.20), 0.22, 0.55, f"{num}. {title}\n{desc}", [BLUE, GREEN, ORANGE, RED][i])
        if i < 3:
            ax.annotate("", xy=(x + 0.245, 0.475), xytext=(x + 0.225, 0.475), arrowprops=dict(arrowstyle="->", lw=1))
    ax.text(0.01, 0.04, "Success at an earlier level does not establish validity at a later level.", fontsize=7, color=GREY)
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
    ax.plot(traj.day, traj.I, marker="o", color=GREEN, label="Identity (I)")
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
    sns.heatmap(sign, cmap=LinearSegmentedColormap.from_list("sign", [RED, "white", GREEN]), vmin=-1, vmax=1,
                annot=age.T, fmt=".3f", cbar=False, linewidths=0.5, ax=ax)
    ax.set_yticklabels([METHOD_LABELS[m] for m in METHODS], rotation=0)
    ax.set_xticklabels(["Mouse fibroblasts", "Human fibroblasts"], rotation=20, ha="right")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Both young-minus-old controls orient correctly")

    ax = axes[1, 0]
    panel_label(ax, "c")
    ordered = list(FAMILY_LABELS)
    raw = family.loc[ordered, METHODS]
    sign = np.sign(raw)
    annotation = raw.map(lambda x: "+" if x > 0 else ("−" if x < 0 else "0"))
    sns.heatmap(sign, cmap=LinearSegmentedColormap.from_list("sign2", [RED, "white", GREEN]), vmin=-1, vmax=1,
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
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 7.2), constrained_layout=True)

    ax = axes[0, 0]
    panel_label(ax, "a")
    sns.heatmap(sign, cmap="Blues", vmin=0.5, vmax=1, annot=True, fmt=".2f", square=True, cbar_kws={"label": "Sign agreement"}, ax=ax)
    ax.set_xticklabels([METHOD_LABELS[m] for m in METHODS], rotation=35, ha="right")
    ax.set_yticklabels([METHOD_LABELS[m] for m in METHODS], rotation=0)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Method agreement across 184 contrasts")

    ax = axes[0, 1]
    panel_label(ax, "b")
    sns.heatmap(gene, cmap="vlag", center=0, vmin=-0.25, vmax=0.55, annot=False, square=True,
                cbar_kws={"label": "Spearman"}, ax=ax)
    ax.set_xticklabels([FAMILY_LABELS.get(x, x) for x in gene.columns], rotation=45, ha="right")
    ax.set_yticklabels([FAMILY_LABELS.get(x, x) for x in gene.index], rotation=0)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Family gene-effect directions are weakly concordant")

    ax = axes[1, 0]
    panel_label(ax, "c")
    labels = [FAMILY_LABELS.get(x, x) for x in loo.held_family]
    ax.barh(range(len(loo)), loo.spearman, color=BLUE)
    ax.axvline(0, color="black", lw=0.7)
    ax.set_yticks(range(len(loo)), labels)
    ax.invert_yaxis()
    ax.set_xlabel("LOFO consensus-to-held-family Spearman")
    ax.set_title("A cross-family consensus transfers only weakly")
    ax.text(0.98, 0.04, "Median = 0.197", transform=ax.transAxes, ha="right", fontsize=7)

    ax = axes[1, 1]
    panel_label(ax, "d")
    ax.set_axis_off()
    ax.set_title("What the benchmark supports", loc="left")
    statements = [
        (GREEN, "Supported", "Age-associated directions orient\nsimple age controls."),
        (ORANGE, "Not portable", "Signs and family calls diverge\nunder reprogramming."),
        (RED, "Not established", "Biological age reversal, causality,\nefficacy or safety."),
    ]
    for i, (color, title, desc) in enumerate(statements):
        y = 0.72 - i * 0.27
        box(ax, (0.02, y), 0.94, 0.18, f"{title}\n{desc}", color)
    save(fig, "figure-4-concordance-and-interpretation")


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
    sns.stripplot(data=long, x="axis", y="delta", hue="axis", palette=[BLUE, GREEN, ORANGE, RED],
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
    supplementary_figures()
    manifest = {}
    for path in sorted(OUT.iterdir()):
        if path.is_file():
            import hashlib

            manifest[path.name] = {"bytes": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
