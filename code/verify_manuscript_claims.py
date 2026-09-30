#!/usr/bin/env python3
"""Verify headline manuscript claims against frozen aggregate result files.

This is deliberately independent of the scientific runners.  It does not refit a
model or infer a preferred interpretation; it fails if the frozen files no longer
support the exact numerical statements used in the manuscript.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
OUT = ROOT / "receipts" / "manuscript-claim-verification.json"


def close(actual: float, expected: float, tolerance: float = 5e-7) -> None:
    if not math.isclose(float(actual), expected, rel_tol=0.0, abs_tol=tolerance):
        raise AssertionError(f"expected {expected}, observed {actual}")


checks: list[str] = []


def checked(label: str) -> None:
    checks.append(label)


with (RESULTS / "adaptive-y-a1" / "summary.json").open() as handle:
    a1 = json.load(handle)
assert a1["n_features"] == 10_532
for cohort, value in {
    "GSE113957": -0.8465936466546055,
    "GSE226189": -0.45769704410949924,
    "pooled": -0.7301624388033175,
}.items():
    close(a1["outer_spearman_y_vs_age"][cohort], value)
checked("A1 feature count and three outer-CV age correlations")

mptr = pd.read_csv(RESULTS / "adaptive-discovery-a1" / "gse165177-donor-bootstrap.csv").set_index("axis")
for axis, expected in {
    "Y_A1": (0.43695921354137407, 0.16284584676077432, 0.6100423457245076),
    "I": (-0.0035575492352255943, -0.008145576875014382, 0.00017404193716202037),
    "P": (0.0025626610126664438, -0.00130492217044723, 0.005073991363721719),
    "D": (-0.0031008266078036972, -0.003923618873401785, -0.0026320821565173347),
}.items():
    for column, value in zip(("estimate", "ci_low", "ci_high"), expected):
        close(mptr.loc[axis, column], value)
checked("MPTR four-axis estimates and intervals")

with (RESULTS / "adaptive-discovery-a1" / "summary.json").open() as handle:
    discovery = json.load(handle)
close(discovery["gse246954_young_minus_old_y_estimate"], 0.24544955206262722)
close(discovery["gse246954_2c_y_estimate"], -0.3958298811427097)
close(discovery["gse246954_7c_y_estimate"], -0.46265714333161073)
checked("chemical-discovery age control and 2c/7c effects")

chemical = pd.read_csv(RESULTS / "validation" / "gse297984-donor-day-deltas.csv")
assert len(chemical) == 8
assert int((chemical["Y_A1"] < 0).sum()) == 6
checked("human chemical validation contains eight effects, six negative")

trajectory = pd.read_csv(RESULTS / "validation" / "gse297234-old-trajectory-deltas.csv").set_index("day")
for day, expected in {
    3: (0.3076324235836989, -0.08229099952065544, 0.09386927746852275),
    7: (0.561118815982927, -0.05796939467274409, 0.09974106338036398),
    10: (-0.9196534658040866, -0.026682370331084893, 0.044875899673946285),
}.items():
    for column, value in zip(("Y_A1", "I", "P"), expected):
        close(trajectory.loc[day, column], value)
with (RESULTS / "validation" / "summary.json").open() as handle:
    validation = json.load(handle)
assert validation["gse297234_first_youthward_day"] == 3
assert validation["gse297234_first_identity_or_pluripotency_risk_day"] == 3
close(validation["gse297234_young_minus_old_baseline_y"], 0.8555904546463162)
checked("human OSKM trajectory and baseline age control")

adverse = pd.read_csv(RESULTS / "validation" / "gse300625-bootstrap.csv").set_index(["tissue", "axis"])
for tissue, expected in {
    "liver": (-0.11741932044454562, -0.19120848234527196, -0.04129964040580258),
    "kidney": (-0.07520447125426533, -0.17312857532977483, 0.024944847400665552),
}.items():
    for column, value in zip(("estimate", "ci_low", "ci_high"), expected):
        close(adverse.loc[(tissue, "Y_A1"), column], value)
rpe = pd.read_csv(RESULTS / "validation" / "gse304043-bootstrap.csv").set_index("axis")
for axis, expected in {
    "Y_A1": (0.03244199743137326, -0.6375884621234577, 0.5566098079685533),
    "D": (0.004187497047330804, -0.0008930224626552755, 0.008592663809320847),
}.items():
    for column, value in zip(("estimate", "ci_low", "ci_high"), expected):
        close(rpe.loc[axis, column], value)
checked("adverse in-vivo and GSTA4 validation estimates")

with (RESULTS / "target-evaluation" / "summary.json").open() as handle:
    targets = json.load(handle)
assert targets["candidate_universe"] == 9_131
assert targets["missing_targets"] == ["GSTA4", "TOMM70A"]
close(targets["observed_mean_percentile_evaluable_only"], 0.5609097214616873)
close(targets["one_sided_empirical_p"], 0.3323666763332367)
target_rows = pd.read_csv(RESULTS / "target-evaluation" / "known-target-ranks.csv").set_index("gene")
for gene, value in {"DDX21": 0.7929032964625999, "SERBP1": 0.6073814478151353, "PHGDH": 0.2824444201073267}.items():
    close(target_rows.loc[gene, "percentile"], value)
checked("known-target universe, missingness, percentiles, and null result")

with (RESULTS / "benchmark-b1" / "summary-r2.json").open() as handle:
    b1 = json.load(handle)
assert b1["contrast_count"] == 184
assert b1["favourable_family_count"] == 8
assert b1["methods_orienting_both_age_controls"] == 5
assert b1["methods_at_or_above_80_percent_families"] == 0
close(b1["median_pairwise_method_sign_agreement"], 0.7445652173913043)
close(b1["median_pairwise_family_gene_effect_spearman"], 0.0538123679329506)
close(b1["median_leave_one_family_out_gene_consensus_spearman"], 0.1969951702374822)
for method, value in {
    "ridge_rank_a1": -0.7301624388033175,
    "pca50_ridge": -0.7130439071027083,
    "meta_effect_projection": -0.3140637681358926,
    "gse113957_only_ridge": -0.6104703897465925,
    "gse226189_only_ridge": -0.36476481338362987,
}.items():
    close(b1["outer_cv_pooled_spearman"][method], value)
expected_fractions = {"ridge_rank_a1": 0.375, "pca50_ridge": 0.375, "meta_effect_projection": 0.625,
                      "gse113957_only_ridge": 0.25, "gse226189_only_ridge": 0.375}
for method, value in expected_fractions.items():
    close(b1["method_favourable_family_fractions"][method], value)
checked("B1 model count, contrasts, controls, gates, and concordance summaries")

hyper = pd.read_csv(RESULTS / "benchmark-b1" / "final-hyperparameters.csv").set_index("method")
close(hyper.loc["elastic_net_rank", "alpha"], 10.0)
close(hyper.loc["elastic_net_rank", "l1_ratio"], 0.1)
models = __import__("joblib").load(RESULTS / "benchmark-b1" / "final-models.joblib")
elastic = models["elastic_net_rank"]
coef = elastic["coefficients"]
assert int((coef != 0).sum()) == 0
checked("elastic-net selected hyperparameters and zero non-zero coefficients")

with (RESULTS / "post-review-c1" / "summary.json").open() as handle:
    c1 = json.load(handle)
assert c1["contrast_count"] == 184
assert c1["benchmark_family_count"] == 14
close(c1["contrast_weighted_median_pairwise_agreement"], 0.7445652173913043)
close(c1["equal_family_median_pairwise_agreement"], 0.6225143903715332)
close(c1["cluster_bootstrap"]["lower_95"], 0.5151998299319728)
close(c1["cluster_bootstrap"]["upper_95"], 0.7377029778257458)
close(c1["leave_one_family_out_minimum"], 0.5934770357847281)
close(c1["leave_one_family_out_maximum"], 0.6521819526627219)
assert c1["favourable_positive_counts"] == {
    "ridge_rank_a1": 3, "pca50_ridge": 3, "meta_effect_projection": 5,
    "gse113957_only_ridge": 2, "gse226189_only_ridge": 3,
}
close(c1["feature_coverage_minimum"], 0.8513102924420812)
close(c1["common_feature_sensitivity"]["contrast_weighted_median_agreement"], 0.779891304347826)
close(c1["common_feature_sensitivity"]["equal_family_median_agreement"], 0.6486214678178963)
assert c1["pasta"]["implementation_validation"]["status"] == "PASS"
assert c1["pasta"]["both_age_controls_positive"] is True
assert c1["pasta"]["positive_favourable_families"] == 4
checked("C1 family weighting, uncertainty, coverage, common-feature, and Pasta results")

manuscript = (ROOT / "manuscript" / "manuscript.md").read_text()
for forbidden in (
    "prospectively validated",
    "reversed biological age",
    "therapeutically effective",
    "causes rejuvenation",
):
    assert forbidden.lower() not in manuscript.lower()
checked("high-risk unsupported affirmative phrases absent")

payload = {
    "status": "PASS",
    "checks": checks,
    "check_count": len(checks),
    "sources": [
        "results/adaptive-y-a1/summary.json",
        "results/adaptive-discovery-a1/",
        "results/validation/",
        "results/target-evaluation/",
        "results/benchmark-b1/summary-r2.json",
        "results/post-review-c1/summary.json",
    ],
}
OUT.write_text(json.dumps(payload, indent=2) + "\n")
print(json.dumps(payload, indent=2))
