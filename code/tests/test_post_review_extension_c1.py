import numpy as np
import pandas as pd

from run_post_review_extension_c1 import (
    equal_family_statistic,
    family_pairwise_agreement,
    pairwise_sign_agreement,
    pasta_predict,
    spearman_brown,
    unique_balanced_splits,
)


def test_contrast_and_equal_family_weighting_are_distinct():
    table = pd.DataFrame(
        {
            "family": ["large"] * 4 + ["small"],
            "a": [1, 1, 1, 1, 1],
            "b": [1, 1, 1, 1, -1],
        }
    )
    contrast = pairwise_sign_agreement(table, ["a", "b"])
    assert contrast.loc[0, "agreement"] == 0.8
    family = family_pairwise_agreement(table, ["a", "b"])
    statistic, per_pair = equal_family_statistic(family)
    assert statistic == 0.5
    assert per_pair.loc[0, "families"] == 2


def test_balanced_splits_remove_even_mirrors():
    splits = unique_balanced_splits(4, maximum=100, seed=1)
    assert len(splits) == 3
    assert all(0 in first for first, _ in splits)
    assert all(set(first).isdisjoint(second) for first, second in splits)


def test_spearman_brown():
    assert np.isclose(spearman_brown(0.5), 2 / 3)
    assert np.isnan(spearman_brown(-1))


def test_pasta_predict_reverses_age_shift_for_display():
    matrix = pd.DataFrame({"sample": [1.0, 2.0]}, index=["g1", "g2"])
    model = {"intercept": 0.5, "coefficients": np.array([1.0, -1.0]), "beta": 2.0}
    # Age-shift link is (0.5 + 1 - 2) * 2 = -1; youth direction is +1.
    assert np.isclose(pasta_predict(matrix, model).loc["sample"], 1.0)
