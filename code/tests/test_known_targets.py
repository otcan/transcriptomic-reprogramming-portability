import numpy as np
import pandas as pd

from evaluate_known_targets import assign_expression_deciles, matched_null


def test_expression_deciles_cover_ten_bins() -> None:
    abundance = pd.Series(np.arange(100), index=[f"g{i:03}" for i in range(100)])
    deciles = assign_expression_deciles(abundance)
    assert sorted(deciles.unique()) == list(range(10))
    assert (deciles.value_counts() == 10).all()


def test_matched_null_is_reproducible_and_excludes_targets() -> None:
    genes = [f"g{i:03}" for i in range(100)]
    percentiles = pd.Series(np.linspace(0, 1, 100), index=genes)
    deciles = assign_expression_deciles(percentiles)
    targets = ["g005", "g015", "g025", "g035", "g045"]
    observed_a, null_a = matched_null(percentiles, deciles, targets, n_resamples=20, seed=7)
    observed_b, null_b = matched_null(percentiles, deciles, targets, n_resamples=20, seed=7)
    assert observed_a == observed_b
    assert np.array_equal(null_a, null_b)


def test_matched_null_rejects_only_targets_passed_to_it() -> None:
    genes = [f"g{i:03}" for i in range(100)]
    percentiles = pd.Series(np.linspace(0, 1, 100), index=genes)
    deciles = assign_expression_deciles(percentiles)
    observed, null = matched_null(percentiles, deciles, ["g005", "g015", "g025"], n_resamples=5, seed=3)
    assert np.isfinite(observed)
    assert np.isfinite(null).all()
