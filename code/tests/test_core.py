import numpy as np
import pandas as pd
import pytest

from core import bh_fdr, percentile_interval, rank_score, select_y_signature


def test_bh_is_bounded_and_monotonic_by_p() -> None:
    p = np.array([0.04, 0.001, 0.02, 0.8])
    q = bh_fdr(p)
    assert np.all((q >= 0) & (q <= 1))
    assert list(q[np.argsort(p)]) == sorted(q)


def test_signed_rank_score_direction() -> None:
    genes = [f"g{i}" for i in range(40)]
    expression = pd.DataFrame({"young": np.arange(40), "old": np.arange(40)[::-1]}, index=genes)
    score = rank_score(expression, genes[20:], genes[:20])
    assert score["young"] > 0
    assert score["old"] < 0


def test_signature_requires_concordant_sign_and_q() -> None:
    a = pd.DataFrame({"beta_age": [-1, 1, 1], "q_age": [0.01, 0.01, 0.01]}, index=list("abc"))
    b = pd.DataFrame({"beta_age": [-2, 2, -2], "q_age": [0.02, 0.02, 0.02]}, index=list("abc"))
    signature = select_y_signature(a, b)
    assert signature == {"young_up": ["a"], "old_up": ["b"]}


def test_bootstrap_returns_columns() -> None:
    values = pd.DataFrame({"Y": [0.1, 0.2, 0.3], "I": [-0.01, 0.0, 0.01]})
    result = percentile_interval(values, n_resamples=100, seed=1)
    assert list(result.columns) == ["estimate", "ci_low", "ci_high"]
    assert result.loc["Y", "estimate"] == pytest.approx(0.2)
