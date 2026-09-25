import numpy as np
import pandas as pd

from ridge_y import fit_ridge_rank_model, score_ridge_rank_model


def test_ridge_score_has_youthward_direction() -> None:
    genes = [f"g{i}" for i in range(20)]
    a = pd.DataFrame({f"a{i}": np.arange(20) + i for i in range(6)}, index=genes)
    b = pd.DataFrame({f"b{i}": np.arange(20) + i for i in range(6)}, index=genes)
    age_a = pd.Series(range(6), index=a.columns, dtype=float)
    age_b = pd.Series(range(10, 16), index=b.columns, dtype=float)
    model, means = fit_ridge_rank_model(a, b, age_a, age_b, genes, np.array([0.1, 1.0]))
    score = score_ridge_rank_model(a, genes, means, pd.Series(model.coef_, index=genes), model.intercept_)
    assert score.index.equals(a.columns)
    assert score.notna().all()


def test_ridge_missingness_gate() -> None:
    features = [f"g{i}" for i in range(20)]
    expression = pd.DataFrame({"s": np.arange(5)}, index=features[:5])
    score = score_ridge_rank_model(
        expression,
        features,
        pd.Series(0.5, index=features),
        pd.Series(0.0, index=features),
        0.0,
        minimum_coverage=0.60,
    )
    assert score.isna().all()

