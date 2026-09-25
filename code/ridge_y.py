"""Adaptive A1 ridge age-direction model."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import RidgeCV


def fit_ridge_rank_model(
    expression_a: pd.DataFrame,
    expression_b: pd.DataFrame,
    age_a: pd.Series,
    age_b: pd.Series,
    features: list[str],
    alphas: np.ndarray,
) -> tuple[RidgeCV, pd.Series]:
    x_a = expression_a.loc[features].rank(axis=0, pct=True).T
    x_b = expression_b.loc[features].rank(axis=0, pct=True).T
    y_a = (age_a - age_a.mean()) / age_a.std(ddof=1)
    y_b = (age_b - age_b.mean()) / age_b.std(ddof=1)
    x = pd.concat([x_a, x_b])
    y = pd.concat([y_a, y_b])
    model = RidgeCV(alphas=alphas, fit_intercept=True).fit(x, y)
    return model, x.mean(axis=0)


def score_ridge_rank_model(
    expression: pd.DataFrame,
    features: list[str],
    feature_means: pd.Series,
    coefficients: pd.Series,
    intercept: float,
    minimum_coverage: float = 0.60,
) -> pd.Series:
    observed = [feature for feature in features if feature in expression.index]
    coverage = len(observed) / len(features)
    if coverage < minimum_coverage:
        return pd.Series(np.nan, index=expression.columns, dtype=float)
    ranks = expression.rank(axis=0, pct=True)
    x = pd.DataFrame(index=features, columns=expression.columns, dtype=float)
    x.loc[observed] = ranks.loc[observed]
    x = x.apply(lambda column: column.fillna(feature_means), axis=0)
    older_prediction = intercept + x.T.dot(coefficients)
    return -older_prediction


def save_model(
    path: Path,
    features: list[str],
    feature_means: pd.Series,
    coefficients: np.ndarray,
    intercept: float,
    alpha: float,
) -> None:
    table = pd.DataFrame(
        {
            "feature": features,
            "training_mean": feature_means.loc[features].to_numpy(),
            "coefficient": coefficients,
        }
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, index=False)
    metadata = {"intercept": float(intercept), "alpha": float(alpha), "n_features": len(features)}
    path.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")

