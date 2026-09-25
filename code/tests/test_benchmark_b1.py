import numpy as np
import pandas as pd

from run_benchmark_b1 import effect_weights, fit_method, predict_y, score_expression


def test_effect_projection_zeroes_discordant_features() -> None:
    index = [f"a{i}" for i in range(6)] + [f"b{i}" for i in range(6)]
    y = pd.Series(list(range(6)) + list(range(6)), index=index, dtype=float)
    study = pd.Series(["GSE113957"] * 6 + ["GSE226189"] * 6, index=index)
    x = pd.DataFrame(
        {
            "same": list(range(6)) + list(range(6)),
            "opposite": list(range(6)) + list(range(5, -1, -1)),
        },
        index=index,
        dtype=float,
    )
    weights = effect_weights(x, y, study)
    assert weights["same"] > 0
    assert weights["opposite"] == 0


def test_linear_model_is_oriented_youthward() -> None:
    index = [f"a{i}" for i in range(10)] + [f"b{i}" for i in range(10)]
    study = pd.Series(["GSE113957"] * 10 + ["GSE226189"] * 10, index=index)
    y = pd.Series(list(range(10)) + list(range(10)), index=index, dtype=float)
    x = pd.DataFrame({f"g{j}": y + j * 0.01 for j in range(20)}, index=index)
    model = fit_method("ridge_rank_a1", x, y, study)
    prediction = predict_y(model, x)
    assert prediction.corr(y, method="spearman") < 0


def test_expression_coverage_gate() -> None:
    features = [f"g{i}" for i in range(20)]
    expression = pd.DataFrame({"s": np.arange(5)}, index=features[:5])
    model = {
        "kind": "linear",
        "means": pd.Series(0.5, index=features),
        "coefficients": pd.Series(0.0, index=features),
        "intercept": 0.0,
    }
    assert score_expression(expression, features, model).isna().all()
