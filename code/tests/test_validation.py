import numpy as np
import pandas as pd
import pytest

from score_validation import independent_bootstrap_delta


def test_independent_bootstrap_delta_direction_and_columns() -> None:
    treated = pd.DataFrame({"Y": [2.0, 3.0, 4.0], "D": [0.0, 0.5, 1.0]})
    control = pd.DataFrame({"Y": [0.0, 1.0, 2.0], "D": [1.0, 1.5, 2.0]})
    result = independent_bootstrap_delta(treated, control, n_resamples=100, seed=4)
    assert list(result.columns) == ["estimate", "ci_low", "ci_high"]
    assert result.loc["Y", "estimate"] == pytest.approx(2.0)
    assert result.loc["D", "estimate"] == pytest.approx(-1.0)
    assert np.isfinite(result.to_numpy()).all()
