import pandas as pd
import pytest

from repair_b1_cv import standardise_training_age


def test_training_age_standardised_within_cohort() -> None:
    age = pd.Series([10.0, 20.0, 30.0, 50.0, 60.0, 70.0], index=list("abcdef"))
    study = pd.Series(["a", "a", "a", "b", "b", "b"], index=age.index)
    result = standardise_training_age(age, study)
    for cohort in ["a", "b"]:
        values = result[study == cohort]
        assert values.mean() == pytest.approx(0.0)
        assert values.std(ddof=1) == pytest.approx(1.0)
