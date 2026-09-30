# Analysis code boundary

This directory contains deterministic command-line stages for manifest validation, preprocessing,
discovery scoring, validation, Adaptive Benchmark B1, post-review Extension C1, figures and claim
verification. Historical runners and repaired outputs are retained separately; do not overwrite an
earlier result directory when reproducing a later stage. Run tests with
`PYTHONPATH=code python -m pytest -q code/tests`.
