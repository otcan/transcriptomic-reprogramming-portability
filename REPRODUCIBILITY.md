# Reproducibility

## Environment

Create a Python 3.12 virtual environment and install
`environment/requirements-frozen.txt`. The release bundle records platform details and hashes.

## Public inputs

Download discovery and validation inputs from `inputs/discovery-manifest.tsv` and
`inputs/validation-manifest.tsv`, preserving the recorded filenames. Verify every SHA-256 before
execution. Annotation files are listed in `inputs/annotation-manifest.tsv`. Raw inputs are not part
of the release archive.

## Execution order

The historical order is important:

1. protocol-v1 discovery construction;
2. A1 fit and outer CV;
3. adaptive discovery scoring;
4. frozen target evaluation;
5. temporal validation;
6. Adaptive Benchmark B1;
7. B1 outer-CV r2 repair;
8. post-review Extension C1 in a separate outcome-aware stage;
9. manuscript figures.

Exact commands and the historical integrity boundaries are recorded in `receipts/` and `protocol/`.
Do not overwrite prior outputs during reproduction. Run new outputs in a separate directory or clean
clone.

Extension C1 requires an independently cloned copy of the official Pasta repository pinned to commit
`58bcc7a69ee97f2dc9e3623ac86538c251ddd498`. Verify the three model hashes in
`readiness/contemporary-clock-audit.md`, then run:

```bash
PYTHONPATH=code python code/run_post_review_extension_c1.py --root . --pasta-dir /path/to/pasta
```

## Tests and figure build

```bash
PYTHONPATH=code python -m pytest -q code/tests
PYTHONPATH=code python code/make_manuscript_figures.py
python code/verify_manuscript_claims.py
```

## Expected headline checks

- A1 pooled outer-CV Spearman: `-0.7301624388033175`.
- Five evaluable methods orient both age controls.
- Median pairwise intervention sign agreement: `0.7445652173913043`.
- Methods at or above 80% favourable families: `0`.
- Median pairwise family gene-effect Spearman: `0.0538123679329506`.
- Median leave-one-family-out consensus Spearman: `0.1969951702374822`.
- Equal-family median pairwise sign agreement: `0.6225143903715332`.
- Family-cluster bootstrap 95% interval: `0.5151998299319728` to `0.7377029778257458`.
- Common-feature refit agreement: `0.779891304347826` contrast-weighted and
  `0.6486214678178963` equal-family.
- Exploratory Pasta: both age controls positive; four of eight favourable-labelled families positive.

These checks establish computational reproduction, not biological validation.
