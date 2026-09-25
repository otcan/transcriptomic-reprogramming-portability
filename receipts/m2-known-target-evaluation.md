# M2 known-target evaluation — frozen-rank result

Date: 2026-09-25  
Discovery release evaluated: tag `discovery-a1-v1`, commit
`cdece0dec297800e7e3ab1aa1688f7f6d8a9d6b9`

The five target names were first looked up only after the complete 9,131-gene candidate rank was
immutable. The complete benchmark could not be evaluated: **GSTA4 and TOMM70A are absent** from the
frozen intersection universe. They were reported as missing and were not replaced.

Exact evaluable ranks are:

- DDX21: 1,892 (percentile 0.79290)
- SERBP1: 3,586 (percentile 0.60738)
- PHGDH: 6,553 (percentile 0.28244)

The mean percentile of these three evaluable genes is 0.56091. Against 100,000 sets matched to their
discovery-reference expression deciles, the one-sided empirical p-value is 0.33237 (null 95%
interval 0.17238 to 0.80345). This does not support enrichment. Because two prespecified targets are
missing, the original five-target endpoint is additionally marked incomplete.

Evidence: `results/target-evaluation/known-target-ranks.csv`, `matched-null.csv`, and `summary.json`.
The null implementation was committed at `9d14d66da94b0e1572ade6d5828490e8aab24d4d`
before this first lookup.
