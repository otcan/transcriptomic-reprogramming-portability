# M3 adaptive benchmark B1 — corrected r2

Date: 2026-09-25  
Inferential status: adaptive, defined after the frozen primary validation failed

## Scope

Benchmark B1 tested whether the failure of the frozen ridge-rank youth direction was specific to one
model or persisted across a locked panel of transcriptomic chronological-age directions. The panel,
contrast set, family units, and two descriptive gates were hashed before comparator outputs were
computed. B1 cannot rescue the failed primary claim and is not blinded, prospective, or confirmatory.

## Execution

- Six methods were attempted on two human skin/fibroblast age references and 184 intervention
  contrasts.
- The initial full run completed locally using single-thread numerical libraries. No paid compute
  was used.
- The selected elastic-net model (`alpha=10`, `l1_ratio=0.1`) had exactly zero non-zero
  coefficients and constant intervention outputs. The implementation failure is retained and the
  method is not replaced.
- Audit found that r1 outer CV standardized age using full-cohort rather than training-fold moments.
  This affected developmental CV target scaling only. It did not affect final models, sample scores,
  contrasts, family summaries, gene-effect summaries, or either gate.
- r1 is preserved. r2 recomputed outer CV for the five evaluable methods using training-fold-only
  within-cohort standardization. Fifteen code tests pass. The corrected ridge-rank predictions match
  the independently developed A1 predictions to maximum absolute difference
  `1.56e-15`.

## Corrected age-reference results

| Method | GSE113957 Spearman | GSE226189 Spearman | Pooled Spearman |
|---|---:|---:|---:|
| Ridge-rank A1 | -0.84659 | -0.45770 | -0.73016 |
| PCA50-ridge | -0.83610 | -0.36583 | -0.71304 |
| Meta-effect projection | -0.59814 | -0.02790 | -0.31406 |
| GSE113957-only ridge | -0.83474 | -0.27645 | -0.61047 |
| GSE226189-only ridge | -0.38042 | -0.40046 | -0.36476 |

Higher model scores were oriented as more youth-associated, so negative score-age correlations are
the expected direction. These are chronological-age associations, not validation of biological age.

## Locked gate outcomes

- All five evaluable methods correctly oriented both independent young-minus-old controls.
- Median pairwise sign agreement across the 184 intervention contrasts was `0.744565`, below the
  locked `0.80` transportability threshold. Broad transportability failed.
- Across eight source-labelled favourable intervention families, method-level positive fractions
  were `0.375`, `0.375`, `0.625`, `0.250`, and `0.375`. No method reached `0.80`; the universal
  direction gate failed.
- Median pairwise Spearman correlation between intervention-family gene-effect vectors was
  `0.053812`.
- Median leave-one-family-out consensus-to-held-family Spearman correlation was `0.196995`.

Source-paper labels summarize the source authors' favourable interpretations; they are not treated
as biological ground truth. Failure of a transcriptomic scalar does not negate functional,
epigenetic, metabolic, or organismal findings in those source studies.

## Evidence

- Lock: `protocol/adaptive-benchmark-b1.md`, `protocol/adaptive-benchmark-b1.json`
- Implementation: `code/run_benchmark_b1.py`, `code/repair_b1_cv.py`
- Tests: `code/tests/test_benchmark_b1.py`, `code/tests/test_b1_cv_repair.py`
- Corrected CV: `results/benchmark-b1/outer-cv-predictions-r2.csv`,
  `results/benchmark-b1/outer-cv-correlations-r2.csv`
- Intervention results: `results/benchmark-b1/contrast-matrix.csv`,
  `results/benchmark-b1/method-sign-agreement.csv`,
  `results/benchmark-b1/method-favourable-fractions.csv`
- Gene effects: `results/benchmark-b1/family-gene-effect-spearman.csv`,
  `results/benchmark-b1/leave-one-family-out-gene-consensus.csv`
- Machine-readable summary: `results/benchmark-b1/summary-r2.json`
- Deviations: `results-and-deviations.md`

