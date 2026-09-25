# Adaptive benchmark B1 — transportability after primary failure

Date defined: 2026-09-25  
Status: lock before fitting or applying any B1 comparator  
Primary frozen validation outcomes: already visible and failed

## Purpose and inferential status

B1 asks whether the failed temporal transport of Y_A1 is model-specific or shared by several
reasonable transcriptomic age directions. It is explicitly post hoc: all Y_A1 discovery and
validation outcomes were known before B1. B1 cannot rescue or replace the preregistered framework
claim and cannot be described as blinded, prospective, or confirmatory.

## Common inputs and representation

- Training inputs remain GSE113957 and GSE226189 only, with the same donors, common 10,532
  protein-coding genes, within-sample percentile ranks, and within-cohort standardized age used by
  Adaptive Amendment A1.
- The unchanged fixed five-fold donor splits are used for outer development estimates.
- Every final comparator is fitted on both complete references. No intervention, validation,
  target, identity, pluripotency, damage, or publication-outcome variable enters fitting.
- New samples require at least 60% feature coverage. Missing features are filled with frozen
  training means.
- Every model is oriented so a higher score is more youth-associated. Outputs are dimensionless;
  they are not years of biological age.

## Frozen comparator panel

1. **Ridge-rank:** the existing Y_A1 model, unchanged.
2. **Elastic-net-rank:** ElasticNetCV over `l1_ratio = [0.1, 0.5, 0.9]` and 25 alphas from
   `10^-4` to `10^1`; fivefold CV occurs inside each outer training fold and uses mean squared
   error. Maximum iterations 50,000; seed 1729.
3. **PCA50-ridge:** center training rank features; retain `min(50, n_train-1)` PCA components;
   choose ridge alpha by generalized leave-one-out from `10**[-2,-1.5,...,6]`; seed 1729.
4. **Meta-effect projection:** fixed weight for each gene is the equally weighted mean of the two
   cohort standardized age-effect z statistics, set to zero when effect signs disagree. The score
   is the negative weighted mean of centered within-sample ranks, divided by the sum of absolute
   nonzero weights. This has no fitted hyperparameter.
5. **GSE113957-only ridge-rank:** A1 ridge procedure trained on GSE113957 only.
6. **GSE226189-only ridge-rank:** A1 ridge procedure trained on GSE226189 only.

All hyperparameters and final coefficients are exported. An implementation failure is reported; a
method is not replaced after seeing intervention outcomes.

## Contrasts

Use the exact biological units, normalisation, orthology, intervention-control pairings, and nested
replicate collapses already recorded for:

- GSE165177 successful and failed MPTR;
- GSE246954 2c, 7c, and young-minus-old untreated control;
- GSE176206 factor subsets by experimental pool;
- GSE297984 donor/day 2c and 7c;
- GSE297233 OSK and O4YRSK;
- GSE297234 96-year donor days 3, 7, and 10 versus day 0, plus young-minus-old day-0 control;
- GSE300625 liver and kidney 7c;
- GSE304042 ARPE OSK;
- GSE304043 GSTA4.

The source papers' favourable interpretations are labels for concordance auditing, not biological
ground truth. GSE300625 is the adverse family. Age controls are separate from intervention labels.

## Fixed summaries

1. Outer-CV Spearman correlation with age for every comparator, pooled and by reference cohort.
2. Complete method-by-contrast Delta youth-score matrix; no contrast is removed for disagreement.
3. Pairwise method sign agreement and Spearman correlation across contrasts.
4. For each method, the fraction of source-labelled favourable study families whose mean Delta is
   positive, with study family—not replicate, cell, or contrast—as the unit.
5. Correct orientation of the GSE246954 and GSE297234 young-minus-old controls.
6. Gene-level effect-vector Spearman correlations across study families on common mapped genes,
   reported descriptively without treating genes as biological replicates.
7. Leave-one-study-family-out consensus direction: average standardized gene effects from the
   remaining favourable discovery families and report its Spearman correlation with the held-out
   family's effect vector. No intervention family may train a predictor that is tested on itself.

## Interpretation rules

- The original positive framework remains failed regardless of B1.
- Broad transportability is supported only if at least four of six methods orient both age controls
  correctly and the median pairwise method sign agreement across intervention contrasts is at least
  0.80. This threshold was selected before B1 outputs.
- A universal rejuvenation-direction claim additionally requires at least four of six methods to
  show positive family-mean effects in at least 80% of favourable-labelled families. This is a
  stringent descriptive gate, not a population hypothesis test.
- Failure supports a benchmark conclusion about non-portability of transcriptomic age directions,
  not a conclusion that source interventions lack functional or epigenetic benefit.
- No adaptive result establishes causality, safety, therapeutic efficacy, or biological-age change.
