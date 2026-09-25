# Adaptive Amendment A1 — regularized youthward axis

Date defined: 2026-09-25
Status: to be locked before intervention-state outcomes or any validation matrix is inspected

## Trigger and disclosure

The protocol-v1.0 Y feature rule required BH q<=0.10 in each of GSE113957 and GSE226189. It selected
zero genes because GSE226189 had no gene meeting that cohort-wide FDR threshold. Identifier,
metadata, model, and FDR calculations were audited and no implementation defect was found. The
protocol-defined Y axis therefore fails and remains reported as a failed primary construction.

Amendment A1 is a discovery-stage model-development response. Its held-out age performance has
already been examined during development and is not confirmatory. It is nevertheless eligible for
later temporal evaluation because its complete model is frozen before intervention-state outcomes
are inspected and before any post-2024 expression matrix is opened.

## Exact A1 model

1. Use the same filtered GSE113957 and GSE226189 healthy-donor expression matrices and covariate-free
   sample inclusion as protocol v1.0.
2. Feature universe: protein-coding gene symbols present after filtering in both references,
   according to Ensembl release 111.
3. Convert expression to within-sample percentile ranks across the complete filtered matrix, then
   retain the common feature universe.
4. Standardize chronological age to mean zero and sample standard deviation one separately within
   each cohort. Concatenate cohorts.
5. Fit ridge regression with an intercept. Choose alpha by generalized leave-one-out cross-validation
   from `10**[-2,-1.5,...,6]` (17 log-spaced values). No validation or intervention data enter alpha
   selection.
6. Report one outer five-fold donor cross-validation using the unchanged v1.0 fold assignments.
   Refit alpha and all coefficients inside each training fold.
7. Fit the final model on all reference donors and export ordered features, training feature means,
   intercept, coefficients, alpha, source/input/code hashes, and outer predictions.
8. For a new sample, rank all observed genes within that sample, select model features, and impute a
   missing feature with its frozen training mean. At least 60% of model features must be observed.
9. Define `Y_A1 = -ridge_prediction`. Higher values are therefore more youth-associated. This is a
   relative transcriptomic age-direction score, not years of age or a biological-age clock.

## Evaluation and claims

- The v1.0 Y gate remains failed and is never overwritten.
- A1 age cross-validation is model-development evidence and is labelled adaptive.
- I, P, D, preservation margins, biological units, intervention contrasts, candidate rank, known
  target list, adverse challenge, and temporal validation datasets remain unchanged.
- The temporal endpoints use `Y_A1` in place of unavailable `Y_v1`. This substitution is disclosed
  in every primary table and the manuscript.
- No new hyperparameter or threshold may be chosen after the A1 tag.

