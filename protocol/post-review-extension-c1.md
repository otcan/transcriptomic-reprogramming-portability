# Post-review extension C1: clustered portability, reliability, coverage and contemporary comparator

**Locked before C1 outcome computation: 30 September 2026.** This amendment responds to author
review of the immutable Paper 2 r2 candidate at baseline commit
`6bb77f63f813e4160f43f1b490d31a3e5694c9f0`. It is outcome-aware with respect to all previously
reported A1/B1 results. Every C1 result is therefore a labelled post-review sensitivity or
exploratory analysis, not a new prospective test and not a revision of the historical locks.

## Fixed claim and title boundary

The revised title is:

> **Five chronological-age transcriptomic directions fail a locked portability benchmark across
> partial-reprogramming interventions**

The paper may conclude that the five evaluated directions recognize chronological-age controls yet
do not show consistent direction under the evaluated interventions. It may recommend that age
prediction alone is insufficient qualification for an intervention biomarker. It may not conclude
that all transcriptomic clocks fail, that chronological clocks cannot measure rejuvenation, or that
the source interventions lack beneficial biological effects.

Target recovery remains secondary corroboration because only three of five named targets were
evaluable. The 184 contrasts are clustered observations and must never be described as 184
independent experiments.

## Historical analysis preservation

- Protocol-v1, A1 and B1 files and accepted outputs remain byte-preserved.
- The degenerate elastic-net attempt remains reported but excluded from all five-method C1
  summaries because it is not an evaluable direction.
- C1 uses the five evaluable B1 directions without refitting except in the explicitly labelled
  common-feature sensitivity.
- Existing favourable-family labels remain source-author interpretations, not rejuvenation truth.

## C1-A: family-aware agreement and uncertainty

The five evaluable methods are `ridge_rank_a1`, `pca50_ridge`, `meta_effect_projection`,
`gse113957_only_ridge` and `gse226189_only_ridge`.

1. Reproduce the contrast-weighted median of the ten pairwise sign-agreement proportions over all
   184 contrasts.
2. For each method pair and benchmark family, compute within-family sign agreement. Average those
   proportions equally across families, then take the median over the ten method pairs. This
   equal-family statistic is the primary C1 clustered sensitivity.
3. Cluster-bootstrap families with replacement for 10,000 draws (seed `3092026`). Each selected
   family contributes its within-family agreement once; families, not contrasts, are the resampling
   units. Report the percentile 95% interval for the median pairwise equal-family statistic.
4. Repeat the equal-family statistic after omitting each family in turn. Report the complete range;
   do not select a favourable omission.
5. Show the transportability decision over agreement thresholds from 0.50 to 1.00 in steps of 0.025
   for both contrast-weighted and equal-family statistics. The earlier 0.80 rule remains a declared
   engineering criterion, not a natural biological boundary.
6. For the eight favourable families, show each method's exact positive count out of eight and the
   number of methods meeting thresholds 4/8 through 8/8. The existing 80% rule corresponds to at
   least 7/8.

## C1-B: prespecified descriptive strata

Use favourable families only and report exact family counts, signs and equal-family method
agreement for strata with at least two families:

- species: human or mouse;
- intervention: genetic or chemical;
- context: fibroblast, mesenchymal/adipogenic screen, or retinal pigment epithelium.

The fixed family-to-stratum mapping is stored in `post-review-extension-c1.json`. Strata are
descriptive and overlapping; no population inference or multiple-testing claim is allowed.

## C1-C: gene-effect reliability / noise ceiling

The original cross-family gene-effect analysis remains unchanged. C1 adds two reliability checks:

1. **Family contrast-unit split halves.** For a favourable family with at least four contrast-level
   gene-effect vectors, make 2,000 seeded balanced random splits of contrast units. Correlate the two
   half-mean vectors by Spearman correlation and report the median raw and Spearman--Brown-adjusted
   value. Seed: `3092027`.
2. **Within-contrast arm splits.** Where both treated and control arms contain at least two biological
   units, enumerate balanced non-empty arm splits when feasible, otherwise sample up to 2,000 unique
   splits with seed `3092028`. Compute independent treatment-minus-control vectors for the two
   halves, their Spearman correlation and Spearman--Brown adjustment.

Report coverage and inestimable families/contrasts. Compare cross-family similarity with within-
context reliability only where estimable; do not call the result a universal noise ceiling.

## C1-D: feature coverage and missing-feature sensitivity

For every expression key, report observed B1 features, total B1 features and coverage. Recompute the
five-method summaries while retaining only expression keys meeting thresholds 0.60, 0.70, 0.80 and
0.90. Thresholding only changes eligibility; it does not change a score already computed with the
frozen training-mean fill.

As a secondary sensitivity, form the intersection of B1 genes observed across every benchmark
expression key. If at least 2,000 genes remain, refit the same five methods on the two age-reference
cohorts using only that intersection and repeat scoring. If fewer remain, report the prespecified
feasibility failure and do not relax the threshold after seeing the count.

## C1-E: contemporary exploratory comparator

Run Pasta as one explicitly contemporary, non-historical comparator:

- repository: `https://github.com/jsalignon/pasta`;
- pinned commit: `58bcc7a69ee97f2dc9e3623ac86538c251ddd498`;
- licence: MIT;
- `v_genes_model.rda` SHA-256:
  `42289278a0a0d5171af23b108b086deea7b413d287559a9ec39cf7c76269cb5f`;
- `beta_Pasta.rda` SHA-256:
  `b4a1891b7ed0737fbb784a89f6e00a9dfa07a15fa53359975f78f42fe92722e8`;
- `cvfit_Pasta.rda` SHA-256:
  `158a9ec9d9e53786c85bfc8af3593bfa99ae4bc949c4ad3f2afb767c14befa56`.

Reproduce the package's published preprocessing: map human symbols to versionless human Ensembl
identifiers; use the fixed 8,113 genes; fill missing rows with the median of observed matrix values;
rank each sample across genes with average ties; apply the pretrained `lambda.min` coefficients;
multiply the link prediction by the supplied Pasta beta; and reverse sign only for display so higher
means more youth-associated. Validate the Python reader/predictor against the package's bundled
example before benchmark interpretation.

Report Pasta feature coverage, both age-control directions, all contrast effects and the eight
favourable-family signs. Pasta is not added retroactively to the B1 gate. Audit exact training and
evaluation accession overlap; any unresolved overlap must be disclosed. Whether Pasta agrees or
disagrees, report the result without changing endpoints.

The Gladyshev-Lab `tAge` package is not executed because its official MGB Open Access License 1.0
limits use to non-commercial academic purposes and requires a separate agreement for commercial
use. This exclusion is a rights boundary, not an adverse scientific judgment.

## C1-F: figures, accessibility and release

- Replace red/green-only semantics with a colourblind-safe blue/orange/grey palette plus signs,
  shapes or direct labels.
- Add a clustered-sensitivity panel and a coverage/reliability/contemporary summary without hiding
  historical failures.
- Regenerate every figure and table from committed source data.
- Update the claim verifier and tests before manuscript packaging.
- Keep the repository private and the release draft until renewed author approval. Public release,
  DOI registration and any preprint-server post remain separate irreversible actions.

