---
title: "Supplementary information: Five chronological-age transcriptomic directions fail a locked portability benchmark across partial-reprogramming interventions"
author: "Oğuzcan Ünver"
date: "30 September 2026"
geometry: margin=1in
fontsize: 9pt
colorlinks: true
linkcolor: blue
urlcolor: blue
---

# Contents

1. Extended study design and temporal firewall  
2. Dataset and biological-unit audit  
3. Frozen framework and deviations  
4. Primary and secondary results  
5. Adaptive Benchmark B1  
6. Post-review Extension C1
7. Supplementary figures and tables
8. Reproducibility and rights boundary

# 1. Extended study design and temporal firewall

The discovery eligibility cutoff was 31 December 2024 at 23:59:59 UTC. Eligibility concerned the
exact data and annotations used, not merely an article date. The protocol, model definition,
validation endpoints, known-target list and analysis configuration were committed as tag
`protocol-v1.0` at commit `8be128cc5ab5f2cfc3147a0ee569fa59ff256b79`. The aggregate protocol
manifest SHA-256 is `c41ab8c013582925d1cd10c81513b4ee3c2a133da8eaafdc32bb8ea9f94c7d18`.

Nine discovery files (6.6 GB) were acquired and recorded in a manifest with aggregate SHA-256
`c7a5ddeadd7b15806313b82265b22e9508cd0fb83fe991f15a9058953d9a589a`. Protocol-v1
failed its Y-axis construction. Adaptive Amendment A1 was then defined and tagged at commit
`b7109948380cfd14acba5cad26f98b5acce7a897` before any intervention score was inspected.

The full discovery output, including all candidate ranks, was frozen as tag `discovery-a1-v1` at
commit `cdece0dec297800e7e3ab1aa1688f7f6d8a9d6b9`. Its 63-file release manifest has SHA-256
`73c7dc5c37d1b7f3170371d0bb7bb629ebbf0b70f07a8b1ac4a0aa266d29d0a7`. Only after this
release were target ranks inspected and validation matrices obtained. The conclusions of the source
papers and the five target names were known throughout; the design is not described as blinded or
prospective.

After the positive framework failed, Adaptive Benchmark B1 was specified and hashed before any
comparator model was fitted or applied. Its lock is tag `adaptive-benchmark-b1` at commit
`8ca3859e4118cbc86f6440069cf8dd697ab56e37`. The first complete benchmark run and every defect
found during audit were preserved. The corrected benchmark is committed at
`c05439da2e9bba4897b9fa312b46ad199afde2f0` [abbreviated in prose as `c05439d`].

# 2. Dataset and biological-unit audit

Supplementary Table 1 lists all registered datasets. GEO SuperSeries sample totals were never used as
replicate counts. Inference used the independent donor, animal, culture or independently prepared
experimental pool. Single cells and technical/culture replicates nested within the same source unit
were not counted as independent observations.

The two chronological-age references contain 133 and 82 healthy human donors. Ten Hutchinson–Gilford
progeria samples in GSE113957 were excluded from the normal-age reference. GSE297984 contains 24
RNA-seq profiles but only two aged donor lines (56 and 83 years); repeated cultures are nested within
donor, day and treatment. GSE297234 contains one 22-year and one 96-year donor in the analysed
trajectory and is descriptive. GSE300625 contains 20 mouse tissue profiles split between liver and
kidney. Exact sample-to-unit mappings are retained in `receipts/sample-metadata.tsv` and the analysis
code.

GSE276656 was preregistered as a secondary sensitivity but was not needed to adjudicate any primary
gate once the positive framework had failed. It was not added selectively to improve the benchmark.
Its omission is explicit rather than being presented as a missing result.

# 3. Frozen framework and deviations

## 3.1 Protocol-v1 failure

Protocol-v1 required concordant age-effect signs and Benjamini–Hochberg `q<=0.10` in each human age
reference. It selected zero genes. GSE226189's minimum adjusted value was approximately 0.132. The
rule and zero-gene output were preserved.

## 3.2 Adaptive Amendment A1

A1 used 10,532 common protein-coding genes, within-sample percentile ranks, and ridge regression on
within-cohort standardized age. The final penalty was 3.1622776601683795. Corrected outer-CV
Spearman correlations were −0.846594 in GSE113957, −0.457697 in GSE226189 and −0.730162 pooled.
Higher `Y_A1` was defined as more youth-associated by negating predicted standardized age. The score
has no unit and was not converted to years.

## 3.3 State axes

- `I` measures expression of a frozen somatic-identity programme, not complete cell function.
- `P` measures an endogenous late-pluripotency programme after removing the introduced Yamanaka
  symbols and vector features; it is not a tumorigenicity score.
- `D` averages response to oxidative stress, cellular response to DNA damage, endoplasmic-reticulum
  stress and cellular senescence programmes. It does not measure every form of damage.

The guarded-youthward call required all four axes to be available and their bootstrap intervals to
pass the frozen margins. An unavailable axis forced “not evaluable,” never “pass.”

## 3.4 Complete deviation register relevant to the manuscript

1. Protocol-v1 Y failed with zero genes. A1 was defined before intervention scores were viewed.
2. An early GSE176206 identity implementation admitted genes not detected across every control pool.
   It was repaired before the discovery release. All released discovery results use the corrected
   rule.
3. Two of five known targets were outside the immutable candidate universe. They were left missing,
   not replaced.
4. Frozen positive validation failed. No threshold, endpoint or model was changed to restore it.
5. B1 was defined after the primary outcomes were visible; it is adaptive and cannot rescue the
   original claim.
6. B1 r1 outer CV used full-cohort age-standardization moments. r1 was preserved, and r2 recomputed
   developmental CV using training-fold moments. Final models, intervention scores, contrasts,
   family summaries, gene-effect vectors and gates were unaffected.
7. B1's selected elastic-net model had zero non-zero coefficients and constant outputs. It was
   retained as an implementation failure and not replaced after results were visible.

The full dated register is `results-and-deviations.md`.

# 4. Primary and secondary results

## 4.1 Discovery-stage MPTR separation

Thirteen matched successful MPTR contrasts spanned three donors. Donor-level mean changes and
percentile-bootstrap intervals were:

| Axis | Mean change | 95% interval |
|---|---:|---:|
| Y_A1 | 0.43696 | 0.16285 to 0.61004 |
| Identity | -0.00356 | -0.00815 to 0.00017 |
| Endogenous pluripotency | 0.00256 | -0.00130 to 0.00507 |
| Aggregate measured damage | -0.00310 | -0.00392 to -0.00263 |

Twenty-one failed-to-reprogram matched contrasts were retained as a negative comparator, but they do
not constitute 21 independent donors. Their mean `ΔY_A1` was 0.06866.

![Discovery-stage donor-level changes for all four state axes. Dashed lines show the fixed ±0.05
identity/pluripotency margins.](figures/supplementary-figure-3-discovery-state-axes.png){width=75%}

## 4.2 Chemical discovery study

The untreated young-minus-old contrast in GSE246954 was positive (`ΔY_A1=0.24545`, interval 0.20114
to 0.29638). In contrast, old-cell 2c and 7c changes were −0.39583 and −0.46266; all four matched
replicates were negative for both cocktails. This discordance was visible before the temporal
validation release and already limited a universal interpretation.

## 4.3 Frozen temporal endpoints

In GSE297984, the eight donor-line/day/cocktail effects were:

| Donor age | Day | 2c ΔY_A1 | 7c ΔY_A1 |
|---:|---:|---:|---:|
| 56 | 6 | -0.2533 | -0.5561 |
| 56 | 14 | -0.5629 | -0.6454 |
| 83 | 6 | 0.0204 | -0.3880 |
| 83 | 14 | 0.1522 | -0.2202 |

Neither cocktail was positive across both donor lines with consistent day-level aggregate direction.
No population p-value was computed from the two donor lines.

In GSE297234, the young-minus-old day-0 control was `+0.85559`. Old-donor day changes relative to day
0 were:

| Day | ΔY_A1 | ΔI | ΔP | ΔD |
|---:|---:|---:|---:|---:|
| 3 | 0.30763 | -0.08229 | 0.09387 | 0.01129 |
| 7 | 0.56112 | -0.05797 | 0.09974 | 0.00311 |
| 10 | -0.91965 | -0.02668 | 0.04488 | 0.00368 |

Day 3 was simultaneously the first youthward point and the first point at which both fixed risk
margins were exceeded.

GSE300625 liver and kidney `ΔY_A1` values were −0.11742 (interval −0.19121 to −0.04130) and −0.07520
(−0.17313 to 0.02494). Tissue identity was not available, making the guarded-state endpoint
structurally non-evaluable. For GSE304043, GSTA4 `ΔY_A1` was 0.03244 (−0.63759 to 0.55661) and
aggregate `ΔD` was 0.00419 (−0.00089 to 0.00859).

## 4.4 Known-target rank

![Frozen percentile ranks for the five known later targets. Missing genes were not assigned a
replacement or artificial low rank.](figures/supplementary-figure-2-target-recovery.png){width=70%}

The evaluable-target mean percentile was 0.56091. In 100,000 expression-decile-matched three-gene
sets, the null mean was 0.48699 and its 95% interval was 0.17238 to 0.80345. The one-sided empirical
`P` value was 0.33237. This analysis is evaluable-only because two targets were missing; the complete
five-target endpoint failed by construction.

# 5. Adaptive Benchmark B1

## 5.1 Model panel and developmental performance

The five evaluable methods and pooled corrected outer-CV correlations were ridge-rank (−0.73016),
PCA50-ridge (−0.71304), meta-effect projection (−0.31406), GSE113957-only ridge (−0.61047), and
GSE226189-only ridge (−0.36476). All five correctly oriented both independent age controls. Final
hyperparameters and cohort-specific correlations are in Supplementary Table 3.

The corrected A1 predictions reproduce the independently generated A1 outer predictions to maximum
absolute difference `1.56e-15`. Fifteen repository tests passed at the B1-r2 checkpoint.

## 5.2 Complete contrast set

The complete method-by-contrast table is distributed as
`tables/supplementary-table-5-contrast-matrix.csv`. It contains all 184 declared contrasts,
including favourable-labelled, negative, trajectory, exploratory, adverse and age-control
categories. Source-unit labels are replaced by stable contrast codes in the redistribution-safe
table. No contrast was removed because its sign was unexpected.

![All 184 contrasts. Each method is divided by its own median absolute effect for visualization;
values are clipped to ±4. This scaling is visual only and is not used by any endpoint.](figures/supplementary-figure-1-full-contrast-matrix.png){width=80%}

Across the five evaluable methods, the pairwise sign-agreement values were:

|  | Ridge-rank | PCA50-ridge | Meta-effect | Cohort A | Cohort B |
|---|---:|---:|---:|---:|---:|
| Ridge-rank | 1.000 | 0.804 | 0.717 | 0.832 | 0.766 |
| PCA50-ridge | 0.804 | 1.000 | 0.783 | 0.842 | 0.701 |
| Meta-effect | 0.717 | 0.783 | 1.000 | 0.723 | 0.668 |
| Cohort A | 0.832 | 0.842 | 0.723 | 1.000 | 0.620 |
| Cohort B | 0.766 | 0.701 | 0.668 | 0.620 | 1.000 |

The upper-triangle median was 0.744565. Thus, the broad transportability gate failed despite all age
controls being oriented correctly.

## 5.3 Family-level direction

Family-mean directions were positive in 3/8 families for ridge-rank, 3/8 for PCA50-ridge, 5/8 for
meta-effect projection, 2/8 for cohort-A ridge and 3/8 for cohort-B ridge. Zero methods reached the
locked 80% threshold. The source-labelled family set is a concordance audit, not a sensitivity and
specificity analysis, because the source conclusions are heterogeneous and not a gold standard.

## 5.4 Gene-level effect concordance

There were 9,543 common genes in the complete family-effect matrix. The 28 off-diagonal family
correlations had median 0.053812. Leave-one-family-out correlations were:

| Held family | Spearman with seven-family consensus |
|---|---:|
| MPTR | 0.06081 |
| Mouse 2c/7c | 0.26497 |
| SOKM screen | 0.16757 |
| Human 2c/7c | 0.22642 |
| Human OSK/O4YRSK | 0.33626 |
| Human early OSKM | 0.23808 |
| RPE OSK | 0.13868 |
| RPE GSTA4 | 0.06008 |

The median was 0.196995. These are descriptive correlations over a fixed gene universe. Gene counts
were not used as biological-replicate sample sizes.

# 6. Post-review Extension C1

C1 was frozen at commit `93c04d7` before its outcomes were calculated. Its Markdown and JSON lock
files have SHA-256 values `7fdf1b4e27b0d16af3b0613cc011b0a91d8f5b80f9b9013430278b9a0166a9dd`
and `01420f7ebf81df459e022c7c87c49f300600da4f2ef1eabfaaae4bee1ec72ca0`. It was requested after
review of the complete manuscript, so every C1 result is explicitly outcome-aware and does not alter
the historical locks.

## 6.1 Hierarchical agreement and gate sensitivity

The contrast-weighted median pairwise agreement of 0.744565 was reproduced. Giving each of 14
families equal weight yielded 0.622514. The family-cluster bootstrap estimate had median 0.630315 and
95% percentile interval 0.515200–0.737703. Leaving out one family at a time produced estimates from
0.593477 to 0.652182. Thus no single family or contrast-rich study accounts for the below-gate
result.

The five model counts across the eight favourable-labelled families were 3, 3, 5, 2 and 3. One model
met a four-of-eight threshold and one met five-of-eight; none met six-, seven- or eight-of-eight.
Across agreement thresholds from 0.50 to 1.00, the observed contrast-weighted result passes only at
thresholds up to 0.725 and the equal-family result only up to 0.600. These curves are descriptive and
do not make any threshold biologically privileged.

## 6.2 Stratified summaries

| Dimension | Stratum | Families | Median agreement | Pairwise range |
|---|---|---:|---:|---:|
| Species | Human | 5 | 0.700 | 0.600–1.000 |
| Species | Mouse | 3 | 0.667 | 0.333–1.000 |
| Intervention | Chemical | 2 | 1.000 | 0.500–1.000 |
| Intervention | Genetic | 6 | 0.667 | 0.500–0.833 |
| Context | Fibroblast | 5 | 0.800 | 0.400–1.000 |
| Context | Retinal pigment epithelium | 2 | 0.500 | 0.000–1.000 |

Small strata and wide method-pair ranges preclude attribution of heterogeneity to one biological or
technical factor.

## 6.3 Gene-effect reliability

Only four favourable families had at least four contrast units. Their raw median split-half
Spearman correlations were 0.558 (MPTR), 0.840 (mouse chemical), 0.323 (SOKM) and 0.678 (human
chemical). Eight individual contrasts also supported independent arm splitting; raw correlations
ranged from 0.188 to 0.912. These values establish a partial, heterogeneous noise ceiling. They do
not support reliability correction for the four ineligible families, but show that the observed
cross-family median of 0.054 is well below within-context reproducibility in several estimable
settings.

## 6.4 Feature coverage

Historical-model coverage across ten expression matrices ranged from 0.8513 to 0.9999. Thresholds of
0.60, 0.70 and 0.80 retained all ten matrices, 184 contrasts and the same agreement statistics. A
0.90 threshold retained five matrices, 49 contrasts and eight total families; contrast-weighted and
equal-family median agreement remained 0.734694 and 0.666123. A separate refit on the 8,427 features
present in every matrix yielded 0.779891 and 0.648621.

## 6.5 Contemporary exploratory comparator

Pasta was pinned to official repository commit
`58bcc7a69ee97f2dc9e3623ac86538c251ddd498`. The three required model files had SHA-256 values:

- `v_genes_model.rda`: `42289278a0a0d5171af23b108b086deea7b413d287559a9ec39cf7c76269cb5f`
- `beta_Pasta.rda`: `b4a1891b7ed0737fbb784a89f6e00a9dfa07a15fa53359975f78f42fe92722e8`
- `cvfit_Pasta.rda`: `158a9ec9d9e53786c85bfc8af3593bfa99ae4bc949c4ad3f2afb767c14befa56`

The Python reconstruction reproduced the official three-sample example with maximum absolute error
`4.73×10^-11`. Pasta oriented the mouse and human age controls youthward by 18.195 and 52.184 units.
Its favourable-family effects were negative for MPTR, SOKM, RPE OSK and RPE GSTA4 and positive for
mouse chemical, human OSK/O4YRSK, early human OSKM and human chemical: four of eight positive. This
result is exploratory, post-review and outcome-aware. No exact benchmark accession was found in the
public article or packaged files, but the complete training sample list could not be audited.

The 2026 tAge model was not run. Its official public repository uses the MGB Open Access License 1.0
for non-commercial academic use and directs commercial users to obtain a separate agreement. The
author's commercial affiliation made exclusion the conservative rights-compliant choice.

# 7. Supplementary tables

- **Supplementary Table 1:** `tables/supplementary-table-1-datasets.csv` — dataset roles, species,
  modality, biological unit and public date.
- **Supplementary Table 2:** `tables/supplementary-table-2-primary-scorecard.csv` — complete frozen
  endpoint scorecard.
- **Supplementary Table 3:** `tables/supplementary-table-3-benchmark-summary.csv` — model
  hyperparameters, CV correlations, control orientation and family fractions.
- **Supplementary Table 4:** `tables/supplementary-table-4-known-targets.csv` — target
  evaluability and ranks.
- **Supplementary Table 5:** `tables/supplementary-table-5-contrast-matrix.csv` — all 184
  contrasts and model effects, with source-unit labels replaced by stable codes.
- **Supplementary Table 6:** `results/benchmark-b1/favourable-family-means.csv` — family means.
- **Supplementary Table 7:** `results/benchmark-b1/family-gene-effect-spearman.csv` — complete
  gene-effect correlation matrix.
- **Supplementary Table 8:** `results/benchmark-b1/leave-one-family-out-gene-consensus.csv` —
  consensus transfer results.
- **Supplementary Tables 9–22:** `results/post-review-c1/` — hierarchical agreement, cluster
  bootstrap, threshold sensitivity, strata, reliability, feature coverage, common-feature refit and
  Pasta outputs.

# 8. Reproducibility and rights boundary

The reproducibility package includes analysis and figure source code, configurations, exact input
manifests, tests, aggregate result tables and audit receipts. It excludes raw public expression
matrices, source H5AD/H5 files, participant-linked records and source content not licensed for
redistribution. Users obtain raw data from the original repositories and verify them against the
provided hashes.

The first B1 run and its corrected CV revision remain separately available. Historical outputs are
never overwritten. The main manuscript uses `outer-cv-correlations-r2.csv` and `summary-r2.json`.
Intervention results are the original B1 results because the repaired target scaling did not enter
final model fitting or application.

The analysis package makes no claim of causal rejuvenation, direct biological age, clinical
validity, therapeutic efficacy, safety, peptide efficacy or target causality. It also does not imply
that transcriptomic non-concordance invalidates orthogonal functional or epigenetic findings in the
source studies.
