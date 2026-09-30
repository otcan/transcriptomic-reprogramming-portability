---
title: "Five chronological-age transcriptomic directions fail a locked portability benchmark across partial-reprogramming interventions"
author:
  - "Oğuzcan Ünver"
date: "30 September 2026"
geometry: margin=1in
fontsize: 10pt
colorlinks: true
linkcolor: blue
urlcolor: blue
---

**Article type:** Analysis<br>
**Target journal:** *Nature Aging*<br>
**Affiliation:** Metastate Bio Inc, Wilmington, Delaware, USA<br>
**Corresponding author:** Oğuzcan Ünver; can@metastate.bio; ORCID:
0009-0007-2023-5084

# Abstract

Chronological-age models are used as scalar measures of cellular rejuvenation, although age
prediction does not establish sensitivity to beneficial intervention. We constructed a temporally
locked four-axis framework from public data available by 2024, challenged it against later studies
without refitting, and compared five evaluable age directions across 184 partial-reprogramming
contrasts. The frozen positive framework failed its human chemical and trajectory-order endpoints.
All five benchmark models oriented two independent young–old controls correctly, but their median
pairwise agreement on intervention direction was 0.745. Equal weighting across 14 study families
reduced agreement to 0.623 (family-cluster bootstrap 95% interval 0.515–0.738). No model was
youthward in more than five of eight source-labelled favourable families. Cross-family gene effects
were weakly concordant relative to the within-family reproducibility estimable in four families. A
post-review exploratory application of the contemporary Pasta model oriented both age controls but
was positive in only four of eight favourable families. These results show that age prediction did
not establish intervention portability for the evaluated directions. Rejuvenation claims require
intervention-specific and orthogonal qualification.

Molecular clocks compress high-dimensional measurements into quantities associated with age,
morbidity, mortality or the pace of ageing [1–3]. This compression is attractive for partial
reprogramming, where transcription-factor or chemical interventions can induce broad and rapid state
changes. A single scalar appears to offer a direct answer to a difficult question: did an old cell
become younger? Yet chronological age is a training label, not an intervention outcome. A predictor
can discriminate young from old samples while responding to perturbation-specific expression,
changes in composition, stress, proliferation or loss of identity in ways that do not generalize
across experiments.

This distinction is especially important for partial reprogramming. Transient expression of
pluripotency factors has been reported to shift transcriptomic and epigenetic features towards
younger states while retaining or reacquiring somatic identity [4–6]. Single-cell analyses also show
that youthful expression can coincide with transient suppression of identity and activation of
pluripotency programmes [7]. Chemical approaches produce extensive metabolic and transcriptional
changes and clock responses [8,9], whereas later in-vivo work has exposed toxicity and an absence of
detectable transcriptomic-age improvement in liver and kidney [10]. Later human trajectories and
retinal studies add further experimental contexts [11,12]. These studies need not be mutually
inconsistent: they interrogate different species, cell types, interventions, schedules, modalities
and functional endpoints. They do, however, create a stringent test of whether a transcriptomic
age-associated direction is a portable scalar measure of intervention-induced rejuvenation.

We therefore ran a failure-preserving computational study. First, we specified a four-axis state
vector before opening the later expression matrices: a youth-associated chronological-age direction
(`Y`), retained somatic identity (`I`), endogenous pluripotency or dedifferentiation (`P`), and four
measured stress/damage programmes (`D`). Second, we froze discovery outputs before application to
datasets released in 2025–2026. Publication conclusions were known, so the exercise was temporally
external and outcome-informed rather than blinded or prospective. Third, after the positive framework
failed, we locked an explicitly adaptive multi-model benchmark before computing comparator outputs.
This design distinguishes prespecified failures from post-failure characterization and retains every
unfavourable result (Fig. 1).

![Study design and qualification logic. **a,** Discovery, locking, temporal validation and the
post-primary benchmark. **b,** The four frozen state axes. **c,** The qualification ladder used for
interpretation.](figures/figure-1-study-design.png){width=100%}

## A frozen multidimensional framework fails its positive validation

The original protocol required genes to show concordant age effects at a false-discovery rate of
0.10 in each of two human fibroblast cohorts. No gene met that rule because the smaller cohort's
minimum adjusted value was 0.132. We preserved that failure and, before viewing any intervention
outcome, defined Adaptive Amendment A1: within-sample percentile ranks for 10,532 common
protein-coding genes and ridge regression against within-cohort standardized chronological age. The
negative of the predicted standardized age was called `Y_A1`, such that higher scores were more
youth-associated. In five-fold outer cross-validation, Spearman correlations between `Y_A1` and age
were −0.847 in GSE113957, −0.458 in GSE226189 and −0.730 pooled (Fig. 2a). These values establish an
age association; they do not validate `Y_A1` as biological age.

In discovery data, maturation-phase transient reprogramming (MPTR) provided the intended separation.
Across three donors, the donor-mean change in `Y_A1` was 0.437 (95% percentile-bootstrap interval
0.163 to 0.610), while mean changes in identity (−0.0036), endogenous pluripotency (0.0026) and
aggregate measured damage (−0.0031) remained within the locked guardrails. However, chemical
reprogramming already exposed context dependence: the untreated young-minus-old control was
correctly oriented (`ΔY_A1=0.245`, interval 0.201 to 0.296), whereas old-mouse 2c and 7c effects were
−0.396 and −0.463, respectively.

We hashed and committed the complete discovery release before obtaining the temporal-validation
matrices. In the two-donor human chemical study, neither 2c nor 7c was positive across both donor
lines and both assessed days; six of eight donor–day–cocktail effects were negative (Fig. 2c). In the
96-year donor OSKM trajectory, `Y_A1` first increased at day 3 (`Δ=0.308`), but identity loss
(`ΔI=−0.082`) and pluripotency gain (`ΔP=0.094`) crossed their locked margins on the same day
(Fig. 2d). Thus, youth-associated movement did not precede the risk axes. The 22-year minus 96-year
day-0 control remained correctly oriented (`ΔY_A1=0.856`).

The adverse in-vivo 7c challenge yielded negative mean `Y_A1` changes in liver (−0.117; 95%
interval −0.191 to −0.041) and kidney (−0.075; interval −0.173 to 0.025). It formally avoided a
guarded-youthward call, but this endpoint was non-informative because a tissue-specific identity axis
was unavailable. GSTA4 overexpression in aged retinal pigment epithelium produced `ΔY_A1=0.032`
(interval −0.638 to 0.557) and did not lower aggregate `D`. These expression results neither validate
the frozen score nor negate functional findings in the source studies.

Finally, a complete candidate list had been frozen before five known later targets were looked up.
GSTA4 and TOMM70A were absent from the immutable 9,131-gene intersection universe. DDX21, SERBP1 and
PHGDH ranked at percentiles 0.793, 0.607 and 0.282. Their evaluable-only mean percentile was 0.561,
with an expression-decile-matched empirical one-sided `P=0.332` over 100,000 draws. The target
endpoint was therefore negative and incomplete. Collectively, the frozen positive framework failed
its chemical-validation and trajectory-order gates and did not recover the later target set
(Fig. 2b).

![Frozen framework results. **a,** Held-out age association in the two reference cohorts. **b,**
Complete primary scorecard. **c,** Human 2c/7c effects by donor line and day. **d,** Old-donor OSKM
trajectory for the youth, identity and endogenous pluripotency axes.](figures/figure-2-frozen-framework.png){width=100%}

## Age models orient controls but disagree under intervention

One failed age direction cannot establish a general limitation. After closing the primary claim, we
therefore specified Adaptive Benchmark B1 and fixed its models, contrasts, family units and
descriptive gates before fitting any comparator. Six approaches were attempted: the existing
ridge-rank model; PCA followed by ridge regression; a sign-concordant cross-cohort age-effect
projection; ridge models fitted to each reference cohort separately; and elastic net. The selected
elastic-net model had zero non-zero coefficients and constant outputs. We retained this failure and
did not choose a replacement after seeing intervention results. Five methods remained evaluable.

All five models associated with chronological age in corrected outer cross-validation (pooled
Spearman correlations −0.730 to −0.314; Fig. 3a). Moreover, all five correctly oriented two
young-minus-old controls, one in mouse fibroblasts and one in human fibroblasts (Fig. 3b). These
checks matter: the subsequent disagreement cannot be attributed simply to arbitrary reversal of
score direction.

We next applied the fixed final models to 184 matched intervention contrasts spanning MPTR,
Yamanaka-factor combinations, 2c/7c chemical treatments, human OSK/O4YRSK and OSKM trajectories,
retinal OSK, GSTA4 and the adverse in-vivo study. For the locked broad-transportability gate, at least
four of six attempted methods had to orient both age controls and median pairwise sign agreement
across intervention contrasts had to reach 0.80. Restricting the summary to the five evaluable
methods, median sign agreement was 0.745. The result remained below threshold; including the failed
constant method would reduce interpretability rather than provide evidence for transportability.

Because the 184 contrasts are clustered within studies and families, we added an explicitly
post-review, outcome-aware sensitivity under locked Amendment C1. Equal weighting of the 14
benchmark families lowered median pairwise agreement to 0.623. A 10,000-draw cluster bootstrap that
resampled families gave a percentile interval of 0.515–0.738, and leave-one-family-out estimates
ranged from 0.593 to 0.652 (Fig. 5a,b). The original contrast-weighted statistic is therefore not
presented as 184 independent experiments; the family-weighted result is the more conservative
summary.

The divergence was visible at the study-family level (Fig. 3c). MPTR was positive for every method,
whereas the human chemical family was negative for every method. Other families split by model. For
example, the meta-effect and cohort-B directions were positive for human OSK/O4YRSK while the other
three models were negative. GSTA4 was positive under three directions and negative under two. We had
defined a stringent universal-direction gate: at least four methods had to be positive in at least
80% of eight source-labelled favourable families. With eight families, that rule requires at least
seven positive calls. The observed counts were two to five; only one method reached four or five,
and none reached six, seven or eight. Thus the conclusion does not depend on treating 80% as a
biologically privileged boundary (Fig. 3d). The labels describe favourable interpretations in source
papers and are not ground-truth rejuvenation outcomes.

![Adaptive multi-model benchmark. **a,** Corrected outer-CV age associations for five evaluable
models. **b,** Independent young-minus-old controls. **c,** Signs of family-mean intervention
effects. **d,** Fraction of source-labelled favourable families with positive mean effects.](figures/figure-3-model-benchmark.png){width=100%}

Descriptive stratification did not identify a broadly consistent subset. Median method agreement was
0.700 across five human families, 0.667 across three mouse families and 0.667 across six genetic
families. The five fibroblast families had a median of 0.800, but individual method pairs ranged from
0.400 to 1.000. The chemical and retinal strata each contained only two families and are too small
for general conclusions. These analyses localize heterogeneity but do not attribute it to species,
tissue or intervention class.

## Intervention effects do not converge on a common gene direction

Model disagreement might arise even if interventions share a common gene-level programme that is
weighted differently. We therefore compared family-mean expression-effect vectors on 9,543 common
mapped genes. Pairwise Spearman correlations were low overall (median 0.054), ranging from −0.085 to
0.517 (Fig. 4b). The largest value linked the two related human OSK/OSKM families. This local
agreement did not imply a universal programme across genetic, chemical, species and retinal
contexts.

We also formed a leave-one-family-out consensus by averaging standardized expression effects from
the other seven favourable-labelled families and correlating it with the held family. Correlations
ranged from 0.060 to 0.336, with a median of 0.197 (Fig. 4c). Genes were not treated as biological
replicates and these correlations are descriptive; no gene-count-based population inference is
claimed.

We estimated a partial noise ceiling wherever the design permitted split-half comparisons. Four of
eight favourable families had at least four contrast units; their median raw split-half correlations
were 0.323, 0.558, 0.678 and 0.840. Eight individual contrasts supported arm-split estimates, ranging
from 0.188 to 0.912. Reliability is therefore heterogeneous and unavailable for half the families,
but in the estimable higher-reliability contexts it substantially exceeds the cross-family median of
0.054. The weak cross-family concordance cannot be assigned entirely to measurement noise, although
the incomplete noise ceiling prevents a universal reliability-corrected claim. Together with the
method sign matrix (Fig. 4a), these results support non-convergence in the evaluated data at both the
scalar-score and gene-effect levels.

![Cross-method and cross-family concordance. **a,** Pairwise method sign agreement across 184
contrasts. **b,** Spearman correlations between family-mean gene effects. **c,** Leave-one-family-out
consensus correlations. **d,** Supported and unsupported interpretation.](figures/figure-4-concordance-and-interpretation.png){width=100%}

## Robustness and a contemporary exploratory comparator

Observed coverage of the 10,532 historical-model features ranged from 85.1% to 100.0%, rather than
approaching the allowed 60% floor. Raising the eligibility threshold to 70% or 80% retained every
dataset and result. At 90%, five expression matrices and 49 contrasts remained; contrast-weighted
and equal-family agreements were 0.735 and 0.666. As a separate stress test, refitting all five
historical methods on the 8,427 features present across every benchmark matrix yielded agreements of
0.780 and 0.649, still below 0.80 (Fig. 5c).

We also applied Pasta, a contemporary rank-based age-shift model released after the historical
cutoff, as a post-review exploratory comparator [14]. A Python reconstruction reproduced the three
bundled reference predictions to a maximum absolute error of `4.8×10^-11`. Pasta oriented both
young-minus-old controls correctly but was youthward in four of the eight favourable-labelled
families (Fig. 5d). Its feature coverage was 83.3–99.5%. This analysis was not part of B1, was added
after the historical outcomes were known, and does not constitute an independent prospective test.
No exact benchmark accession overlap was identified in the public article or model package, but a
complete sample-level audit of all 21 Pasta training datasets was not possible from the released
materials. We excluded the contemporary tAge model from this commercial-affiliation release because
its public licence is restricted to non-commercial academic use [13].

![Post-review robustness analyses. **a,** Contrast-weighted and equal-family method agreement; the
point and whisker show the family-cluster bootstrap median and 95% interval. **b,** Leave-one-family-
out agreement. **c,** Feature-coverage thresholds and the common-feature refit. **d,** Historical
models and exploratory Pasta signs across eight source-labelled favourable families. Blue denotes a
positive youth-oriented effect and orange a negative effect; signs provide redundant encoding.](figures/figure-5-post-review-robustness.png){width=100%}

# Discussion

Our central result is a qualification failure with a useful boundary: five transcriptomic
chronological-age directions recognized young–old contrasts yet were not sufficiently portable
across the evaluated partial-reprogramming interventions. This is not evidence that the source
interventions lack biological or functional effects. It shows that age prediction and intervention
measurement are different validation problems.

Recent large-scale work reaches a compatible conclusion from a different direction. Universal
transcriptomic clocks trained on more than 11,000 samples found that chronological clocks correlated
well with age but poorly with lifespan-modulating intervention effects, while mortality-oriented
objectives performed better [13]. Pasta demonstrates that rank-based, multi-tissue transcriptomic
models can generalize across platforms and recover experimental perturbations [14]. Our exploratory
Pasta result shows that a stronger contemporary chronological model does not automatically resolve
the present portability problem, but one model cannot represent the full modern clock class. The
historical result and contemporary extension remain analytically separate. Our contribution is a
partial-reprogramming-specific stress test with temporally separated data, identity and pluripotency
axes, an adverse challenge, preserved failures and complete method-by-contrast disclosure.

The findings argue against treating “transcriptomically younger” as a self-sufficient rejuvenation
claim. At minimum, a candidate score should demonstrate held-out age association, correct control
orientation, portability to interventions that were not used to define it, and concordance with an
outcome-relevant endpoint. Identity, pluripotency, stress, toxicity and function should remain
separate measurements rather than being hidden inside a single scalar. This distinction also avoids
the inverse error: a beneficial intervention need not reverse a chronological-age clock. For
example, elamipretide improved cardiac and skeletal-muscle function in aged mice without detectable
changes in tissue epigenetic or transcriptomic age [15].

Several limitations constrain the generality of our benchmark. The age references were cultured
human fibroblasts, whereas interventions spanned species, cell types, platforms and tissues. That
heterogeneity is the intended portability challenge, but it also prevents attribution of failure to
one cause. Some validation studies had only one or two donor lines, and several contrasts are
descriptive. Source-labelled favourable families are not a gold standard. Within-sample ranks reduce
scale dependence but cannot remove composition and state confounding. The four `D` programmes cover
selected transcriptional responses rather than total damage, and `P` is not a tumorigenicity score.
The benchmark comparator was adaptive after primary failure, and its numerical gates are transparent
engineering criteria rather than natural biological thresholds or population hypothesis tests. The
post-review analyses were outcome-aware, Pasta's complete training-sample overlap could not be
audited from public materials, and the tAge licence did not permit inclusion under the present
commercial affiliation. Newly published mortality- and outcome-trained clocks should therefore be
evaluated in a prospectively locked future benchmark with compatible inputs and rights.

The practical implication is methodological. A rejuvenation biomarker should be qualified for its
intended context, not promoted from chronological-age prediction by analogy. Public reprogramming
data now permit a tiered evaluation in which construction, intervention transport, state preservation
and orthogonal function are separately auditable. In the datasets assessed here, none of the five
historical directions met the locked portability criteria, and the exploratory Pasta model remained
split across favourable-labelled families. This does not establish that chronological transcriptomic
clocks can never measure an outcome-relevant component of rejuvenation.

# Methods

## Study design and temporal firewall

The discovery cutoff was 31 December 2024. Discovery comprised human fibroblast age references
GSE113957 (133 healthy donors, age 1–94 years) and GSE226189 (82 healthy donors, age 22–89 years),
GSE165176/GSE165177 (human Sendai and MPTR trajectories), GSE176206 (mouse single-cell factor
screens), and GSE246954/GSE247199 (mouse chemical reprogramming). Temporally external evaluation
used GSE297984, GSE297233/GSE297234, GSE300625 and GSE304042/GSE304043/GSE304044. Exact source
URLs, retrieval times, file sizes and SHA-256 hashes are recorded in the input manifests.

The protocol, code configuration and named-target list were committed and tagged before discovery
expression matrices were opened. Protocol-v1 failed, after which A1 was defined and frozen before any
intervention score was examined. Discovery results and the full 9,131-gene candidate rank were then
hashed and committed before validation matrices were obtained. B1 was specified only after primary
validation failed and is labelled adaptive throughout. Published study conclusions and named targets
were known; “temporally external” does not mean prospective or blinded.

Reviewer-requested Extension C1 was locked before its outcomes were computed. It prespecified family
weighting, cluster bootstrap, leave-one-family-out, threshold, stratum, reliability and feature-
coverage analyses, plus one contemporary comparator. C1 is explicitly post-review and outcome-aware;
it tests robustness but does not retroactively alter or rescue B1.

## Expression representation and orthology

Gene identifiers were reduced to approved symbols using Ensembl release 111 annotations. Duplicate
symbols were collapsed by median expression, except raw counts mapped before normalization were
summed. Cross-species analyses used NCBI HomoloGene build 68 one-to-one human–mouse mappings; groups
with ambiguous or one-to-many mappings were excluded. Counts were converted to log2 counts per
million plus 0.5 after source-level filtering. Deposited normalized expression was used when a common
raw matrix was unavailable. Within each sample, eligible genes were percentile-ranked. New datasets
required at least 60% feature coverage; unobserved features were filled with frozen training-feature
means. C1 reported actual coverage for every matrix, repeated eligibility at 70%, 80% and 90%, and
refitted the historical models on the 8,427-feature intersection shared by all benchmark matrices.

## Construction of Y_A1

We retained 10,532 protein-coding genes common to both age references. Chronological age was
standardized within cohort using training-fold means and standard deviations. Ridge regression was
fitted to within-sample rank features, with the penalty selected by generalized cross-validation from
17 values spanning `10^-2` to `10^6`. Five outer folds were assigned within study, sex and age
quintile with seed 1729. Feature transformation, age standardization and fitting were repeated inside
each training fold. `Y_A1` was the negative predicted standardized age so that higher values were
youth-associated. Final model coefficients were fitted on both complete references. Scores are
dimensionless and were never converted to years.

## Identity, pluripotency and damage axes

The fibroblast identity and endogenous pluripotency programmes were derived from six donor lines in
GSE165176 by paired baseline-to-late-state consistency rules. Identity genes decreased by at least
one deposited log2-expression unit in at least five of six lines; pluripotency genes increased by the
same criterion. POU5F1, SOX2, KLF4, MYC, vector features and non-standard identifiers were excluded
from the primary pluripotency programme. Source-specific adipogenic and mesenchymal identity
programmes were used only within GSE176206.

Damage components were Gene Ontology descendants of response to oxidative stress (GO:0006979),
cellular response to DNA damage (GO:0006974), response to endoplasmic-reticulum stress (GO:0034976)
and cellular senescence (GO:0090398), using the 17 January 2024 ontology and matching GOA/MGI
annotations [16]. `NOT` and `ND` annotations were excluded. Each one-sided programme score was mean
within-sample percentile rank minus 0.5; aggregate `D` was the unweighted mean of all four available
components. Programme coverage required at least 15 genes and 50% of the frozen set.

## Biological units and intervention contrasts

The donor, animal, culture or independently prepared experimental pool was the unit. Cells and
technical replicates were nested and never counted as independent replicates. Treatment effects were
treatment minus matched contemporaneous control. Primary intervals were percentile bootstrap
intervals over biological units (10,000 draws, seed 1729) when each arm contained at least three
units; smaller designs were descriptive. For GSE297984, culture replicates were collapsed within the
56- and 83-year donor lines, day and treatment. For GSE176206, cells were pseudobulked by source
experimental pool and factor combination.

The frozen guarded-youthward rule required the 95% interval for `ΔY` to be above zero, the `ΔI`
interval above −0.05, the `ΔP` interval below +0.05 and the `ΔD` interval below zero. If any axis
was unavailable, guarded status was not evaluable. The human chemical gate required one cocktail to
be positive in both donor lines with the same aggregate direction at both days. The trajectory gate
required youthward movement before either identity loss below −0.05 or pluripotency gain above
+0.05.

## Candidate ranking and matched null

Candidates were one-to-one protein-coding orthologues detected in both age references, MPTR and the
mouse chemical study. The frozen score combined age-effect magnitude, reversal of age direction in
successful discovery contrasts, and reversal consistency relative to negative controls with weights
0.40, 0.40 and 0.20. The five known later targets were GSTA4, DDX21, TOMM70A, SERBP1 and PHGDH [12,17].
Their mean frozen-rank percentile was compared with 100,000 gene sets matched to discovery-reference
expression deciles (seed 271828). Missing targets were reported and not replaced.
This underpowered three-of-five evaluable analysis was treated as secondary corroboration, not as a
major basis for the portability conclusion.

## Adaptive Benchmark B1

B1 used the same 10,532 ranked features and age references. Its six attempted methods were:

1. A1 ridge-rank;
2. elastic-net rank regression over `l1_ratio={0.1,0.5,0.9}` and 25 penalties from `10^-4` to `10^1`;
3. 50-component PCA followed by ridge regression;
4. a projection on the mean of cohort-specific standardized age effects, with weights set to zero
   when cohort signs differed;
5. ridge-rank fitted only to GSE113957; and
6. ridge-rank fitted only to GSE226189.

The final elastic-net model selected penalty 10 and `l1_ratio=0.1`, producing zero non-zero
coefficients. It was retained as an implementation failure and excluded from evaluable-method
summaries without replacement. Corrected outer CV standardized chronological age using training-fold
moments. Final intervention scores, contrasts and gates were unaffected by that correction.

The complete benchmark contained 184 matched contrasts. Eight source-labelled favourable families
were used for the family gate: MPTR, GSE176206 SOKM, mouse 2c/7c, human 2c/7c, human OSK/O4YRSK,
early human OSKM, RPE OSK and RPE GSTA4. Source labels summarized authors' interpretations and were
not ground truth. Broad transportability required at least four of six attempted methods to orient
both age controls and median pairwise sign agreement of at least 0.80. A universal direction
additionally required at least four methods to have positive family means in at least 80% of the
eight families. These are locked descriptive gates, not significance tests.

For gene-level concordance, family effects were mean treatment-minus-control log-expression vectors,
standardized within family on 9,543 common mapped genes. We computed all pairwise Spearman
correlations and leave-one-family-out correlations between each held family and the mean of the other
seven. Genes were not treated as biological replicates.

## Post-review Extension C1

For each method pair, sign agreement was first computed within each of the 14 benchmark families and
then averaged across families; the reported statistic is the median across ten method pairs. The
cluster bootstrap resampled 14 families with replacement 10,000 times (seed 3092026). The
leave-one-family-out analysis recomputed the equal-family statistic after omitting each family. Gate
sensitivity evaluated agreement thresholds from 0.50 to 1.00 in increments of 0.025 and favourable-
family requirements from four to eight of eight.

Strata were fixed by species, intervention class and cellular context. Split-half reliability used
balanced partitions of contrast units for families with at least four units. Where both treatment
and control contained at least two biological units, a separate arm-split analysis reconstructed two
independent effect vectors. Spearman–Brown values were recorded, but the manuscript reports raw
split-half correlations to avoid implying full-study reliability.

Pasta was reconstructed from the authors' MIT-licensed model package at commit
`58bcc7a69ee97f2dc9e3623ac86538c251ddd498`, using its 8,113-gene list, median imputation,
within-sample ranks and frozen ridge coefficients [14]. Age-shift predictions were sign-reversed so
that positive values denote a youth-associated direction. Package files and their SHA-256 hashes,
implementation validation and the licensing decision for tAge are recorded in the release audit.

## Software and reproducibility

Analyses used Python 3.12 with NumPy, pandas, SciPy, scikit-learn, h5py, anndata, rdata, matplotlib
and seaborn in a frozen local environment. All stochastic procedures used recorded seeds. Tests cover
metadata parsing, fold construction, rank scoring, orthology, normalization, contrasts, method
orientation and the CV repair. The release bundle contains code, configurations, hashes, aggregate
results, complete contrast tables and figure sources, but excludes raw public matrices and
participant-linked records.

# Data availability

All expression data are public under the GEO accessions listed above. The submission bundle includes
exact source URLs, retrieval dates, file sizes and SHA-256 hashes. Raw matrices are not redistributed.
The reviewed manuscript repository will be made publicly available at
<https://github.com/otcan/transcriptomic-reprogramming-portability>. The immutable archival
DOI will be inserted before submission: **[ZENODO DOI — REQUIRED BEFORE SUBMISSION]**.

# Code availability

The complete analysis and figure code, frozen configurations, tests and aggregate result tables will
be made publicly available at
<https://github.com/otcan/transcriptomic-reprogramming-portability>. Immutable commits and
manifests are included. No language model or post-cutoff biological database was used as a model
feature or to rank targets.

# Acknowledgements

This research received no external funding. Computational resources were provided by Metastate Bio
Inc.

OpenAI Codex was used under human direction to assist with code drafting, analysis orchestration,
literature organization and manuscript drafting. The authors inspected the source data, code,
computations, citations, figures and claims and retain full responsibility for the work and its
conclusions.

# Author contributions

O.Ü. conceived the study, defined the research questions and claim boundaries, directed the
computational programme, provided resources, interpreted the results through iterative review,
revised the manuscript and accepts responsibility for the integrity of the work. CRediT roles:
Conceptualization, Methodology, Resources, Supervision, Project administration, and Writing – review
& editing.

# Competing interests

O.Ü. is the founder of Metastate and is affiliated with Metastate Bio Inc, which develops commercial
computational biology, biomarker and modelling products and services that could benefit from the
publication of this work. The author declares no other competing interests.

# Ethics statement

This study reanalysed public, de-identified molecular datasets and involved no new human or animal
recruitment or experimentation by the author. The original studies' approvals govern source-data
collection.

# References

1. López-Otín, C. *et al.* Hallmarks of aging: An expanding universe. *Cell* **186**, 243–278 (2023). https://doi.org/10.1016/j.cell.2022.11.001
2. Rutledge, J., Oh, H. & Wyss-Coray, T. Measuring biological age using omics data. *Nature Reviews Genetics* **23**, 715–727 (2022). https://doi.org/10.1038/s41576-022-00511-7
3. Teschendorff, A. E. & Horvath, S. Epigenetic ageing clocks: statistical methods and emerging computational challenges. *Nature Reviews Genetics* **26**, 350–368 (2025). https://doi.org/10.1038/s41576-024-00807-w
4. Fleischer, J. G. *et al.* Predicting age from the transcriptome of human dermal fibroblasts. *Genome Biology* **19**, 221 (2018). https://doi.org/10.1186/s13059-018-1599-6
5. Tsitsipatis, D. *et al.* Transcriptomes of human primary skin fibroblasts of healthy individuals reveal age-associated mRNAs and long noncoding RNAs. *Aging Cell* **22**, e13915 (2023). https://doi.org/10.1111/acel.13915
6. Gill, D. *et al.* Multi-omic rejuvenation of human cells by maturation phase transient reprogramming. *eLife* **11**, e71624 (2022). https://doi.org/10.7554/eLife.71624
7. Roux, A. E. *et al.* Diverse partial reprogramming strategies restore youthful gene expression and transiently suppress cell identity. *Cell Systems* **13**, 574–587.e11 (2022). https://doi.org/10.1016/j.cels.2022.05.002
8. Mitchell, W. *et al.* Multi-omics characterization of partial chemical reprogramming reveals evidence of cell rejuvenation. *eLife* **13**, e90579 (2024). https://doi.org/10.7554/eLife.90579
9. Schoenfeldt, L. *et al.* Chemical reprogramming ameliorates cellular hallmarks of aging and extends lifespan. *EMBO Molecular Medicine* **17**, 9 (2025). https://doi.org/10.1038/s44321-025-00265-9
10. Mitchell, W. *et al.* In vivo chemical reprogramming is associated with a toxic accumulation of lipid droplets hindering rejuvenation. *Aging Cell* **25**, e70390 (2026). https://doi.org/10.1111/acel.70390
11. Lu, A. T. *et al.* Prevalent mesenchymal drift in aging and disease is reversed by partial reprogramming. *Cell* **188**, 5895–5911.e17 (2025). https://doi.org/10.1016/j.cell.2025.07.031
12. Lu, Y. R. *et al.* Reprogramming factors activate a non-canonical oxidative resilience pathway that can rejuvenate RPEs and restore vision. *bioRxiv* (2025). https://doi.org/10.1101/2025.08.30.673239
13. Tyshkovskiy, A. *et al.* Universal transcriptomic hallmarks of mammalian ageing and mortality. *Nature* **654**, 173–188 (2026). https://doi.org/10.1038/s41586-026-10542-3
14. Salignon, J. *et al.* Pasta, a versatile transcriptomic clock, maps the chemical and genetic determinants of aging and rejuvenation. *Advanced Science*, e76740 (2026). https://doi.org/10.1002/advs.76740
15. Mitchell, W. *et al.* The mitochondria-targeted peptide therapeutic elamipretide improves cardiac and skeletal muscle function during aging without detectable changes in tissue epigenetic or transcriptomic age. *Aging Cell* **24**, e70026 (2025). https://doi.org/10.1111/acel.70026
16. Gene Ontology Consortium. The Gene Ontology knowledgebase in 2023. *Genetics* **224**, iyad031 (2023). https://doi.org/10.1093/genetics/iyad031
17. Min, B. *et al.* Rejuvenation potential of developmental genes downregulated in aging. *International Journal of Stem Cells* (2026). https://doi.org/10.15283/ijsc25144
18. Benjamini, Y. & Hochberg, Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society B* **57**, 289–300 (1995). https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
19. Pedregosa, F. *et al.* Scikit-learn: machine learning in Python. *Journal of Machine Learning Research* **12**, 2825–2830 (2011).
20. de Lima Camillo, L. P. *et al.* pyaging: a Python-based compendium of GPU-optimized aging clocks. *Bioinformatics* **40**, btae200 (2024). https://doi.org/10.1093/bioinformatics/btae200
