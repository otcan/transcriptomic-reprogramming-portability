# State and target-ranking specification — release candidate 1.0

This is an interpretable, rank-based model. It is not a foundation-model embedding, causal graph,
clinical clock, or pooled batch classifier. The algorithm is frozen before any expression matrix is
opened. Gene membership produced by a discovery-only selection rule is an output of that frozen
algorithm, not manual curation.

## 1. Common gene representation

- Human gene identifiers are reduced to approved gene symbols using the source annotation and an
  archived identifier table recorded in the manifest. Duplicate symbols are collapsed by their
  median expression after source-level normalisation.
- Primary mouse-to-human mapping uses NCBI HomoloGene build 68 (file dated 2014-05-06). Retain a
  group only when it contains exactly one human (`tax_id=9606`) and one mouse (`tax_id=10090`) gene.
  All one-to-many and ambiguous mappings are excluded from primary analyses.
- Ensembl Compara release 111 one-to-one mappings are a labelled sensitivity only. A current or
  post-2024 mapping may not replace the primary mapping.
- Protein-coding status follows source annotations. If a source does not carry biotype, the frozen
  annotation table determines eligibility.

## 2. Within-sample rank score

For each sample or pseudobulk biological unit, rank all eligible detected genes by expression using
average ranks for ties and convert to percentile ranks in `[0,1]`. For a one-sided programme `G`,

`R(G) = mean(percentile_rank(g), g in G) - 0.5`.

For a signed programme with positive set `G+` and negative set `G-`,

`R(G+,G-) = 0.5 * [mean(rank(G+)) - mean(rank(G-))]`.

At least 15 genes and 50% of the frozen programme must be observed; otherwise the score is missing.
No missing gene is imputed. Higher scores always mean more of the named construct.

## 3. State vector

For biological unit `b` report `S(b) = (Y, I, P, D)`.

### Y — youthward age-associated state

1. Use only healthy donors in GSE113957 and all eligible donors in GSE226189.
2. Convert each reference to `log2(CPM + 0.5)` when counts are available; use deposited FPKM only
   for GSE113957 because no common raw-count matrix is supplied. Filter genes expressed above the
   source minimum (CPM >= 1 or FPKM >= 1) in at least 20% of biological units.
3. Within each study and gene, fit ordinary least squares on standardised expression:
   `z(expression) ~ z(age) + sex`. Unknown sex is an explicit indicator; no age-derived sequencing
   group is adjusted away.
4. A final Y gene must have the same age-effect sign and Benjamini-Hochberg `q <= 0.10` in both
   studies. Young-up genes have negative age coefficients; old-up genes have positive coefficients.
   At least 25 genes per direction are required or Y construction fails.
5. `Y = R(young-up, old-up)`. Higher Y is greater concordance with the normal human fibroblast
   youth-associated direction. Genes receive equal weight.
6. Construct validity is estimated with one fixed five-fold cross-validation (`seed=1729`),
   stratified by study, sex, and age quintile. Feature selection is rerun inside each training fold.
   The prespecified pass is pooled held-out Spearman `rho(Y, age) <= -0.35`, with both study-specific
   correlations negative and no study-specific correlation greater than `-0.15`.

Y is not a biological-age estimate. No result is reported in “years younger.”

### I — retained somatic identity

Two source-specific programmes are allowed; they are never pooled into a universal identity score.

1. **Human/mouse fibroblast identity.** In GSE165176, separately for donor lines N2, N3, Y1, Y2,
   O1, and O2, compare each day-0 dermal-fibroblast baseline with the median of that line's late
   SSEA4-positive day-47/day-54 samples. Select genes whose median decrease is at least 1 deposited
   log2-expression unit, whose decrease is at least 1 unit in five of six donor lines, that are
   expressed >=1 in at least five of six baselines, and that are not selected for P or a D
   programme. This paired consistency rule avoids treating repeated time points as independent
   donors. Map to mouse only through the primary orthology table.
2. **GSE176206 cell-type identity.** For each adipogenic and MSC source, compare untreated/control
   pseudobulk expression with the other cell type using only control pools. Select the 100 largest
   positive log-fold-change genes subject to log2 fold change >=1 and expression in at least 25% of
   cells in each available pool. This programme is used only inside GSE176206. If fewer than two
   independent experimental pools exist, results are descriptive and carry no inferential p-value.

`I = R(identity programme)`. Higher I means retained expression of the specified baseline programme,
not preservation of all cell functions.

### P — endogenous pluripotency/dedifferentiation programme

In GSE165176, use the same six donor-line paired late-versus-baseline contrast as I. Select genes
whose median increase is at least 1 deposited log2-expression unit, whose increase is at least 1
unit in five of six lines, and that are expressed >=1 in at least five of six late-state medians.
Exclude `POU5F1`, `SOX2`, `KLF4`, `MYC`, vector/transgene features, and genes without a standard
endogenous gene identifier. `P = R(P programme)`.

The primary P score therefore cannot rise merely because introduced OSKM transcripts are present.
An unfiltered score including endogenous OSKM symbols is reported only as a sensitivity. P is not a
tumorigenicity or clinical-safety score.

### D — measured stress/damage programmes

Use Gene Ontology release 2024-01-17 (`go-basic.obo`) with the corresponding human GOA and mouse MGI
GAF files. Include annotations to the term or any `is_a`/`part_of` descendant, exclude `NOT` and
evidence code `ND`, and freeze these four components:

- GO:0006979 — response to oxidative stress;
- GO:0006974 — cellular response to DNA damage stimulus;
- GO:0034976 — response to endoplasmic reticulum stress;
- GO:0090398 — cellular senescence.

Each component is `R(G)`. `D` is their unweighted mean only when all four are available; components
remain visible in every primary table. Higher D means greater concordance with these programmes,
not total cellular damage.

## 4. Intervention contrasts and guarded state

For a treatment and its matched contemporaneous control, compute `Delta S = treatment - control`.
Primary effects are differences in biological-unit means with percentile-bootstrap 95%
intervals when each arm has at least three biological units (`seed=1729`, 10,000 resamples). With
fewer than three units, report the effect descriptively without a confidence interval or p-value.

A contrast is called **guarded youthward** only when all conditions hold:

- the 95% interval for `Delta Y` lies above 0;
- the 95% interval for `Delta I` lies above `-0.05` rank-score units;
- the 95% interval for `Delta P` lies below `+0.05` rank-score units;
- the 95% interval for `Delta D` lies below 0.

If an axis is unavailable, guarded status is “not evaluable,” never passed. The fixed 0.05 margins
are 5% of the complete rank-score range and are not estimated from validation data.

## 5. Discovery tests

- GSE165177: successful MPTR versus matched negative-control fibroblasts by donor, experiment, and
  duration; failed-to-reprogram and full iPSC states are negative/trajectory comparators.
- GSE176206: factor-subset and SOKM effects by deposited experimental pool. Single cells are nested
  observations only. Any factor result with fewer than two independent pools is descriptive.
- GSE246954: old 2c and old 7c versus old untreated; young untreated versus old untreated is the
  age-direction control. Replicate labels are biological units unless source documentation proves
  they are technical.

The separation claim requires at least one replicated discovery contrast with `Delta Y > 0` by its
95% interval while meeting both I and P preservation margins. D is reported but is not required for
this narrower separation test. A guarded-youthward claim requires all four conditions.

## 6. Frozen candidate ranking

The primary rank deliberately uses no interaction network and no language-model knowledge.
Candidates are one-to-one human/mouse protein-coding orthologues detected in both age references,
GSE165177, and GSE246954.

For gene `g`:

- `A_g`: percentile rank of the absolute fixed-effect meta-z for age, set to zero unless age-effect
  signs agree across the two human references.
- `R_g`: percentile rank of the median signed reversal `-sign(beta_age) * z(delta_expression)` over
  successful MPTR contrasts and old-mouse 2c/7c contrasts. A component is zero if fewer than two
  discovery families contribute or their median reversal is non-positive.
- `C_g`: proportion of eligible discovery contrasts that reverse the age direction, minus the
  proportion that reverse it in failed-to-reprogram/negative comparator contrasts, truncated to
  `[0,1]` and converted to a percentile rank.

The immutable priority is `T_g = 0.40 A_g + 0.40 R_g + 0.20 C_g`. Ties receive the worst shared
rank. The complete table, components, candidate universe, code/config hashes, and input hashes are
exported once before validation matrices are opened.

GSTA4, DDX21, TOMM70A, SERBP1, and PHGDH are evaluated as an explicitly known five-gene historical
benchmark. Report exact ranks and percentiles. The single prespecified enrichment test compares
their mean percentile to 100,000 expression-decile-matched five-gene sets (`seed=271828`). This
secondary test does not establish causality and cannot rescue a failed state framework.

## 7. Temporal validation endpoints

- **GSE297984 positive challenge:** collapse repeated cultures within donor line, day, and treatment;
  then compare 2c and 7c with the matching control separately for the 56- and 83-year donor lines.
  The study has only two donor lines, so it receives no population-level p-value or confidence
  interval. The family passes if both donor lines have Delta Y > 0 for at least one cocktail and the
  same cocktail's aggregate direction agrees at both days. Guarded status is descriptive only.
- **GSE297233:** OSK and O4YRSK induction versus matched no-doxycycline controls. Because n=2 per
  arm, report exact effects without inferential claims.
- **GSE297234 trajectory:** score old and young donor libraries at days 0,3,7,10. Because there is
  one donor per age and one library per time, report trajectories descriptively. A favourable order
  requires old-donor Y to increase before either P increases by >0.05 or I decreases by >0.05.
- **GSE300625 adverse challenge:** within kidney and liver, 7c versus vehicle. A specific framework
  must not call either tissue guarded youthward. This is the prespecified false-positive endpoint.
- **GSE304043/GSE304044:** GSTA4 effects are tested on Y and each D component; RPE identity is not
  inferred from the fibroblast I programme. GSTA4 rank is reported independently of expression
  effects.

The positive-validation and adverse-challenge families are co-primary and interpreted together.
GSE276656 and the 2026 RIF experiment are secondary only.

## 8. Multiplicity and sensitivity

- BH adjustment is performed separately for discovery separation, positive temporal validation,
  adverse validation, and target recovery. Every tested contrast remains in exported tables.
- Leave-one-discovery-study-out Y application, alternative normalisation, Ensembl-111 orthology,
  inclusion of IEA GO evidence, cell-cycle exclusion, and unfiltered P are sensitivity analyses.
- No sensitivity result replaces a failed primary endpoint.
