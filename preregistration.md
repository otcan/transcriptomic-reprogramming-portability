# Frozen-analysis protocol — release candidate 1.0

Freeze date: 2026-09-25; hash and local tag recorded in the M0 lock receipt
Owner: Metastate Research
Execution authorization: user instruction of 2026-09-25

## Question

Can an interpretable model built only from public resources available by 2024-12-31 distinguish a
human-fibroblast youth-associated transcriptional direction from somatic-identity loss, endogenous
pluripotency, and measured stress during partial reprogramming, then behave coherently in later
positive and adverse public datasets without refitting?

## Eligible discovery data

Discovery is limited to GSE113957, GSE226189, GSE165176/GSE165177, GSE176206, and
GSE246954/GSE247199. Exact files, source dates, checksums, biological units, and exclusions are
frozen in the discovery manifest. No dataset may be substituted after lock without a dated
deviation.

## Sealed temporal data

Expression matrices for GSE297984, GSE297233/GSE297234, GSE300625, GSE304043/GSE304044, and
GSE276656 remain unopened until the discovery score definitions and complete target rank are hashed
and committed. Metadata and primary-paper methods may be inspected now to define endpoints. The
published conclusions and target names are known to the team; validation is historical,
outcome-informed, and temporally external—not blinded or prospective.

## Unit of analysis

The independent donor, animal, culture, or independently prepared experimental pool is the unit.
Cells and technical replicates are nested observations. If source metadata cannot establish a
biological unit, that comparison becomes descriptive and cannot satisfy an inferential gate.

## Model and ranking

`model-specification.md` release 1.0 defines all transformations, selection rules, scores, margins,
contrasts, candidate rank, seeds, and validation endpoints and is incorporated by reference. No LLM,
post-2024 database, current network, or manual target knowledge contributes a model feature.

## Primary hypotheses and gates

1. **Age-axis construct validity:** five-fold held-out Spearman `rho(Y, age) <= -0.35` pooled, with
   negative study-specific correlations no greater than `-0.15`.
2. **State separation:** at least one replicated discovery intervention has a 95% interval for
   `Delta Y > 0` while its Delta I and Delta P intervals remain inside the fixed preservation
   margins.
3. **Positive temporal validation:** both GSE297984 donor lines pass the frozen directional Delta Y
   endpoint without refitting; no culture replicate is counted as an additional donor.
4. **Adverse specificity:** neither GSE300625 tissue is called guarded youthward.

The positive-framework claim requires all four gates. Failure leads to a negative/benchmark paper
or project stop; thresholds and axes will not be redefined to recover the claim.

## Secondary hypotheses

- GSE297233/234 follows a youthward-before-dedifferentiation trajectory.
- GSTA4 moves Y/D in the expected direction and/or ranks above matched controls.
- The known five-gene later-target set is enriched in the immutable candidate rank.
- GSE247199 transcript/protein/metabolite directions are coherent at the pathway level.

Secondary outcomes cannot override a failed primary gate.

## Controls

- Young versus old fibroblasts anchor Y direction.
- Negative-control and failed-to-reprogram arms anchor nonspecific perturbation.
- Full iPSC states anchor I loss and endogenous P rise.
- Label permutations occur only within study blocks and preserve replicate structure.
- GSE300625 is a real adverse-treatment challenge, not a simulated negative.

## Statistics

Primary effects use biological-unit mean differences and 10,000-resample bias-corrected bootstrap
intervals when n>=3 per arm. Paired designs resample pairs. Smaller experiments are descriptive.
Spearman uncertainty is donor-level. BH correction is applied within the families listed in the
model specification. Cross-study summaries show study-specific effects and a random-effects model
only when at least three independent study families estimate the same construct.

## Missingness and exclusions

- Score only programmes with at least 15 observed genes and 50% coverage.
- Exclude ambiguous orthologues and unresolved biological units from primary inference.
- Apply source-independent expression/QC rules without reference to score direction.
- Record every exclusion before computing state scores.

## Stopping and wording rules

No output may claim causal rejuvenation, direct biological age, safety, clinical efficacy,
therapeutic efficacy, peptide efficacy, or target causality. A high rank is an expression-derived
priority only. If licence or participant-data constraints prevent a reproducible lawful release,
release code/manifests/aggregates only or stop.

## Deviations

After the lock hash, every change requires a dated entry in `results-and-deviations.md` with reason,
scope, whether it was made before or after seeing outcomes, and whether the result is primary,
sensitivity, or exploratory. Historical files are preserved.
