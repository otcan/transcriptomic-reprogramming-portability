# Project charter — Separating Cellular Rejuvenation from Dedifferentiation Using a Multidimensional State-Transition Model of Partial Reprogramming

## Identity

- Project slug: `metastate-paper2-rejuvenation`
- Owner: Metastate Research (corresponding author to be assigned before submission)
- Start date: scaffolded 2026-09-13; execution authorized 2026-09-25
- Working title: *Separating Cellular Rejuvenation from Dedifferentiation Using a Multidimensional State-Transition Model of Partial Reprogramming*
- Intended article class: computational systems-biology / methods-and-resource article
- Project state: M0 protocol/source lock in progress; no expression matrix has been downloaded or analysed by this project.

## Decision question

Can public transcriptomic and multi-omic data support an auditable state-transition framework that
distinguishes youthward movement from identity loss, pluripotency, and stress, and that rejects an
adverse reprogramming intervention rather than scoring every large perturbation as rejuvenation?

## Primary scientific question

Using only data and biological resources publicly available on or before 2024-12-31, can a pre-specified state-transition model quantify movement toward an age-associated reference state independently from somatic-identity loss, pluripotency risk, and stress/damage; and does that frozen model generalize to later independent partial-reprogramming datasets?

## Claim boundary

**Allowed claim:**

The frozen, data-only framework provides temporally out-of-sample evidence that its multidimensional state scores and pre-committed ranking are associated with later observed partial-reprogramming trajectories and independently tested mechanisms, conditional on all prespecified gates being met.

**Explicitly blocked claims:**

- The framework has discovered a clinically effective rejuvenation therapy, peptide, or safe human intervention.
- A score is a direct measurement of biological age, causal rejuvenation, cell identity, or cancer risk.
- A target is validated because it ranks highly in a retrospective analysis.
- A later-dated dataset is a truly prospective experiment performed after this project began. It is a **temporal external validation**, not a prospective wet-lab validation.

## Falsification and stopping rule

- The claim is weakened if youthward movement is inseparable from somatic-identity loss, if effects reverse in leave-one-study-out analysis, or if validation effects are indistinguishable from negative/control signatures.
- Target recovery is secondary and is reported whether positive or negative; it cannot rescue a
  failed state-transition result.
- Stop the positive-framework manuscript if the age axis fails cross-validation, if no discovery
  contrast separates youthward movement from identity/pluripotency movement, or if the framework
  labels the prespecified adverse in-vivo challenge as guarded rejuvenation. A useful negative
  methods/benchmark report remains permitted.

## Planned evidence

- Data sources and access class: public discovery datasets released by 2024-12-31; public 2025-26 temporal validation datasets isolated from model construction. See `dataset-registry.md` and `data-access-ledger.csv`.
- Comparator/control: young versus old within matched study/cell type; vehicle/control versus
  partial-reprogramming intervention; failed/full-reprogramming states; and expression-decile-
  matched gene sets for the secondary target-rank benchmark.
- Primary outcome and metric: donor/animal/sample-level change in a four-axis score vector `S = (Y, I, P, D)` and its prespecified composite utility. Single cells are observations nested within biological replicates, not independent replicates.
- Validation unit and leakage controls: biological donor/animal/sample is the inferential unit. The post-2024 data, author conclusions, target names, and outcome labels may not enter feature selection, parameter tuning, candidate selection, or manual curation.
- Required review: computational-statistics, ageing-biology, reporting-consistency, and
  release/rights review before author-review readiness. Automated/AI review is not independent human
  peer review and will be labelled accurately.

## Next decision

Complete and hash protocol release 1.0, then acquire discovery matrices. The user's 2026-09-25
instruction to proceed authorizes internal execution through the pre-submission gate; it does not
authorize journal submission or public release.
