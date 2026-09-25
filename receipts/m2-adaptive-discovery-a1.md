# M2 adaptive discovery checkpoint — before release tag

Date: 2026-09-25  
Validation matrices opened: **no**  
Named historical target ranks inspected: **no**

## Fixed model provenance

- Protocol tag: `protocol-v1.0`, commit `8be128cc5ab5f2cfc3147a0ee569fa59ff256b79`.
- Protocol-v1 Y construction: failed with zero selected genes; preserved unchanged.
- Adaptive model tag: `adaptive-y-a1`, commit `b7109948380cfd14acba5cad26f98b5acce7a897`.
- A1 model: 10,532 protein-coding within-sample rank features; frozen alpha
  `3.1622776601683795`; higher `Y_A1` is more youth-associated.

## Discovery outcomes

- GSE165177: 13 successful matched MPTR contrasts across three donors. Donor-mean Delta Y_A1 is
  0.43696 with percentile-bootstrap interval 0.16285 to 0.61004. Identity Delta is -0.00356
  (-0.00815 to 0.00017), pluripotency Delta is 0.00256 (-0.00130 to 0.00507), and aggregate damage
  Delta is -0.00310 (-0.00392 to -0.00263). The prespecified separation and guarded rules pass.
- GSE165177: 21 failed-to-reprogram matched contrasts are exported as a negative comparator; their
  mean Delta Y_A1 is 0.06866 and is not treated as an independent 21-donor estimate.
- GSE246954: young-minus-old untreated controls orient correctly (Delta Y_A1 0.24545; interval
  0.20114 to 0.29638). Contrary to a general youthward interpretation, old 2c and 7c effects are
  -0.39583 and -0.46266, respectively, with all four paired replicates negative for both cocktails.
- GSE176206: 64 adipogenic and 64 MSC pseudobulks yield 120 treated-versus-NT factor contrasts.
  The source-specific identity programmes contain 61 and 75 genes after enforcing detection in at
  least 25% of cells in every NT pool. Delta Y_A1 and Delta I have Spearman rho 0.01344. Factor
  contrasts are descriptive across the limited independent pools.

## Interpretation boundary

Discovery supports a narrow claim that the frozen axes separate youth-associated movement from
identity/pluripotency movement in the human MPTR dataset. It does **not** support a general detector
of rejuvenating interventions: the chemical and most factor-screen results expose substantial
transportability limits. No output establishes biological age reversal, causality, safety,
functional rejuvenation, or therapeutic efficacy.

The next irreversible action is a complete discovery-release hash/commit/tag. Only after that may
the predeclared target ranks be inspected or validation expression matrices be acquired.

## Release manifest

The 63-file release manifest is `receipts/m2-discovery-a1-release-files.sha256`; its SHA-256 is
`73c7dc5c37d1b7f3170371d0bb7bb629ebbf0b70f07a8b1ac4a0aa266d29d0a7`. The intended immutable tag is
`discovery-a1-v1`.
