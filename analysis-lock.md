# Analysis-lock procedure

No expression matrix is opened, normalised, plotted, or fitted until sections 1-3 are complete and
the resulting local protocol tag is created. After that tag, discovery files may be downloaded as
opaque bytes and hashed; they still may not be opened until their complete acquisition manifest is
committed. Validation matrices remain unavailable until the discovery release is immutable. The
user's 2026-09-25 instruction authorizes this internal lock; it does not authorize external release
or submission.

## 1. Freeze inputs

- [x] Complete `dataset-registry.md` with the exact normal-ageing references.
- [x] Save a machine-readable pre-acquisition `inputs/source-plan.tsv` with exact discovery URLs and
  roles. After protocol tag, create `inputs/discovery-manifest.tsv` with retrieval time, SHA-256,
  byte size, access terms, and biological-unit mapping before opening a file.
- [x] Reserve all post-2024 files for `inputs/validation-data/`, with a distinct manifest. Do not
  create that directory until the discovery release is immutable.
- [x] Freeze versions/releases of annotation, orthology, ontology, and gene sets. No primary network
  is used.

## 2. Freeze model specification

- [x] Define exactly how each of `Y`, `I`, `P`, and `D` is computed, including discovery-derived gene rules, transforms, missing-gene rule, and score direction.
- [x] Complete `model-specification.md`; changes after the protocol tag are recorded deviations.
- [x] Define margins, primary endpoints, target rank, null matching, seeds, and multiplicity families.
- [x] Write environment lockfiles; deterministic pipeline entry points are built before the first matrix is opened.
- [x] Specify biological-unit uncertainty and the conditions under which cross-study synthesis is allowed.

## 3. Freeze evaluation

- [x] Define GSE297984 and GSE297233/GSE297234 endpoints before matrix access.
- [x] Define GSE300625 and GSE304043/GSE304044 endpoints before matrix access.
- [x] Record the already-known target-evaluation list (GSTA4, DDX21, TOMM70A, SERBP1, PHGDH) and
  commit its hash. It is not presented as concealed from the team.
- [ ] Commit and tag the full specification, then archive its hash before any discovery matrix is opened.

## 4. Honest wording rule

The algorithm can be temporally out-of-sample, but the research team already knows the later studies,
reported directions, GSTA4, and the 2026 RIFs. The paper must not use "blinded", "prospective", or
"predicted before publication" for this exercise. It may say "evaluated using a frozen pre-2025
data-only model against later public datasets".
