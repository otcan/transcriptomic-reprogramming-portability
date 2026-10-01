# Transcriptomic partial-reprogramming portability benchmark

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23078079.svg)](https://doi.org/10.5281/zenodo.23078079)

- Archived release v1.0.0: <https://doi.org/10.5281/zenodo.23078079>
- Preprint (not peer reviewed): <https://doi.org/10.5281/zenodo.23078081>

This repository is the reproducibility release for:

> **Five chronological-age transcriptomic directions fail a locked portability benchmark across
> partial-reprogramming interventions**

The study asks whether transcriptomic directions trained to predict chronological age provide a
portable scalar readout across public partial-reprogramming experiments. The original positive
four-axis framework failed its frozen temporal-validation and target-recovery gates. A subsequent,
explicitly adaptive benchmark compared six attempted age-direction models (five evaluable) across
184 intervention contrasts.

All five evaluable models correctly oriented two independent young-minus-old controls. Their median
pairwise agreement on intervention direction was 0.744565, below the locked 0.80 gate, and no model
was youthward in at least 80% of eight source-labelled favourable families. The supported conclusion
is therefore non-portability for these evaluated directions—not absence of functional, epigenetic or
other biological effects in the source studies. Post-review family weighting, reliability and
feature-coverage analyses are retained separately under Extension C1; exploratory Pasta oriented both
age controls but was positive in four of eight favourable-labelled families.

## Start here

- `paper/manuscript.pdf` — submission manuscript.
- `paper/supplement.pdf` — supplementary information.
- `manuscript/` — editable source, figures and publication tables.
- `REPRODUCIBILITY.md` — environment, inputs and execution order.
- `protocol/` — frozen protocol and adaptive amendments.
- `code/` — acquisition, analysis, benchmark, tests and claim verification.
- `results/` — redistribution-safe aggregate and gene-level outputs.
- `inputs/` — source URLs and cryptographic hashes; raw matrices are not redistributed.
- `receipts/` — temporal-lock and verification evidence.

## Quick verification

Use Python 3.12:

```bash
python -m venv .venv
.venv/bin/pip install -r environment/requirements-frozen.txt
PYTHONPATH=code .venv/bin/python -m pytest -q code/tests
.venv/bin/python code/verify_manuscript_claims.py
```

Expected result: 19 tests pass and eleven manuscript-claim groups pass.

## Claim boundary

The release does not establish biological-age reversal, causal rejuvenation, clinical validity,
therapeutic efficacy, safety, peptide efficacy or target causality. Later-study conclusions were
known during design, so the historical evaluation is temporally external but neither blinded nor
prospective. Adaptive Benchmark B1 was defined after the positive framework failed and is labelled
adaptive throughout.

## Data and privacy boundary

All source expression data are available from the GEO accessions and URLs in `inputs/`. This release
does not redistribute upstream raw matrices, source H5AD/H5 files, sample-level expression,
participant-linked records or the original source-unit labels in the 184-contrast table. The public
contrast table uses stable `C0001`–`C0184` identifiers.

## Licence

Source code is Apache-2.0. Manuscripts, protocols, figures and newly authored aggregate tables are
CC BY 4.0. Third-party datasets are not included or relicensed. See `LICENSE` and `LICENSES/`.

## Citation

See `CITATION.cff`. A versioned archival DOI will be added to the manuscript and citation metadata
after the reviewed release is deposited.
