# Contemporary-clock audit for Extension C1

Date: 30 September 2026

## Pasta

- Primary article: Salignon et al., *Advanced Science* (2026), DOI
  `10.1002/advs.76740`.
- Official code: `https://github.com/jsalignon/pasta`.
- Pinned commit: `58bcc7a69ee97f2dc9e3623ac86538c251ddd498`.
- Repository licence: MIT.
- Model inputs: 8,113 gene symbols, median imputation from the packaged training reference,
  within-sample ranks and a frozen ridge model.
- Independent implementation check: the reconstructed Python scorer reproduced all three packaged
  README example predictions; maximum absolute error was `4.73e-11`.
- Training overlap: searches of the public article and packaged repository found no exact benchmark
  GEO accession. The article describes 21 healthy-donor training datasets, but the complete
  sample-level accession list was not available in the audited materials. Disjointness is therefore
  not claimed.
- Interpretation: post-review and outcome-aware; not part of the historical B1 gate.

## tAge

- Primary article: Tyshkovskiy et al., *Nature* (2026), DOI
  `10.1038/s41586-026-10542-3`.
- Official code: `https://github.com/Gladyshev-Lab/tAge`.
- The official repository states that the MGB Open Access License 1.0 permits non-commercial
  academic use and requires a separate agreement for commercial use.
- Because the author is commercially affiliated and no commercial agreement was obtained, tAge was
  excluded from C1. This is a rights boundary, not a scientific selection based on its outputs.

