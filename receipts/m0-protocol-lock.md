# M0 protocol lock — Paper 2

Date: 2026-09-25
State: ready for local commit and tag `protocol-v1.0`

The state axes, biological units, discovery cutoff, candidate rank, positive and adverse validation
endpoints, missingness rules, margins, multiplicity families, seeds, orthology, ontology, and claim
limits are fixed before any expression matrix is opened. Metadata and primary-paper methods were
audited; published validation conclusions and target names were already known and are disclosed.

The protocol manifest is `receipts/m0-protocol-lock-manifest.sha256`. Its SHA-256 is recorded in the
tag annotation. The known-target file hash is
`7065869ea326158719afd8fce5931f33022798693af14ff78a1aea52ae103fd7`.

Pre-lock checks:

- JSON configuration parses.
- Metadata parser smoke test resolves all 18 GSE176206 libraries.
- One automated parser test passes in the frozen Python environment.
- Annotation files total 149 MB and match `receipts/annotation-sha256.txt`.
- No expression matrix has been downloaded or opened.
- No validation-data directory exists.

After the tag, discovery files may be downloaded and hashed as opaque inputs. They may be opened
only after the acquisition manifest and deterministic pipeline entry points are committed.

