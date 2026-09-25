# M1 discovery acquisition checkpoint

Date: 2026-09-25

Nine predeclared discovery files were downloaded as opaque inputs: 6.6 GB total. All independent
SHA-256 checks passed. The acquisition manifest is `inputs/discovery-manifest.tsv`, SHA-256
`c7a5ddeadd7b15806313b82265b22e9508cd0fb83fe991f15a9058953d9a589a`.

The manifest includes source URL, public date, retrieval time, byte size, SHA-256, biological unit,
condition mapping, inclusion decision, access boundary, and `content_opened=no`. Discovery inputs
remain outside Git. No validation-data directory exists.

The two selected GSE176206 screen files account for most storage (3,387,048,968 and 3,370,578,849
bytes). Single-pool SOKM H5ADs were not acquired because they cannot support biological-unit
inference and are not required by the locked primary factor-screen endpoint.

Next gate: commit this manifest and the deterministic acquisition code; then parsers may inspect
discovery file structure and produce QC summaries. Validation remains sealed.

