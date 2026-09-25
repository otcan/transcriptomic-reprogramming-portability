# Dataset registry and temporal split

This registry governs dataset role, not scientific quality. A study can be useful yet still be
ineligible for discovery because it became public after the cutoff.

## Eligibility clock

The discovery cutoff is **2024-12-31 at 23:59:59 UTC**. Eligibility is determined by the first
public availability of the exact data, metadata, and external biological annotations used. A current
web page, updated record, or later database release cannot silently contribute information to the
discovery model. Every downloaded file will receive its source URL, retrieval time, checksum, file
size, and licence/terms entry before reading its content.

## Fixed roles

| Role | Accession / record | Public date verified | Species / modality | Use | Boundary |
| --- | --- | ---: | --- | --- | --- |
| Discovery anchor | GSE165180 | 2022-04-19 | Human; RNA-seq and methylation SuperSeries | Maturation-phase transient reprogramming trajectory | Enumerate its eligible subseries and biological-unit metadata in the input manifest before use. |
| Discovery anchor | GSE176206 | 2022-06-10 | Mouse; scRNA-seq | Factor subsets, identity/rejuvenation separation | Inference uses animal-level pseudobulk or hierarchy, never cell counts as replication. |
| Discovery anchor | GSE247199 | 2024-03-21 | Mouse; transcriptome, methylation, proteome, phosphoproteome, metabolome | Orthogonal chemical-reprogramming and cross-modal check | No 2025+ annotations or conclusions may be imported. |
| Discovery reference | GSE113957 | 2018-11-16 | Human; bulk RNA-seq | Age-associated reference direction | Include the 133 healthy donors aged 1-94; exclude the 10 HGPS samples from the primary reference. |
| Discovery reference | GSE226189 | 2023-07-06 | Human; bulk RNA-seq | Independent age-associated reference direction | Include the 82 healthy donors aged 22-89. |
| Temporal validation | GSE297233 / GSE297234 | 2025-08-14 | Human; RNA-seq / scRNA-seq | Trajectory validation | Quarantined from construction and tuning. |
| Temporal validation | GSE297984 | 2025-06-10 | Human dermal fibroblasts; bulk RNA-seq | Two-donor-line 2c/7c positive challenge at days 6 and 14 | Quarantined; culture repetitions are nested within 56- and 83-year donor lines. |
| Temporal validation | GSE300625 | 2026-06-22 | Mouse liver/kidney; bulk RNA-seq | Prespecified adverse in-vivo challenge | Quarantined; must not be labelled guarded rejuvenation if the framework is specific. |
| Temporal validation | GSE304043 / GSE304044 | 2026-06-11 | Mouse/human RPE; RNA-seq and functional-genomics SuperSeries | State-direction and mechanism validation | Quarantined from construction and tuning. |
| Secondary temporal validation | GSE276656 | 2026-02-10 | Mouse mPFC; single-nucleus RNA/ATAC | Cell-type generalization sensitivity | Not a primary gate; submitted in 2024 but exact data were not public until 2026. |
| Literature target benchmark | Min et al., 2026, DOI 10.15283/ijsc25144 | 2026 | Experimental fibroblast study; public expression matrix not located | Exact ranks for DDX21, TOMM70A, SERBP1, PHGDH | Secondary target-name benchmark only; no transcriptomic validation claim without a primary data deposit. |

## Mandatory registry completion before lock

1. Record the exact GSE165180 subseries/files used; do not substitute a related study unnoticed.
2. Record the frozen releases, download dates, licences, and checksums for gene annotation,
   orthology, ontology, identity, pluripotency, and damage programmes. No interaction network is used
   in the primary model or rank.
3. Record biological-unit mapping and all excluded samples before looking at score results.
4. Create separate `discovery/` and `validation/` data locations with permissions and manifests.

## Source verification notes

The accessions and public dates above were checked directly from GEO brief records on 2026-09-13.
GSE113957 and GSE226189 were selected before data retrieval by a documented eligibility rule: human
primary/dermal fibroblast transcriptomics, publicly released by the cutoff, known donor age, and
multiple healthy donors spanning adulthood. They are independent source studies; GSE113957's HGPS
samples are excluded from the primary age reference because the primary construct is normal ageing.
GSE176206's record describes 18 source samples and large processed H5AD files; the public record
must be used to determine the true experimental replication rather than the number of cells. The
GSE304043 series record contains seven source samples. These counts are metadata facts, not a
substitute for a biological-replicate audit. GEO metadata list 337 samples in GSE165180, 60 in
GSE247199, and 120 in GSE304044; SuperSeries totals do not override the biological-unit audit.

GSE297984 contains 24 human fibroblast RNA-seq samples but only two aged donor lines (56 and 83
years); culture repetitions and time points are not independent donors. GSE300625 contains 20 mouse
tissue RNA-seq samples from a 7c-versus-vehicle in-vivo experiment. They were added before protocol
lock because together they provide positive and adverse temporal challenges. Their published
conclusions are already known, so these are outcome-informed historical benchmarks, never
prospective tests.
