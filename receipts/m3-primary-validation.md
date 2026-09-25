# M3 frozen temporal validation — primary claim failed

Date: 2026-09-25  
Discovery release: `discovery-a1-v1` at `cdece0dec297800e7e3ab1aa1688f7f6d8a9d6b9`  
Validation refit: **none**

## Co-primary and trajectory results

- **GSE297984 positive human chemical challenge: failed.** Neither 2c nor 7c produced positive
  donor-mean and day-mean Delta Y_A1 across both aged donor lines and both days. Six of eight
  donor/day/cocktail effects were negative; only the 83-year line's 2c contrasts were positive.
- **GSE297234 trajectory order: failed.** GEO identifies GM00731 as the 96-year donor and GM23815 as
  the 22-year donor. Their baseline Y_A1 difference is correctly oriented (+0.85559 young-minus-old).
  In the old donor, day 3 is the first positive Delta Y_A1 (+0.30763), but it is also the first day
  exceeding both fixed risk margins (Delta I -0.08229; Delta P +0.09387).
- **GSE300625 adverse specificity: formal pass, non-informative.** Liver Delta Y_A1 is -0.11742
  (-0.19121 to -0.04130); kidney is -0.07520 (-0.17313 to 0.02494). Neither tissue can be called
  guarded because a tissue-identity axis was not defined. This is non-evaluability, not demonstrated
  specificity.

## Secondary results

- GSE297233: OSK Delta Y_A1 -0.08780 and O4YRSK -0.27928; both also lose fibroblast identity and
  gain the endogenous pluripotency programme.
- GSE304042: ARPE OSK-versus-GFP Delta Y_A1 -0.00054.
- GSE304043: GSTA4-versus-GFP Delta Y_A1 0.03244 (-0.63759 to 0.55661); aggregate D 0.00419
  (-0.00089 to 0.00859). No transcriptomic validation is supported by this endpoint.

## Decision

The preregistered positive-framework claim is closed as failed. The project may continue only as an
explicitly adaptive multi-method benchmark asking whether apparent rejuvenation directions are
portable across interventions, species, and datasets. It may not be presented as a successful
Metastate validation or target-discovery paper.

Evidence: `inputs/validation-manifest.tsv`, `results/validation/`, and
`code/score_validation.py`. All 11 repository tests passed before this receipt.
