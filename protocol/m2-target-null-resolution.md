# M2 implementation resolution — known-target matched null

Date: 2026-09-25  
Status: defined after the immutable discovery rank, before inspecting any named-target rank

The frozen specification requires 100,000 expression-decile-matched five-gene sets but did not
fully define the abundance statistic or repeated-decile sampling. This resolution closes that
implementation ambiguity without changing the candidate rank or target list.

1. Reload the two discovery age-reference matrices with the frozen sample and gene filters.
2. Convert each sample to within-sample percentile ranks across all retained genes.
3. For each candidate gene, calculate the median percentile rank across donors separately in each
   cohort, then average the two cohort medians. This is the matching abundance statistic.
4. Divide the immutable candidate universe into ten equal-count bins by abundance rank. Stable gene
   symbol order breaks an exact abundance tie only for bin assignment.
5. Remove all five known historical targets from the null pool. For each of 100,000 sets, sample one
   gene from the same abundance decile as each evaluable target. Sampling is without replacement
   within a set, including when multiple targets share a decile.
6. The statistic is the mean immutable candidate-rank percentile. The one-sided empirical p-value
   is `(1 + number(null >= observed)) / 100001`. Seed is the frozen `271828`.
7. Report every target's exact rank and percentile. A missing target is reported as not evaluable
   and is not replaced.

This is an outcome-informed historical recovery analysis, not a blinded or prospective target
discovery claim.
