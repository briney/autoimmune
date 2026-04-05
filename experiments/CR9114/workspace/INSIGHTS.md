# CR9114 Insights

## Current best: Iteration 1 — Genotype RF (Spearman ρ = 0.9591)

A RandomForestRegressor trained directly on the 16-bit genotype vectors achieves
outstanding performance: ρ = 0.959, top-11 precision = 0.55, pairwise accuracy = 0.75.
CV Spearman = 0.949 (3×5-fold), so no overfitting.

### Feature importances
| Bit | HC Pos | Mutation | Importance |
|-----|--------|----------|------------|
| 10  | 75     | S→F      | 42.3%      |
| 9   | 74     | K→I      | 16.5%      |
| 4   | 52     | I→S      | 16.3%      |
| 6   | 58     | A→T      | 8.2%       |
| 1   | 29     | F→S      | 6.4%       |
| 7   | 59     | N→A      | 3.2%       |
| others | — | — | <2.2% each |

Bits 10, 9, 4 alone account for ~75% of predictive power. These are all FR3/CDR2
positions likely at or near the paratope.

### ddG tool calibration (10-variant subset)
- `baddg` and `stabddg` were run on a 10-sample subset.
- Both tools showed poor rank correlation with true affinity on this small sample;
  their predicted_score values don't clearly track h1_mean.
- Not yet incorporated into a model — unclear if they add signal beyond genotype.

## Open questions / next experiments

1. **Can ddG tools add signal on top of genotype RF?**  
   Run baddg/stabddg on all ~645 variants, then train a stacked model
   (genotype_RF prediction + ddG score → affinity). Given the genotype RF already
   captures ~96% rank correlation, any structural signal must be very precise to help.

2. **Pairwise accuracy (0.75) is weak relative to Spearman (0.959).**  
   The pairs are single-mutation differences. The RF may be uncertain in that
   fine-grained regime. A model that better captures single-mutation effects
   (e.g., additive linear model, or ddG deltas) might improve pairwise accuracy.

3. **Top-11 precision (0.55) leaves room to improve.**  
   Investigate which top-predicted variants are wrong. Are they consistently off in
   a correctable way (e.g., all missing a particular interaction)?

4. **Sequence embeddings (ESM2, AbLang2) as features?**  
   ESM2 or AbLang2 embeddings of the variant heavy chain, extracted at the mutated
   positions, could capture higher-order epistatic effects the RF misses.

## Strategy
The genotype RF is already excellent. Incremental improvements are likely. Focus on:
- Improving pairwise accuracy (most actionable gap)
- Stacking structural signals on top of the RF
- Watch for overfitting — the training set is only 530 variants
