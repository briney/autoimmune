# CR9114 Insights

## Current best models

### Primary (best Spearman): `rf_esm2_fine_94`
- Pipeline: RF=0.94 × RF(bits) + 0.06 × RF(bits+ESM2_PCA50)
- Eval: ρ=0.9596, top-11=0.55, pairwise=0.72 (29/40)

### Best overall (balanced): `fine2_77_8_15` ← **recommended for test**
- Pipeline: RF=0.77 × RF(bits) + 0.08 × GBM_deep + 0.15 × RF(bits+CurrAb_PCA16)
- Eval: ρ=0.9585, top-11=0.64, pairwise=0.88 (35/40)
- Confirmed stable: 6 configurations in (±5% radius) all hit pairwise=0.88

## Feature importances (RF on bits)
| Bit | HC Pos | Mutation | Importance |
|-----|--------|----------|------------|
| 10  | 75     | S→F      | 42.3%      |
| 9   | 74     | K→I      | 16.5%      |
| 4   | 52     | I→S      | 16.3%      |
| 6   | 58     | A→T      | 8.2%       |
| 1   | 29     | F→S      | 6.4%       |
| 7   | 59     | N→A      | 3.2%       |
| others | — | — | <2.2% each |

Bits 10, 9, 4 account for ~75% of signal. All in FR3/CDR2.

## Key findings across iterations

### Models tried and their role
- **RF on bits** (core): ρ=0.9591, CV=0.949. Near-optimal nonparametric fit over 16-bit space.
- **Ridge regression**: Much worse (ρ=0.914 main, 0.921 degree-2). Confirms high-order epistasis.
- **GBM (shallow, depth=3)**: ρ=0.927 alone; improves pairwise when blended (weights 5–15%).
- **GBM (deep, depth=5)**: ρ=0.938 alone but top-11=0.64. Best pairwise component.
- **Multi-seed RF ensemble**: No improvement. RF variance already low at 600 trees.
- **Stacked meta-model (OOF)**: Worse. OOF noise on 530 samples degrades signal.

### Embedding models tried
| Model | Dim | CV (bits+model) | Key contribution |
|-------|-----|-----------------|-----------------|
| ESM2 (650M) | 1280 | 0.917 | +Spearman in small blend; mean-pool better than per-residue |
| AbLang2 (45M) | 480 | 0.927 | top-11=0.64 AND pairwise=0.82 in one model |
| CurrAb (650M) | 1280 | 0.921 | pairwise=0.85 in standalone blend; best for local discrimination |
| ESM1b (650M) | 1280 | 0.903 | No unique contribution over ESM2 |
| ESM2 per-residue | 16×1280 | 0.906-0.919 | Worse than mean-pool (too noisy) |

**Key insight**: Embedding deltas (variant − mature) are very low-dimensional (~99% variance in first 16 PCs for a 1280-dim model). This reflects the small number of mutations (6-8/121 residues).

### ddG tools analysis
- `baddg` and `stabddg` on 10-variant calibration subset:
  - Both tools correlate with truth (rho≈0.84-0.89) but poorly with RF residuals (r≈0.18-0.22)
  - Not worth running at scale — marginal information gain over RF
  - Physics tools may not differentiate within the narrow 6-8 mutation window well

### Why the RF is so good
- 16-bit input, 530 samples → very well-specified problem for tree-based models
- RF with 600 trees essentially computes nonparametric E[affinity | genotype]
- Dominant bits (10, 9, 4) explain ~75% variance; RF handles the nonlinearity perfectly
- The remaining epistatic effects are captured by deep trees

## Pareto frontier (eval set)
| Config | Spearman | Top-11 | Pairwise |
|--------|---------|--------|---------|
| rf_esm2_fine_94 | **0.9596** | 0.55 | 0.72 |
| rf_all3emb_90_4_3_3 | 0.9593 | 0.55 | 0.72 |
| rfab_gbmd_85_5_10 | 0.9592 | 0.55 | 0.80 |
| rf_cr_gbmd_87_5_8 | 0.9589 | 0.55 | 0.80 |
| 4way_1b_5 | 0.9587 | 0.64 | 0.85 |
| fine2_77_8_15 | 0.9585 | **0.64** | **0.88** |

Spearman differences < 0.002 are within noise (95% CI ≈ ±0.018 for n=115).
`fine2_77_8_15` dominates on secondary metrics with negligible Spearman loss.

## Hypothesis for what pairwise=0.88 means
- 35/40 single-mutation pairs correctly ranked
- Missing 5 pairs likely involve context-dependent effects of low-importance bits
- CurrAb16 captures antibody-specific local context that helps discriminate these

## TEST RESULTS (held-out)
- Model: `predictions_final_best.csv` (RF=0.77, GBMdeep=0.08, CurrAb16RF=0.15)
- **Spearman ρ = 0.9524** (eval was 0.9585 — drop of 0.006, within noise)
- Top-11 precision = 0.55 (eval was 0.64 — some drop, may reflect eval tuning bias)
- Pairwise accuracy = **0.89** (31/35) — better than eval's 0.88 (35/40)!
- Excellent generalization. The model is robust.

## Open questions (deprioritized given convergence)
1. CurrAb per-residue at mutation sites — might add marginal Spearman signal
2. Larger RF (1000+ trees) — unlikely to improve beyond noise
3. Neural network MLP on bits — might capture different epistasis patterns
4. Running ddG tools on all variants — low ROI based on calibration data
