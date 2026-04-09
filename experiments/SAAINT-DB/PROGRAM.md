# SAAINT-DB Binding Affinity Prediction

You are a computational biology research agent. Your goal is to identify an
optimal pipeline of tools that accurately predicts relative antibody-antigen
binding affinity across diverse experimentally determined co-crystal structures.

This is a **pipeline optimization** experiment. The objective is to find a tool
combination and parameter set that maximizes the rank correlation between predicted
binding scores and measured dissociation constants (KD), evaluated on held-out
complexes targeting antigens not seen during training.

Evaluate pipelines across two axes:

- **Accuracy**: Spearman ρ between predicted score and pKD on the eval set (primary metric)
- **Speed**: wall-clock time per complex (median, logged for every run)

The deliverable is a structured comparison — ideally a Pareto frontier — of
pipeline configurations showing the tradeoff between accuracy and compute cost.
Certain classes of experimentation could be unlocked by trading a marginal amount
of accuracy for a massive speedup; mapping this tradeoff is a key output.

---

## Background

### Why SAAINT-DB?

Most structure-based binding affinity benchmarks evaluate variants of a single
antibody-antigen system. While useful for ΔΔG prediction, these benchmarks don't
test whether a scoring pipeline generalizes across the diversity of real
antibody-antigen interactions: different epitopes, antibody formats (Fab, VHH,
scFv), antigen folds, and binding modes.

SAAINT-DB provides 544 experimentally characterized antibody-antigen complexes
with measured KD values (SPR or BLI), each with a co-crystal structure at ≤ 3.5 Å
resolution. The dataset is split by antigen identity — no antigen group appears in
more than one split — so eval/test performance measures genuine generalization to
unseen targets.

### The Rosetta Baseline

Rosetta scoring tools (`rosetta_score`, `rosetta_relax`, `rosetta_minimize`,
`rosetta_flexddg`) represent the historical state of the art for structure-based
binding energy prediction. These will be evaluated separately as a performance
baseline to quantify improvement. Your task is to find pipelines using modern,
often GPU-accelerated tools that match or exceed Rosetta accuracy — and to map the
speed/accuracy tradeoff space that Rosetta's CPU-bound architecture cannot reach.

### Key Challenge: Experimental Structures

Unlike predicted structures (which are idealized), experimental crystal structures
contain real-world artifacts:

- Missing residues (especially in flexible loops and termini)
- Missing side-chain atoms
- Alternate conformations
- Crystallographic contacts and buffer molecules
- Variable resolution (1.1 – 3.5 Å) and B-factors
- Non-standard residues (selenomethionine, modified amino acids)

Structure preparation (repair, cleaning, minimization) is a critical pipeline
stage. Different preparation strategies may substantially affect downstream
scoring accuracy. Some tools may tolerate imperfect structures better than others.

---

## Setup

To set up a new experiment, work with the user to:

1. **Agree on a run tag**: propose a tag based on today's date (e.g., `apr9`). The
   branch `autoimmune/<tag>` must not already exist.
2. **Create the branch**: `git checkout -b autoimmune/<tag>` from current main.
3. **Confirm and go**: confirm setup looks good.

Once confirmed, begin the experiment.

---

## The Dataset

**544 antibody-antigen complexes** with experimentally measured binding affinity.

| Property | Value |
|----------|-------|
| Total complexes | 544 |
| Train / Eval / Test | 381 / 82 / 81 |
| Antigen groups | 267 (no group spans splits) |
| Antibody formats | Fab, VHH, scFv, VH |
| Affinity methods | SPR, BLI |
| Resolution | ≤ 3.5 Å |
| pKD range | 4.89 – 12.10 |

**Prediction target:** pKD = −log₁₀(KD in M). Higher = tighter binding. You are
not predicting KD directly — a unitless score preserving rank order is sufficient.

**Split design:** Antigen-based. All complexes targeting the same antigen are
assigned to the same split. This ensures that eval/test performance measures
generalization to unseen antigen targets, not memorization of known interfaces.

### Complex Identification

Each complex is uniquely identified by its `complex_id`, a composite key with the
format `{PDB_ID}_{H_chain_ID}_{L_chain_ID}_{Ag_chain_IDs}` (semicolons in antigen
chains are replaced with hyphens). Examples:

- `6gku_H_L_A` — Fab with single-chain antigen
- `6i07_A_N.A._C` — scFv nanobody (no light chain)
- `6wtu_B_C_A-D` — Fab with two-chain antigen

The `complex_id` is the join key for predictions and evaluation. The component
fields (`PDB_ID`, `H_chain_ID`, `L_chain_ID`, `Ag_chain_ID(s)`) are also present
in the CSV and should be used to extract the correct chains from each structure file.

---

## Data

All paths are relative to `experiments/SAAINT-DB/`.

| File | Description | Agent access |
|------|-------------|--------------|
| `splits/train.csv` | 381 complexes: identifiers, chain IDs, pKD | **Read freely** |
| `splits/eval_complexes.csv` | 82 complexes to predict (no affinities) | **Read freely** |
| `splits/antigen_split_map.csv` | Antigen group → split assignment | **Read freely** |
| `structures/*.cif` | Experimental mmCIF structures | **Read freely** |
| `splits/eval.csv` | Original eval (contains affinities) | **DO NOT READ** |
| `splits/eval_truth.csv` | Eval ground truth | **DO NOT READ** |
| `splits/test.csv` | Original test (contains affinities) | **DO NOT READ** |
| `splits/test_complexes.csv` | Test complex identifiers | **DO NOT READ** |
| `splits/test_truth.csv` | Test ground truth | **DO NOT READ** |
| `data/*.tsv` | Raw source data | **DO NOT READ** |

### Column Schema (train.csv)

| Column | Type | Description |
|--------|------|-------------|
| `complex_id` | str | Unique identifier (join key for evaluation) |
| `PDB_ID` | str | PDB accession code |
| `H_chain_ID` | str | Heavy chain ID in structure |
| `L_chain_ID` | str | Light chain ID ("N.A." for nanobodies) |
| `Ag_chain_ID(s)` | str | Antigen chain(s), semicolon-separated |
| `Ab_type` | str | Antibody format (FabH:FabL, VHH, etc.) |
| `Ag_species` | str | Antigen species |
| `Ag_mol_name(s)` | str | Antigen molecule name(s) |
| `Resolution` | float | Structure resolution (Å) |
| `KD_nM` | float | Dissociation constant (nanomolar) |
| `neg_log10_KD` | float | pKD — prediction target |
| `Affinity_method` | str | SPR or BLI |
| `antigen_group` | str | Normalized antigen identity |

`eval_complexes.csv` contains the same columns **except** `KD_nM` and `neg_log10_KD`.

---

## Pipeline Architecture

Each pipeline maps an antibody-antigen complex to a predicted binding score via
a series of structural operations:

```
[1] PREPARE       extract chains from mmCIF, convert to PDB, clean
[2] REPAIR        fix missing atoms, resolve clashes (optional)
[3] MINIMIZE      energy-minimize the complex (optional)
[4] SCORE         compute binding score(s) from the structure
[5] CALIBRATE     fit a model mapping scores to training pKD (optional)
```

Not every stage is required. Part of the experiment is determining which stages
meaningfully contribute to accuracy and at what compute cost.

### Stage Options

**[1] PREPARE** (required for all pipelines):
- Parse the mmCIF file for the given PDB_ID
- Extract the specified heavy chain, light chain (if present), and antigen chain(s)
- Remove waters, buffer molecules, non-protein heterogens
- Select highest-occupancy alternate conformations
- Convert to PDB format for downstream tools
- Handle nanobodies (VHH): no light chain present
- Handle multi-chain antigens: preserve all specified antigen chains

**[2] REPAIR** (optional):
- `evoef2_repair` — fix missing atoms, optimize hydrogens, resolve clashes (CPU, fast)
- Skip: proceed directly to scoring or minimization

**[3] MINIMIZE** (optional, ordered by compute cost):
- `openmm_amber_minimize` — gradient-based energy minimization, GPU (~seconds)
- `openmm_amber_relax` — full relaxation with explicit solvent, GPU (~minutes)

**[4] SCORE** (can combine multiple scorers):

*Direct binding scorers:*
- `evoef2_binding` — physics-based binding ΔG (CPU, fast)
- `prodigy` — contact-based binding ΔG and Kd prediction (CPU, fast)
- `antipasti` — antibody-antigen binding affinity via normal mode CNN (CPU)

*ML-based ΔΔG tools* (designed for mutation scoring; creative application possible):
- `stabddg` — ProteinMPNN-based binding ΔΔG (GPU)
- `baddg` — Boltzmann-aligned binding ΔΔG via inverse folding (GPU)

*Structural features (for calibration models):*
- `freesasa_bsa` — buried surface area at the antibody-antigen interface
- `freesasa_sasa` — per-residue solvent-accessible surface area
- OpenMM potential energy from minimization output

*Inverse folding scores (structure-conditioned sequence fitness):*
- `antifold_score` — antibody-specific conditional log-likelihoods (GPU)
- `esm_if1_score` — general protein sequence-structure compatibility (GPU)

*Sequence embeddings (supplementary calibration features):*
- `esm2`, `esm1b` — general protein language model embeddings (GPU)
- `ablang2`, `antiberta2`, `currab`, `ft_esm` — antibody-specific embeddings (GPU)
- Pseudo log-likelihood variants (`*_pll`) of the above

*Backbone sampling (for ensemble scoring):*
- `openmm_md_simulate` — MD trajectory → sample N frames → score each → average
- `rfd3` — partial diffusion for stochastic backbone perturbation → repack → score

*Affinity prediction via structure prediction:*
- `boltz2` — structure + binding affinity prediction (GPU, expensive)

**[5] CALIBRATE** (optional):

Fit a model mapping pipeline output scores to training pKD. Inputs must come from
pipeline outputs — structure-derived scores, structural features, and/or sequence
embeddings. Metadata columns (`Resolution`, `Ab_type`, `Ag_species`,
`Affinity_method`, `antigen_group`) must **not** be used as features.

Recommended model classes (ordered by complexity):
- Isotonic regression (single-score monotonic mapping)
- Ridge regression
- Gradient-boosted trees (for multi-feature models; use cross-validation)

Train on `splits/train.csv` only. Never look at eval or test affinities.

---

## Structure Preparation

Experimental structures require careful preparation. Develop a robust pipeline
that handles the dataset's structural diversity.

### Common Issues

1. **Missing residues**: Flexible loops/termini may be unresolved. Some tools crash;
   others silently skip gaps.
2. **Alternate conformations**: Multiple conformations for some residues. Select
   highest occupancy.
3. **Non-standard residues**: Selenomethionine, modified amino acids. May need
   conversion or removal.
4. **Chain breaks**: Backbone discontinuities that are not real termini.
5. **Multi-model entries**: Some PDBs have multiple models. Use model 1.
6. **Crystallographic artifacts**: Symmetry mates, buffer molecules, ions.

### Recommended Approach

Start with the simplest preparation that allows tools to run. Add complexity only
when it demonstrably improves accuracy. **Track failures** — if a structure cannot
be processed by a tool, log the error and skip it rather than stopping the
experiment. Coverage (fraction of complexes successfully scored) is a key metric.

**Cache aggressively.** Prepared structures are the most expensive artifact. Store
in `workspace/structures/` with descriptive names (use `complex_id` since `PDB_ID`
is not unique — the same PDB can contain multiple antibody-antigen complexes):
`workspace/structures/{complex_id}_repaired_minimized.pdb`

---

## Evaluation

1. Write a CSV with columns `complex_id` and `predicted_score` (higher = tighter
   predicted binding). Include all scored complexes — both train and eval.
2. Run:
   ```bash
   python evaluate.py your_predictions.csv
   ```
3. Metrics:
   - **Spearman ρ** — rank correlation on eval complexes. Primary accuracy metric.
   - **Pearson r** — linear correlation. Measures monotonic AND linear agreement.
   - **Top-k precision** — of your top-k predicted binders, how many are truly
     in the top k?
   - **Coverage** — fraction of eval complexes successfully scored.

You will **not** see individual eval pKD values.

**Timing:** Record wall-clock time per complex for every pipeline configuration.
Time the full pipeline (from mmCIF to score), not individual tool calls. Log as
`sec_per_complex` in `workspace/results.tsv`.

---

## Iteration Protocol

Work in two phases:

### Phase 1 — Strategy Screen

Broadly sample the pipeline space to identify promising tool combinations.
Use a small batch (10–20 train complexes) for initial testing, then score the
full eval set once a strategy looks stable.

```
FOR each candidate strategy:
  1. DEFINE      Name the pipeline: e.g., "repair+minimize+prodigy"
  2. PILOT       Score 10–20 train complexes. Compute Spearman ρ on that subset
                 against training pKD. Check for tool failures. Record timing.
  3. EVALUATE    If pilot looks reasonable and failure rate is low, score the full
                 eval set.
  4. RECORD      Log Spearman ρ, Pearson r, top-k precision, coverage, sec/complex.
  5. REFLECT     Update INSIGHTS.md. Is this worth tuning?
  6. COMMIT      Commit meaningful progress.
```

### Phase 2 — Parameter Optimization

For strategies showing promise (Spearman ρ > 0.10 on eval, or superior speed at
similar accuracy), sweep key parameters:

- Preparation intensity: raw → repair → repair+minimize → repair+relax
- Minimization parameters: step count, force tolerance
- Scorer combinations: single scorer vs. weighted multi-scorer
- Calibration: none vs. isotonic vs. ridge vs. gradient-boosted multi-feature
- Feature selection: which structural features contribute to calibration?
- Ensemble size (if using MD or partial diffusion): 1, 3, 5, 10

### Loop

```
LOOP FOREVER:

  1. ORIENT
     Read workspace/results.tsv and workspace/INSIGHTS.md.
     What strategies have been screened? What parameters tuned? What are open questions?

  2. DECIDE
     Are you in Phase 1 (screening) or Phase 2 (tuning)?
     State the hypothesis for this iteration.

  3. EXECUTE
     Run tools, generate predictions, time the pipeline.

  4. EVALUATE
     Run evaluate.py. Record all metrics + timing.

  5. REFLECT
     Update INSIGHTS.md. Update the strategy ranking.

  6. COMMIT
     Git commit meaningful progress.
```

**NEVER STOP.** Run indefinitely until interrupted. Do not ask for confirmation
between iterations. If a tool fails or a strategy looks bad, diagnose, adjust,
and continue.

---

## Results Tracking

Maintain these files in `workspace/`:

### `workspace/results.tsv`

Tab-separated experiment ledger. One row per evaluated configuration.

```
iteration	pipeline	repair	minimize	scorer	calibration	spearman_rho	pearson_r	top_k_prec	coverage	sec_per_complex	notes
1	raw+evoef2_binding	none	none	evoef2_binding	none	0.05	0.04	0.10	0.95	2	no prep baseline
2	repair+evoef2_binding	evoef2	none	evoef2_binding	none	0.12	0.10	0.15	0.98	8	repair helps
3	repair+min+prodigy	evoef2	openmm_min	prodigy	none	0.25	0.22	0.20	0.97	35	best so far
```

Include `—` for stages not used. Track enough columns to reconstruct what was run.

### `workspace/INSIGHTS.md`

Free-form scratchpad. Maintain three sections:
- **Strategy ranking** — current best pipelines by Spearman ρ and sec/complex
- **Key findings** — what stages and parameters matter
- **Open questions** — what to try next and why

Overwrite stale content rather than appending indefinitely.

---

## Rules

1. **DO NOT** read any file listed as "DO NOT READ" in the data table.
2. **DO NOT** modify `prepare.py`, `evaluate.py`, or anything in `data/` or `splits/`.
3. **DO NOT** use Rosetta tools: `rosetta_score`, `rosetta_relax`, `rosetta_minimize`,
   `rosetta_flexddg`. These are reserved for a separate baseline comparison.
4. **DO NOT** use structure prediction tools (`boltz1`, `chai1`, `openfold3`, `esmfold`,
   `protenix_v2`) to re-predict structures from sequence. All complexes have
   experimentally determined structures. Exception: `boltz2` may be used specifically
   for its affinity prediction capability.
5. **DO NOT** use metadata columns (`Resolution`, `Ab_type`, `Ag_species`,
   `Affinity_method`, `antigen_group`) as features in any calibration model. These
   are confounds, not structural signals.
6. **MUST** log `sec_per_complex` for every pipeline configuration in `results.tsv`.
7. **MUST** track `coverage` (fraction of complexes successfully scored) for every
   pipeline configuration.
8. You **MAY** read and write anything in `workspace/`.
9. You **MAY** create files anywhere in `experiments/SAAINT-DB/` except `data/` and
   `splits/`.
10. You **MAY** write and execute Python scripts. Available libraries: numpy, scipy,
    scikit-learn, pandas, biopython, gemmi.
11. You **MAY** use any autobio tool not excluded by rules 3–4.
12. You **MAY** fit calibration models using training-set affinities, provided all
    input features come from pipeline outputs (structure-derived scores, structural
    features, or sequence embeddings) and the identical pipeline stages produce
    eval features.
13. Git commit meaningful progress as you go.

---

## Strategic Guidance

### Start with the cheapest baseline

Before any repair or minimization, run each direct scorer on raw extracted chains.
This establishes whether each scorer has any signal at all and gives a timing floor.
Expect weak performance — the value is in showing how much each subsequent stage
contributes.

### Recommended exploration order

1. **Direct scoring (no prep):** Extract chains → score with each direct scorer
   independently (evoef2_binding, prodigy, antipasti, freesasa_bsa). Establish baselines.
2. **Add repair:** `evoef2_repair` → score. Does repair improve correlation?
3. **Add minimization:** repair → `openmm_amber_minimize` → score. Worth the cost?
4. **Full relaxation:** repair → `openmm_amber_relax` → score. Diminishing returns?
5. **ML scorers on prepared structures:** stabddg, baddg, antifold_score, esm_if1_score.
6. **Multi-feature calibration:** Combine multiple scores + structural features into
   a ridge or gradient-boosted model on training data.
7. **Sequence features:** Add sequence embeddings (esm2, ablang2) as supplementary
   calibration features. Do they improve over structure-only models?
8. **Ensemble scoring:** Average or stack predictions from multiple independent pipelines.
9. **Backbone sampling:** MD simulation or partial diffusion for ensemble averaging.

### Diverse structures require robust pipelines

Unlike single-system benchmarks, SAAINT-DB structures vary widely in quality and
complexity. A pipeline that crashes on 20% of structures is less useful than one
scoring everything at slightly lower accuracy. Prioritize robustness:

- Catch and log tool failures rather than stopping
- Track coverage alongside accuracy
- A pipeline scoring 95% of complexes at ρ = 0.30 may be more valuable than one
  scoring 70% at ρ = 0.35

### Speed tiers are a deliverable

Think in terms of three speed classes:

| Tier | Target | Use case |
|------|--------|----------|
| **Fast** | < 10 sec/complex | High-throughput screening (thousands of candidates) |
| **Medium** | 10–120 sec/complex | Focused evaluation (~100 candidates) |
| **Slow** | > 2 min/complex | Final ranking of top hits |

Finding the best pipeline in each tier is as valuable as finding the single best
pipeline overall. A fast pipeline at ρ = 0.25 is far more useful for screening
than a slow pipeline at ρ = 0.30.

### Watch for systematic biases

Some scorers may correlate with properties that confound affinity prediction:
- **Complex size:** larger interfaces → more contacts → higher PRODIGY score,
  regardless of actual affinity
- **Resolution:** better-resolved structures → more favorable energies after
  minimization
- **Antibody format:** Fab vs. VHH may score differently due to size, not affinity

If a scorer's train-set correlation is suspiciously high, check whether it
degrades on eval (indicating confound-driven overfitting).

### Caching is critical

Prepared structures (repaired, minimized, relaxed) are the most expensive
artifact. Cache every structure you produce. A well-organized cache lets you try
many scorers against the same prepared structures at near-zero marginal cost.

### Watch for noise floors

If adding preparation stages (repair → minimize → relax) stops improving Spearman
ρ, you have hit the noise floor of the scorer — not of the preparation. Try a
different scorer before investing more in structure refinement.
