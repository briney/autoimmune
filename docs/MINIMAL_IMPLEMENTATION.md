# Minimal Implementation Plan: CR9114 Proof of Concept

Version 0.1 | April 2026

---

## 1. Philosophy

The agent IS the framework. `PROGRAM.md` is the entire orchestration layer —
a markdown file that Claude Code reads and follows autonomously, modeled on
Karpathy's [autoresearch](https://github.com/karpathy/autoresearch) pattern.

Build only what the agent literally cannot do for itself:

1. **Split data with hidden labels** — the agent must not see eval/test ground truth.
2. **Evaluate predictions against hidden ground truth** — the agent submits
   predictions, receives only aggregate metrics.

Everything else — tool invocation, result tracking, analysis scripts, insight
logging, deciding what to try next — is the agent's job. No query API, no
credit system, no results cache framework, no CLI wrapper, no config system,
no tool registry, no metric registry, no iteration counter, no termination
conditions. The agent calls `autobio` directly and manages its own workspace.

---

## 2. Components

The entire implementation is **two Python scripts and one markdown file**.

### 2.1 `prepare.py`

Run once before the experiment starts. Deterministic (fixed random seed).

**Input:** `data/cr9114_h1_binding_data.csv` (65,536 rows)

**Output:** Seven CSV files in `splits/`:

| File | Rows | Columns | Agent access |
|------|------|---------|--------------|
| `train.csv` | ~530 | `genotype`, `h1_mean` | Full read |
| `eval_genotypes.csv` | ~115 | `genotype` | Full read |
| `eval_truth.csv` | ~115 | `genotype`, `h1_mean` | **DO NOT READ** |
| `eval_pairs.csv` | ~30–35 | `genotype_a`, `genotype_b` | **DO NOT READ** |
| `test_genotypes.csv` | ~115 | `genotype` | **DO NOT READ** until final eval |
| `test_truth.csv` | ~115 | `genotype`, `h1_mean` | **DO NOT READ** |
| `test_pairs.csv` | ~25–30 | `genotype_a`, `genotype_b` | **DO NOT READ** |

Sizes are approximate because high-impact pair injection (phase 2 below) adds
a variable number of variants on top of the fixed base split.

**Splitting strategy (two phases):**

*Phase 1 — Stratified random sampling:*

1. Drop any rows with missing `h1_mean`.
2. Bin variants into four affinity strata by `h1_mean` quartiles.
3. Within each stratum, randomly sample 125 train, 25 eval, 25 test
   (500 / 100 / 100 total).
4. Verify that all 16 mutation positions appear in both somatic and germline
   states in every partition.

*Phase 2 — High-impact pair injection:*

Random sampling produces relatively few cross-boundary Hamming-1 pairs
(~8–16), and most of those have modest affinity differences. To ensure
the pairwise accuracy metric tests something meaningful, we supplement
the splits with "small change, big effect" pairs: Hamming-1 neighbors
whose binding affinities differ by ≥10-fold (|Δ(−log₁₀ KD)| ≥ 1.0).

5. Scan the full 65K dataset for all Hamming-1 pairs with |Δh1_mean| ≥ 1.0.
   This yields ~28,000 candidate pairs. Deduplicate by requiring the
   lexicographically smaller genotype to appear first, so each pair is
   counted once.
6. Select 15 pairs for eval and 15 for test using a diversity-aware
   strategy: first pick one pair per mutation position (highest delta
   first), then fill remaining slots regardless of position. This ensures
   the pairwise metric covers the full range of mutation positions rather
   than concentrating on the few highest-leverage positions. Pairs are
   also filtered to avoid sharing genotypes between selected pairs.
7. For each selected pair, randomly assign one member to train and the
   other to eval (or test). Both members are new variants not present
   in the phase-1 splits.
8. Re-verify uniqueness — no genotype appears in more than one partition.
9. Find ALL cross-boundary Hamming-1 pairs (both naturally occurring from
   phase 1 and injected in phase 2). Write to `eval_pairs.csv` and
   `test_pairs.csv` for use by `evaluate.py`.

The result: ~530 train, ~115 eval, ~115 test, with ~30–35 eval pairs and
~25–30 test pairs, roughly double what the random split produces alone.

Dependencies: `pandas`, `numpy`. No other libraries.

### 2.2 `evaluate.py`

The agent's only interface to ground truth. Called from the command line:

```bash
python evaluate.py predictions.csv [--partition eval|test]
```

**Input:** A CSV with columns `genotype` and `predicted_score` (higher =
tighter predicted binding). The agent writes this file in whatever way it
sees fit.

**Output (printed to stdout):**

```
=== Evaluation (eval, 100 variants) ===
Spearman rho:      0.42
Top-10 precision:  0.30
Pairwise accuracy: 0.68 (23/34 pairs correct)
```

**Metrics:**

| Metric | Definition |
|--------|------------|
| Spearman ρ | Rank correlation between `predicted_score` and true `h1_mean` |
| Top-k precision | Of the top k variants by predicted score, what fraction are in the true top k? Reported for k ≈ 10% of the partition size. |
| Pairwise accuracy | For each single-mutation neighbor pair in `eval_pairs.csv` where both members appear in predictions: did the predicted ranking match the true ranking? |

**Behavior:**

- Reads truth from `splits/eval_truth.csv` (default) or `splits/test_truth.csv`
  (when `--partition test`).
- If predictions are missing genotypes present in the truth set, warns and
  evaluates on the intersection.
- Exits with code 0 on success, 1 on error.
- Prints only the metrics above. No individual predictions or KD values are
  revealed.

Dependencies: `scipy` (for `spearmanr`), `pandas`. No other libraries.

~60–80 lines of code.

### 2.3 `PROGRAM.md`

The agent's complete instruction sheet. Claude Code reads this file and
follows it autonomously. This is the most important component — it replaces
all the orchestration code that a traditional framework would provide.

**Structure:**

```
1. YOUR TASK
   - Discover a pipeline of autobio tools that predicts relative binding
     affinity of CR9114 variants for H1.
   - Performance measured by Spearman ρ on the held-out eval set.

2. THE SYSTEM
   - CR9114 antibody, 16-position binary genotype, H1 antigen.
   - Mutation key (inline or reference to file).
   - Mature CR9114-H1 structure available at structures/cr9114_mature_h1.pdb.

3. DATA
   - train.csv: 500 variants with genotypes and h1_mean (in splits/).
     Read this file — it fits in context.
   - eval_genotypes.csv: 100 genotypes to predict (in splits/).
   - Sequences: cr9114_h1_sequences.fasta for building tool inputs.
   - Mutation key: cr9114_mutation_key.csv for decoding genotypes.

4. TOOLS
   - autobio CLI. Run `autobio list` to see tools, `autobio info <tool>`
     for details.
   - Relevant tools for this task:
     * evoef2_build_mutant — introduce mutations into a PDB
     * rosetta_score, evoef2_binding — score structures
     * stabddg, baddg — ML-based binding ΔΔG prediction
     * rosetta_minimize, openmm_amber_minimize — energy minimization
     * openmm_amber_relax — molecular relaxation
     * esm2, esm1b, ablang2, currab — sequence embeddings
   - Invoke: `autobio run <tool> --config config.json [--format json]`
   - Query: `autobio info <tool> --format json`

5. IN SILICO MUTAGENESIS
   - The reference structure is the mature (fully somatic) CR9114-H1 complex.
   - To model a variant, determine which positions revert to germline
     (0 bits in the genotype) and introduce those mutations using
     evoef2_build_mutant.
   - Example: genotype 1111111111111110 means position 16 (Y106) reverts
     from somatic Ser back to germline Tyr → mutation is S106Y.
   - After building the mutant, optionally minimize/relax before scoring.

6. EVALUATION
   - Write predictions to a CSV: genotype,predicted_score
   - Run: python evaluate.py <predictions.csv>
   - You will see Spearman ρ, top-10 precision, and pairwise accuracy.
   - You will NOT see individual KD values for the eval set.

7. ITERATION PROTOCOL
   a. Orient — read train.csv, review your prior results and notes.
   b. Hypothesize — state what you want to try and why.
   c. Execute — run autobio tools, write analysis scripts, generate predictions.
   d. Evaluate — run evaluate.py, record metrics.
   e. Reflect — update INSIGHTS.md with what you learned.
   f. Loop — go to (a). Never stop. Run until interrupted by the operator.

8. RECORDING RESULTS
   - Maintain workspace/results.tsv with columns:
     iteration, pipeline_description, spearman_rho, top10_precision,
     pairwise_accuracy, notes
   - Maintain workspace/INSIGHTS.md as a running scratchpad of what
     you've learned about the problem.
   - You may create any other files you need in workspace/.

9. RULES
   - DO NOT read eval_truth.csv, test_truth.csv, test_genotypes.csv,
     eval_pairs.csv, or test_pairs.csv.
   - DO NOT modify prepare.py, evaluate.py, or anything in data/.
   - You may read and write anything in workspace/.
   - You may write and execute Python scripts for analysis (numpy, scipy,
     scikit-learn, pandas are available).
   - You may use any autobio tool.
   - Commit meaningful progress to git as you go.

10. STRATEGIC HINTS
    - Start simple. Try a single cheap scoring tool on a small subset before
      building complex pipelines.
    - The training data is small enough to read in full. Study it — which
      positions correlate most with affinity?
    - Sequence-based approaches (embeddings + regression) and structure-based
      approaches (scoring after mutagenesis) may complement each other.
    - The affinity range is only ~3 orders of magnitude. Signal is subtle.
    - Consider ensemble methods that combine multiple tool outputs.
```

The exact wording will be refined during implementation, but this is the
structural skeleton.

---

## 3. Directory Layout

```
experiments/CR9114/
├── data/                              # already exists
│   ├── cr9114_h1_binding_data.csv     # full 65K dataset
│   ├── cr9114_mutation_key.csv        # genotype bit → AA position/residue
│   ├── cr9114_h1_sequences.fasta      # mature heavy, light, H1 sequences
│   └── README.md
├── structures/
│   └── cr9114_mature_h1.pdb           # pre-computed Boltz-2 structure (TO PLACE)
├── splits/                            # generated by prepare.py
│   ├── train.csv                      # ~530 rows (500 base + injected pair members)
│   ├── eval_genotypes.csv             # ~115 rows
│   ├── eval_truth.csv
│   ├── eval_pairs.csv                 # ~30–35 cross-boundary Hamming-1 pairs
│   ├── test_genotypes.csv
│   ├── test_truth.csv
│   └── test_pairs.csv
├── prepare.py                         # one-time data splitting
├── evaluate.py                        # blind evaluation
├── PROGRAM.md                         # agent instruction sheet
└── workspace/                         # agent's scratch space (initially empty)
    ├── results.tsv                    # created by agent
    ├── INSIGHTS.md                    # created by agent
    └── ...                            # whatever the agent needs
```

---

## 4. Prerequisites

Before the experiment can run:

1. **Place the pre-computed structure.** Copy the mature CR9114–H1 PDB into
   `experiments/CR9114/structures/cr9114_mature_h1.pdb`.

2. **Install Python dependencies.** The two scripts need only `pandas`,
   `numpy`, and `scipy`. These should already be available in the `ai`
   conda environment. If not: `pip install pandas numpy scipy`.

3. **Verify autobio.** Run `autobio list` to confirm tools are accessible.
   Run `autobio images` to confirm Docker images are cached for the tools
   the agent is likely to use first (rosetta_score, evoef2_build_mutant,
   evoef2_binding, esm2).

4. **Run prepare.py.** `cd experiments/CR9114 && python prepare.py`. Verify
   the split files look correct (row counts, column names, no data leakage).

---

## 5. Implementation Sequence

| Step | What | Est. effort |
|------|------|-------------|
| 1 | Write `prepare.py` | ~1 hour |
| 2 | Write `evaluate.py` | ~1 hour |
| 3 | Write `PROGRAM.md` | ~2 hours |
| 4 | Place structure, run `prepare.py`, smoke-test `evaluate.py` with dummy predictions | ~30 min |
| 5 | Launch Claude Code with PROGRAM.md, monitor first 2–3 iterations | ~1–2 hours |

Total: roughly half a day from blank slate to running experiment.

---

## 6. PDB Mutagenesis

The CR9114 pilot requires in silico mutagenesis: given the mature CR9114–H1
structure, introduce germline reversions to model each variant. `autobio`
already provides `evoef2_build_mutant`, which takes a PDB and a list of
mutations, swaps residues, optimizes local rotamers via EvoEF2, and outputs
a new PDB. This covers the PoC need.

A leaner standalone mutagenesis tool (residue swap + best rotamer placement,
no energy function coupling) would be a good future addition to autobio for
general use, but building it before the PoC reveals whether `evoef2_build_mutant`
is a bottleneck would be premature.

---

## 7. What Success Looks Like

The PoC succeeds if:

1. The agent autonomously runs multiple iterations without human intervention.
2. The agent discovers at least one pipeline that outperforms any single tool
   in isolation (measured by eval Spearman ρ).
3. The total infrastructure code (prepare.py + evaluate.py) stays under 200
   lines combined.
4. PROGRAM.md is sufficient to drive the experiment — no additional
   orchestration code is needed.

The PoC does NOT need to achieve a specific Spearman ρ threshold. It is a
proof of concept for the autonomous pipeline discovery paradigm, not a
benchmark for antibody affinity prediction.

---

## 8. Future Directions (Out of Scope for PoC)

These are explicitly deferred:

- Scaling to the full 65K dataset (requires query API, credit system)
- Multiple antigens (H3, Flu-B)
- Generalized experiment framework (`src/autoimmune/` package code)
- Germline reference structure
- Formal iteration tracking / Tier 2 summaries
- Termination conditions
