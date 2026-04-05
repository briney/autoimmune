# CR9114 Autonomous Pipeline Discovery

You are a computational biology research agent. Your goal is to discover a
pipeline — or ensemble — of structural biology and sequence analysis tools
that accurately predicts the relative binding affinity of antibody variants
for an antigen.

No single tool is expected to perform well alone. The value is in discovering
how to combine them.

---

## Setup

To set up a new experiment, work with the user to:

1. **Agree on a run tag**: propose a tag based on today's date (e.g. `mar5`). The 
   branch `autoimmune/<tag>` must not already exist — this is a fresh run.
2. **Create the branch**: `git checkout -b autoimmune/<tag>` from current master.
3. **Confirm and go**: Confirm setup looks good.

Once you get confirmation, kick off the experimentation.

---

## The System

**Antibody:** The antibody heavy chain has 16 positions where the mature 
(affinity-matured) residue differs from the inferred germline residue. 
Each variant is defined by a 16-bit genotype string: `1` = somatic residue, 
`0` = germline residue. The light chain is constant across all variants. 

**Affinities:** Measured as −log₁₀(KD). Higher = tighter binding. Range:
~7.0 (weak, ~100 nM) to ~9.8 (strong, ~150 pM). You do not need to predict
exact KD — a unitless score that preserves the rank ordering is sufficient.

### Dataset Restriction

All variants in the train/eval/test splits have **6, 7, or 8 somatic
mutations** (out of 16 possible). This is a deliberate restriction: in the
full dataset, mutation count correlates with affinity (more mutations →
generally tighter binding). Restricting to a narrow mutation-count window
removes that shortcut. Your pipeline must infer affinity from the actual
sequence and structural properties of each variant — not from how many
mutations it carries.

### Mutation Key

| Bit | HC Pos | Region | Germline → Somatic |
|-----|--------|--------|--------------------|
| 1   | 29     | CDR1   | F → S              |
| 2   | 30     | CDR1   | S → N              |
| 3   | 31     | CDR1   | S → N              |
| 4   | 52     | CDR2   | I → S              |
| 5   | 57     | CDR2   | T → S              |
| 6   | 58     | FR3    | A → T              |
| 7   | 59     | FR3    | N → A              |
| 8   | 71     | FR3    | T → S              |
| 9   | 74     | FR3    | K → I              |
| 10  | 75     | FR3    | S → F              |
| 11  | 76     | FR3    | T → S              |
| 12  | 77     | FR3    | S → N              |
| 13  | 84     | FR3    | S → N              |
| 14  | 87     | FR3    | R → T              |
| 15  | 95     | CDR3   | Y → F              |
| 16  | 106    | FR4    | Y → S              |

Genotype `1111111111111111` = fully mature. `0000000000000000` = fully germline.

---

## Data

All paths are relative to `experiments/CR9114/`.

| File | Description | Agent access |
|------|-------------|--------------|
| `splits/train.csv` | ~530 variants: `genotype`, `h1_mean` | **Read freely** |
| `splits/eval_genotypes.csv` | ~115 genotypes to predict | **Read freely** |
| `data/cr9114_mutation_key.csv` | Bit → HC position → residues | **Read freely** |
| `data/cr9114_h1_sequences.fasta` | Mature heavy, light, antigen sequences | **Read freely** |
| `structures/cr9114_mature_h1.pdb` | Predicted structure of CR9114 bound to its antigen | **Read freely** |
| `splits/eval_truth.csv` | Eval ground truth | **DO NOT READ** |
| `splits/test_truth.csv` | Test ground truth | **DO NOT READ** |
| `splits/test_genotypes.csv` | Test genotypes | **DO NOT READ** |
| `splits/eval_pairs.csv` | Eval pair indices | **DO NOT READ** |
| `splits/test_pairs.csv` | Test pair indices | **DO NOT READ** |
| `data/cr9114_h1_binding_data.csv` | Full 65K dataset | **DO NOT READ** |

The training set (~530 variants) is small enough to fit in context. Read it
in full. Study it before reaching for tools.

---

## Tools

You have the full `autobio` toolkit.

```bash
autobio list                          # see all tools
autobio info <tool> --format json     # input schema, parameters, notes
autobio run <tool> --config cfg.json  # run a tool
```

### Key tools for this task

**Mutagenesis:**
- `evoef2_build_mutant` — introduce mutations into a PDB. Swaps residues and
  optimizes local rotamers.

**Scoring:**
- `rosetta_score` — Rosetta energy function (total, per-residue, interface).
- `evoef2_binding` — EvoEF2 physics-based binding energy.
- `stabddg` — ML-based binding ΔΔG (ProteinMPNN architecture).
- `baddg` — Boltzmann-aligned binding ΔΔG.

**Minimization / Relaxation:**
- `rosetta_minimize` — gradient-based energy minimization.
- `openmm_amber_minimize` — OpenMM + Amber force field.
- `openmm_amber_relax` — full relaxation with backbone flexibility.

**Sequence embeddings:**
- `esm2` — ESM-2 (650M parameters).
- `esm1b` — ESM-1b (650M parameters).
- `ablang2` — antibody-specific (45M parameters).
- `currab` — antibody-specific (650M parameters).

Run `autobio info <tool> --format json` before first use to understand
required inputs and parameters.

---

## In Silico Mutagenesis

To model a variant from its genotype:

1. Start from the mature structure (`structures/cr9114_mature_h1.pdb`).
   The mature genotype is all-1s — every position has the somatic residue.
   In the PDB, the antigen (a trimer) is chains A, B, and C. The CR9114 heavy
   chain is chain D, and the CR9114 light chain is chain E.
2. Identify which bits are `0` — those positions revert to germline.
3. For each `0`-bit at position k, construct a mutation string:
   `{chain}{somatic_1letter}{hc_position}{germline_1letter}`.
   Example: bit 1 = 0 means Ser→Phe at HC position 29 → `"DS29F"`
   (reminder, the heavy chain is chain D; **all mutations in the dataset have been
   performed on the CR9114 heavy chain, so only chain D in the PDB will 
   be mutated during your experiments**).
4. Run `evoef2_build_mutant` with the mutation list.
5. Optionally minimize or relax the result before scoring.

Whether to minimize/relax (and for how many iterations) as well as the 
scoring function(s) to be used are totally up to you. Any allowed tools 
in the `autobio` toolkit are possibilities (including ensembling tools or 
selectively applying one tool to certain mutations and another tool to others) 
or you can even build your own scoring functions from scracth and run them. 
Anything is fair game.

For the fully mature genotype (`1111111111111111`), no mutagenesis is needed —
the reference structure is already that variant.

---

## Evaluation

1. Write a CSV with columns `genotype` and `predicted_score` (higher = tighter
   predicted binding). Include scores for **all** genotypes you have scored —
   both train and eval — to enable the pairwise accuracy metric.
2. Run:
   ```bash
   python evaluate.py your_predictions.csv
   ```
3. You will see three metrics:
   - **Spearman ρ** — rank correlation on eval variants. Primary metric.
   - **Top-10 precision** — of your predicted 10 tightest binders, how many
     are truly in the top 10?
   - **Pairwise accuracy** — for single-mutation variant pairs, did you
     correctly predict which binds tighter?

You will **not** see individual eval KD values.

For final test evaluation (only when you believe you have converged):
```bash
python evaluate.py your_predictions.csv --partition test
```

---

## Iteration Protocol

While there is no hard time limit on each iteration, the priority is on rapid 
iteration rather than trying everything possible before submitting your first 
evaluation. Failure is totally fine; a series of safe, low-risk experiments that 
produce only marginal improvements is much less interesting than high-risk "swing for
the fences" iterations that frequntly perform poorly but occasionally uncover 
paradigm-shifting improvements. 

This is a difficult problem, and one that is likely to be best attacked by 
learning from the feedback recieved across many experimental iterations. It is 
typically considered best practice to test only one experimental variable at a 
time, in order to isolate its effects. 

```
LOOP FOREVER:

  1. ORIENT
     Read workspace/results.tsv and workspace/INSIGHTS.md.
     What have you tried? What worked? What are the open questions?

  2. HYPOTHESIZE
     State what you will try this iteration and why, based on prior results.

  3. EXECUTE
     Run autobio tools, write analysis scripts, generate predictions. You are 
     welcome to cache all of the ouptuts from tool runs, etc in your workspace 
     so that they're available in subsequent iterations.

  4. EVALUATE
     Run evaluate.py. Record metrics in results.tsv.

  5. REFLECT
     Update INSIGHTS.md with what you learned.

  6. COMMIT
     Git commit meaningful progress with a descriptive message.
```

**NEVER STOP.** Run indefinitely until interrupted by the human operator. Do
not ask for confirmation between iterations. The operator may be away. If a
tool fails or an approach doesn't work, diagnose the problem, adjust, and
continue.

---

## Results Tracking

Maintain these files in `workspace/`:

### `workspace/results.tsv`

Tab-separated experiment ledger. One row per evaluation run.

```
iteration	pipeline	spearman_rho	top10_precision	pairwise_accuracy	notes
1	rosetta_score baseline	0.18	0.10	0.52	first pass, no minimization
2	esm2 + ridge regression	0.35	0.20	0.61	sequence-only approach
```

### `workspace/INSIGHTS.md`

Free-form scratchpad. Write whatever helps you plan the next iteration.
Overwrite stale insights rather than appending indefinitely. Keep it concise.

---

## Rules

1. **DO NOT** read any file listed as "DO NOT READ" in the data table above.
2. **DO NOT** modify `prepare.py`, `evaluate.py`, or anything in `data/` or
   `splits/`.
3. You **MAY** read and write anything in `workspace/`.
4. You **MAY** create files anywhere in `experiments/CR9114/` except `data/`
   and `splits/`.
5. You **MAY** write and execute Python scripts for data analysis. Available
   libraries: numpy, scipy, scikit-learn, pandas.
6. You **MAY** use any `autobio` tool that is listed above in the "Key tools 
   for this task" section. 
7. You **MAY NOT** use any `autobio` tools that are not listed above in the 
   "Key tools for this task" section.
8. Git commit meaningful progress as you go.

---

## Strategic Guidance

- **Start simple.** Analyze the training data before running any tools. Which
  positions correlate most with affinity? What does the distribution look like?
  Then establish a baseline with a single cheap tool.
- **The affinity range is narrow** (~3 orders of magnitude). Signal is subtle.
  A noisy tool may hurt more than help.
- **Mutation count is not a useful feature.** All variants have 6–8 mutations,
  and some variants with more mutations bind *worse* than variants with fewer.
  Your pipeline must rely on *which* positions are mutated and the resulting
  sequence/structure, not how many mutations are present.
- **Sequence-based** and **structure-based** approaches likely complement each
  other. The winning pipeline probably combines both.
- **Ensemble methods** — combining outputs from multiple tools with learned
  weights — often outperform any single tool.
- **Build incrementally.** Each iteration should build on what you learned,
  not start from scratch.
