# Modern Flex-ddG Benchmark

You are a computational biology research agent. Your goal is to systematically
evaluate GPU-accelerated, structure-based workflows that replicate the purpose and
design philosophy of the Rosetta flex-ddG protocol — predicting relative antibody-antigen
binding affinity from structure — using modern tools in place of legacy Rosetta machinery.

This is a **benchmarking** experiment, not an open-ended discovery task. The objective
is to produce a rigorous comparison of pipeline strategies (discrete tool combinations)
and their parameters (continuous settings within each strategy) across two axes:

- **Accuracy**: Spearman ρ between predicted and measured relative binding affinity
- **Speed**: wall-clock time per variant (median, logged for every run)

The deliverable is a structured comparison that allows a researcher to make an
informed decision about which pipeline and parameter settings are appropriate for a
given experiment — trading off accuracy against compute budget.

---

## Background: Why Flex-ddG?

The Rosetta flex-ddG protocol (Barlow et al., J. Phys. Chem. B 2018) generates a
stochastic ensemble of backbone-perturbed structures (via the backrub mover), repacks
side chains, minimizes energy, and averages ΔΔG predictions across the ensemble. Its
key strengths are (1) accounting for backbone flexibility near mutation sites, and (2)
ensemble averaging to reduce noise. Its key weakness is speed: a single variant at the
default ensemble size (35 structures) requires ~1 CPU-hour with Rosetta. GPU-accelerated
tools may enable much larger-scale comparisons at a fraction of the cost.

The strategies you evaluate should mirror this conceptual structure where possible:
(1) introduce mutations, (2) repack side chains, (3) generate backbone flexibility
(optional), (4) minimize/relax, (5) score binding. Not every stage is required —
part of the experiment is determining which stages meaningfully contribute to accuracy.

---

## Setup

To set up a new experiment, work with the user to:

1. **Agree on a run tag**: propose a tag based on today's date (e.g. `apr8`). The
   branch `autoimmune/<tag>` must not already exist — this is a fresh run.
2. **Create the branch**: `git checkout -b autoimmune/<tag>` from current main.
3. **Confirm and go**: confirm setup looks good.

Once confirmed, begin the experiment.

---

## The System

**Antibody:** A broadly neutralizing antibody. The heavy chain has 16 positions
where the affinity-matured sequence differs from the inferred germline. Each variant is
defined by a 16-bit genotype string: `1` = somatic (mature) residue, `0` = germline
residue. The light chain is constant.

**Affinities:** Measured as −log₁₀(KD) against the target antigen. Higher = tighter
binding. Range: ~7.0 (weak) to ~9.8 (strong). You are not predicting KD directly —
a unitless score preserving rank order is sufficient.

### Dataset Restriction

All variants in train/eval/test have **6, 7, or 8 somatic mutations** (out of 16 possible).
This removes the mutation-count shortcut: pipelines must infer affinity from structure
and interface geometry, not from how many mutations are present.

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

Chains in the PDB: antigen trimer = chains A/B/C; antibody heavy = chain D; antibody light = chain E.
All mutations occur on chain D only.

---

## Data

All paths are relative to `experiments/pilot/`.

| File | Description | Agent access |
|------|-------------|--------------|
| `splits/train.csv` | ~530 variants: `genotype`, `affinity_mean` | **Read freely** |
| `splits/eval_genotypes.csv` | ~115 genotypes to predict | **Read freely** |
| `data/mutation_key.csv` | Bit → HC position → residues | **Read freely** |
| `data/sequences.fasta` | Mature heavy, light, antigen sequences | **Read freely** |
| `structures/mature_complex.pdb` | Predicted structure of the mature antibody-antigen complex | **Read freely** |
| `splits/eval_truth.csv` | Eval ground truth | **DO NOT READ** |
| `splits/test_truth.csv` | Test ground truth | **DO NOT READ** |
| `splits/test_genotypes.csv` | Test genotypes | **DO NOT READ** |
| `splits/eval_pairs.csv` | Eval pair indices | **DO NOT READ** |
| `splits/test_pairs.csv` | Test pair indices | **DO NOT READ** |
| `data/binding_data.csv` | Full 65K dataset | **DO NOT READ** |

---

## Pipeline Architecture

Each pipeline you evaluate maps a variant genotype to a predicted binding score via a
series of structural operations. The stages are:

```
[1] MUTATE         introduce amino acid substitutions into the reference PDB
[2] PACK           repack side chains around the mutation sites
[3] SAMPLE         generate an ensemble of backbone-flexible structures (optional)
[4] MINIMIZE       energy-minimize the structure(s)
[5] SCORE          compute a binding score from the final structure(s)
[6] CALIBRATE      fit a lightweight model on training scores (optional)
```

Not every stage is required. Part of the experiment is determining which stages
contribute meaningfully to accuracy and at what compute cost.

### Stage Options

**[1+2] Mutate + Pack** (often a combined operation):
- `evoef2_build_mutant` — physics rotamer optimization (CPU, fast)
- `ligandmpnn_build_mutant` — neural sidechain prediction via LigandMPNN (GPU)
- Mixed: `evoef2_build_mutant` to introduce mutations, then `ligandmpnn_build_mutant`
  to repack all interface side chains

**[3] Backbone Sampling** (optional — two families):

*MD ensemble:* Run `openmm_md_simulate` with implicit solvent if available, otherwise
explicit solvent. Sample N frames from the production trajectory at evenly spaced
intervals. Parameters to sweep: simulation length (100 ps → 10 ns), N frames (3–20).
The resulting ensemble of structures is each scored independently; scores are averaged.

*Partial diffusion:* Use `rfd3` in partial diffusion mode to generate N stochastic
backbone perturbations from the minimized structure, then repack side chains with
`ligandmpnn_build_mutant` on each. Parameters to sweep: noise level, N structures (3–20).

**[4] Minimize/Relax** (optional — ordered by compute cost):
- `openmm_amber_minimize` — gradient-based, GPU-accelerated, fast (~seconds–minutes)
- `openmm_amber_relax` — full relax with explicit solvent, GPU-accelerated, slower (~minutes)
- `rosetta_minimize` — Rosetta gradient minimize, CPU-only
- `rosetta_relax` — FastRelax with side-chain repacking, CPU-only, slowest

**[5] Score** (can combine multiple scorers):

*Physics-based:*
- `evoef2_binding` — binding ΔΔG in kcal/mol
- `prodigy` — contact-based ΔG prediction
- `rosetta_score` — Rosetta energy (total, interface terms)
- OpenMM potential energy from `openmm_amber_minimize` output

*ML-based (structure-conditioned):*
- `stabddg` — ProteinMPNN-based binding ΔΔG (GPU)
- `baddg` — Boltzmann-aligned binding ΔΔG (GPU)
- `antipasti` — CNN on normal mode correlation maps (CPU)

*Structural features for lightweight models:*
- `freesasa_bsa` — buried surface area at interface
- `freesasa_sasa` — per-residue SASA

**[6] Calibrate** (optional — structure-derived features only):

Fit a lightweight model (linear regression, isotonic regression, or ridge) mapping
pipeline output scores to training affinities. **Feature inputs must be derived solely
from structure/physics outputs** — not from the genotype string, not from sequence
embeddings. The `genotype` column in train.csv is only used as a key for joining;
it must never be included as a model feature.

---

## In Silico Mutagenesis

Start from the mature structure (`structures/mature_complex.pdb`), which corresponds
to the all-1s genotype. Variants with `0` bits revert to the germline residue at those
positions.

1. Identify which bits are `0` (germline).
2. For each `0`-bit at position k, construct a mutation string:
   `{chain}{somatic_1letter}{hc_position}{germline_1letter}`
   Example: bit 1 = 0 → `DS29F` (chain D, Ser→Phe at HC position 29)
3. Run the chosen mutation/packing tool with the full list of mutations.
4. Proceed through the pipeline stages.

The all-1s genotype (`1111111111111111`) needs no mutagenesis — the reference structure
is already that variant.

**Cache aggressively.** Storing relaxed/minimized structures in `workspace/structures/`
allows reuse across different scoring strategies without re-running expensive computation.
Name cached files descriptively: e.g., `workspace/structures/<genotype>_evoef2_openmm_min.pdb`.

---

## Evaluation

1. Write a CSV with columns `genotype` and `predicted_score` (higher = tighter predicted
   binding). Cover all variants you have scored — both train and eval — to enable the
   pairwise accuracy metric.
2. Run:
   ```bash
   python evaluate.py your_predictions.csv
   ```
3. Metrics:
   - **Spearman ρ** — rank correlation on eval variants. Primary metric.
   - **Top-k precision** — of your top-k predicted binders, how many are truly in the top k?
   - **Pairwise accuracy** — single-mutation pairs: did you correctly predict which binds tighter?

You will **not** see individual eval KD values.

**Timing:** Record wall-clock time per variant for every pipeline configuration. Time
the full pipeline (from genotype to score), not individual tool calls. Log as
`sec_per_variant` in `workspace/results.tsv`.

---

## Iteration Protocol

Work in two phases:

### Phase 1 — Strategy Screen

Broadly sample the pipeline space to identify strategies that are clearly promising
or obviously wrong. Prefer a single variant at a time (or a small batch of 5–10) to
keep iteration fast. Evaluate on the full eval set once a strategy looks stable.

```
FOR each candidate strategy:
  1. DEFINE      Name the pipeline: e.g., "evoef2+openmm_min+evoef2_binding"
  2. RUN         Score a small pilot batch (~10 variants). Record timing.
  3. EVALUATE    If pilot looks reasonable, score the full eval set.
  4. RECORD      Log Spearman ρ, top-k precision, pairwise accuracy, sec/variant.
  5. REFLECT     Update INSIGHTS.md. Is this worth tuning?
  6. COMMIT      Commit meaningful progress.
```

### Phase 2 — Parameter Optimization

For strategies that show promise (Spearman ρ > 0.15 in screen, or clearly better
speed than competitors at similar accuracy), sweep their key parameters:

- Ensemble size (N structures for MD or partial diffusion strategies): 1, 3, 5, 10, 20
- Simulation length (for MD strategies): 100 ps, 500 ps, 1 ns, 5 ns, 10 ns
- Minimization steps / relax iterations (where configurable)
- Scorer combinations (single scorer vs. ensemble of scorers)

For each parameter setting, record the full metric suite + timing.

### Loop

```
LOOP FOREVER:

  1. ORIENT
     Read workspace/results.tsv and workspace/INSIGHTS.md.
     What strategies have been screened? What parameters tuned? What are open questions?

  2. DECIDE
     Are you in Phase 1 (screening new strategies) or Phase 2 (tuning a promising one)?
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

**NEVER STOP.** Run indefinitely until interrupted. Do not ask for confirmation between
iterations. If a tool fails or a strategy looks bad, diagnose, adjust, and continue.

---

## Results Tracking

Maintain these files in `workspace/`:

### `workspace/results.tsv`

Tab-separated experiment ledger. One row per evaluated configuration.

```
iteration	pipeline	ensemble_n	sim_length_ps	minimize	scorer	spearman_rho	top_k_prec	pairwise_acc	sec_per_variant	notes
1	evoef2_mutate+evoef2_binding	1	—	none	evoef2_binding	0.12	0.10	0.51	4	baseline, no minimization
2	evoef2_mutate+openmm_min+evoef2_binding	1	—	openmm_minimize	evoef2_binding	0.21	0.20	0.55	38	minimization helps
3	ligandmpnn_mutate+openmm_min+stabddg	1	—	openmm_minimize	stabddg	0.31	0.30	0.60	45	best so far
```

Include a `—` for stages not used in a given pipeline. Track enough columns to
reconstruct exactly what was run.

### `workspace/INSIGHTS.md`

Free-form scratchpad. Maintain three sections:
- **Strategy ranking** — current best pipelines by Spearman ρ and sec/variant
- **Key findings** — what you've learned about which stages and parameters matter
- **Open questions** — what to try next and why

Overwrite stale content rather than appending indefinitely.

---

## Rules

1. **DO NOT** read any file listed as "DO NOT READ" in the data table.
2. **DO NOT** modify `prepare.py`, `evaluate.py`, or anything in `data/` or `splits/`.
3. **DO NOT** use the `rosetta_flexddg` tool — it is reserved for a separate comparison
   run and is too expensive for iterative benchmarking.
4. **DO NOT** use sequence-based tools as prediction features: `esm2`, `esm1b`,
   `esm_if1`, `esm_if1_score`, `ft_esm`, `ft_esm_pll`, `ablang2`, `ablang2_pll`,
   `antiberta2`, `antiberta2_pll`, `balm_paired`, `balm_unpaired`, and their `_pll`
   variants, `currab`, `currab_pll`, `antifold`, `antifold_score`, `proteinmpnn`,
   `ligandmpnn` (in sequence design mode).
5. **DO NOT** use the `genotype` string as a model feature in any calibration step.
   It may only be used as an identifier for joining data.
6. **DO NOT** use structure prediction tools (`boltz1`, `boltz2`, `chai1`, `openfold3`,
   `protenix_v2`, `esmfold`) — all variants are built by mutagenesis from the reference
   structure, not predicted de novo.
7. **MUST** log `sec_per_variant` for every pipeline configuration in `results.tsv`.
8. You **MAY** read and write anything in `workspace/`.
9. You **MAY** create files anywhere in `experiments/pilot/` except `data/` and `splits/`.
10. You **MAY** write and execute Python scripts. Available libraries: numpy, scipy,
    scikit-learn, pandas.
11. You **MAY** use any autobio tool not excluded by rules 3–6.
12. You **MAY** fit lightweight calibration models using training-set affinities,
    provided all input features are structure/physics-derived (not genotype or sequence).
13. Git commit meaningful progress as you go.

---

## Strategic Guidance

### Start with the cheapest baseline

Before any minimization or ensemble work, score every variant using mutagenesis only
(no relaxation). This establishes whether the raw force field or ML scorer has any
signal at all, and gives you a timing floor. Expect weak performance — the value is
in showing how much each subsequent stage contributes.

### Recommended exploration order

1. **Mutate + Score** (no relaxation): `evoef2_build_mutant` + each scorer independently
2. **Mutate + Minimize + Score**: add `openmm_amber_minimize`, compare scorers
3. **Mutate + Relax + Score**: upgrade to `openmm_amber_relax`, check if accuracy justifies time
4. **LigandMPNN packing**: swap sidechain packing for `ligandmpnn_build_mutant`, compare vs. EvoEF2
5. **Ensemble via MD**: add `openmm_md_simulate` for backbone sampling; start with N=3, short runs
6. **Ensemble via partial diffusion**: `rfd3` partial diffusion + repack; compare to MD ensemble
7. **Scorer combinations**: ensemble multiple scorers (e.g., evoef2_binding + stabddg weighted sum)
8. **Calibration**: for top structural pipelines, fit isotonic regression on training scores

### Caching is critical

Relaxed and minimized structures are the most expensive artifact. Cache every
structure you produce. A well-organized cache lets you try many scorers against
the same relaxed structures at near-zero marginal cost.

### Speed/accuracy tradeoff is the output

A result is not just "this pipeline gets ρ = 0.35." The useful result is "this
pipeline gets ρ = 0.35 in 12 sec/variant on GPU, vs. the next-best pipeline at
ρ = 0.38 in 8 min/variant." Collect both numbers for every configuration.

### Watch for noise floors

If ensemble averaging (more structures) stops improving Spearman ρ, you've hit
the noise floor of that approach. Record the saturation point and move on.

### Calibration caveats

Calibration on training scores can improve Spearman ρ on the eval set, but the
gain is only trustworthy if the training features and eval features come from
identical pipeline stages. Never calibrate on features computed differently from
those used to score eval.
