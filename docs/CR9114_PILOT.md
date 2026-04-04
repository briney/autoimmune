# AutoImmune: Antibody Affinity Prediction

**Design Specification for Autonomous Pipeline Discovery**

CR9114 × HA1 Pilot Experiment

Version 1.0 | April 2026

---

## 1. Overview

### 1.1 Motivation

Predicting how mutations in an antibody affect its binding affinity for a target antigen is a central challenge in computational immunology. Existing tools span a wide range of approaches, from physics-based energy functions (Rosetta, OpenMM/Amber) to machine-learned scoring models (StaB-ddG) and sequence embeddings (ESM-2). Each tool has different strengths, weaknesses, computational costs, and parameter sensitivities. No single tool is universally best, and the optimal workflow — including the choice of tools, their ordering, parameter settings, and how to combine their outputs — is an open question that has traditionally been explored through manual, labor-intensive benchmarking.

This project applies an autonomous research paradigm, inspired by Karpathy's autoresearch framework, to systematically discover an optimized computational pipeline for predicting antibody-antigen binding affinity. Rather than training neural networks in short iterative runs, the agent composes and evaluates pipelines of structural biology tools from the autobio toolkit, iteratively refining its approach based on quantitative feedback.

### 1.2 Scope

This document specifies a pilot experiment restricted to a single antibody-antigen system: the broadly neutralizing anti-influenza antibody CR9114 binding to influenza HA1 (H1N1 hemagglutinin). The dataset comprises a combinatorial library of CR9114 heavy-chain variants with experimentally measured binding affinities from Phillips et al. 2021 (eLife 10:e71393). The design is intended to be extensible to additional antibodies (CR6261), antigens (H3, Flu-B), and eventually to novel antibody-antigen pairs.

### 1.3 Success Criteria

The primary success metric is the Spearman rank correlation (ρ) between the agent's predicted binding scores and the experimentally measured binding affinities on a held-out final test set. The agent is not required to predict actual KD values — a unitless binding score that preserves the rank ordering of affinities is the target. Secondary metrics include top-k precision (ability to identify the strongest binders from a pool) and pairwise accuracy on matched single-mutation variant pairs (ability to determine which member of a pair binds more tightly).

---

## 2. Dataset

### 2.1 Source Data

The Phillips et al. 2021 dataset provides a complete combinatorial mutagenesis library for CR9114. The antibody heavy chain has 16 positions where the mature (somatic) residue differs from the inferred germline residue. Each of the 2^16 = 65,536 possible combinations of germline/somatic residues at these 16 positions was constructed and tested for binding against three influenza HA variants. The light chain is held constant across all variants.

The 16 mutable positions in the CR9114 heavy chain are:

| Bit | HC Position | Germline | Somatic | H1 Impact Rank |
|-----|-------------|----------|---------|----------------|
| 1   | 29 (CDR1)   | Phe (F)  | Ser (S) | 4th            |
| 2   | 30 (CDR1)   | Ser (S)  | Asn (N) | Low            |
| 3   | 31 (CDR1)   | Ser (S)  | Asn (N) | Low            |
| 4   | 52 (CDR2)   | Ile (I)  | Ser (S) | 2nd            |
| 5   | 57 (CDR2)   | Thr (T)  | Ser (S) | 6th            |
| 6   | 58 (FR3)    | Ala (A)  | Thr (T) | 5th            |
| 7   | 59 (FR3)    | Asn (N)  | Ala (A) | 7th            |
| 8   | 71 (FR3)    | Thr (T)  | Ser (S) | Low            |
| 9   | 74 (FR3)    | Lys (K)  | Ile (I) | 3rd            |
| 10  | 75 (FR3)    | Ser (S)  | Phe (F) | 1st            |
| 11  | 76 (FR3)    | Thr (T)  | Ser (S) | Low            |
| 12  | 77 (FR3)    | Ser (S)  | Asn (N) | Low            |
| 13  | 84 (FR3)    | Ser (S)  | Asn (N) | Low            |
| 14  | 87 (FR3)    | Arg (R)  | Thr (T) | Low            |
| 15  | 95 (CDR3)   | Tyr (Y)  | Phe (F) | Low            |
| 16  | 106 (FR4)   | Tyr (Y)  | Ser (S) | Low            |

### 2.2 H1 Binding Landscape Characteristics

Restricting to H1 binding only, the dataset contains 65,094 variants with measured affinities (442 of the 65,536 genotypes lack H1 data). Key characteristics of the H1 binding landscape:

- **Dynamic range:** KD values span from 1.46 × 10⁻¹⁰ M (146 pM, the tightest binder) to 1.00 × 10⁻⁷ M (100 nM, the weakest), a ~685-fold range. This is roughly 2.8 orders of magnitude in log-space.
- **Distribution:** The landscape is dominated by functional binders. 68.8% of variants are strong binders (KD < 1 nM), 28.6% are moderate (1–100 nM), and only 2.6% are weak (>100 nM). There is no hard non-binding ceiling; 100 nM appears to be the assay floor.
- **Mutation count vs. affinity:** Mean affinity increases monotonically with mutation count (more somatic mutations = tighter binding), from –log₁₀(KD) = 8.42 for the germline (0 mutations) to 9.59 for the fully mature antibody (16 mutations). However, there is substantial variance at every mutation count.
- **High-leverage positions:** Position 10 (S75F) is by far the most impactful single mutation for H1 binding, appearing in 12,925 single-mutation pairs with >10-fold affinity change. Positions 4 (I52S) and 9 (K74I) are the next most impactful. Most positions (2, 3, 8, 11–16) have individually modest effects.
- **Epistasis:** The monotonic mutation-count trend combined with large single-position effects suggests significant epistasis — the effect of one mutation depends on which other mutations are present. This is exactly the phenomenon that a structure-aware pipeline should be able to capture.

### 2.3 Downsampling Strategy

The full 65K-variant landscape is too large for efficient iteration. Each agent iteration must run its pipeline on the entire eval set and selectively on training variants, so dataset size directly impacts per-iteration compute cost. We downsample to a target of approximately 3,000 training variants, 500 eval variants, and 500 final test variants using a stratified approach that preserves the informative structure of the landscape.

#### 2.3.1 Stratification Procedure

1. **Affinity binning.** Divide variants into four strata by –log₁₀(KD): strong binders (≥ 9.5), good binders (9.0–9.5), moderate binders (8.5–9.0), and weak binders (< 8.5). Because the landscape is skewed toward strong binders, sample roughly equal numbers from each stratum rather than proportionally.
2. **Mutation-count coverage.** Within each affinity stratum, ensure representation across the full range of mutation counts (0–16). Variants with very few (0–3) or very many (13–16) mutations are rare and should be preferentially retained.
3. **Positional diversity.** Verify that every one of the 16 mutation positions appears in both its somatic and germline state across a range of genetic backgrounds in every dataset partition. This ensures the agent can learn position-specific effects.
4. **Informative pair inclusion.** Identify all single-mutation neighbor pairs where the two members differ by >10-fold in KD. Include ~100–150 such pairs, distributed so that one member is in train and the other in eval or final test. These pairs are the hardest and most informative test cases.

#### 2.3.2 Three-Way Split

| Partition   | Size   | Agent Access                                                                 | Purpose                                      |
|-------------|--------|-----------------------------------------------------------------------------|----------------------------------------------|
| Train       | ~3,000 | Full: sequences, genotypes, KD values, queryable via API                    | Pipeline development and calibration         |
| Eval        | ~500   | Sequences and genotypes only; aggregate metrics returned after each iteration | Iterative feedback signal for pipeline refinement |
| Final Test  | ~500   | None until experiment concludes                                              | True held-out generalization measure         |

The eval set functions as a validation set: the agent uses its aggregate metrics to guide pipeline development, making it a contaminated signal to some degree. The gap between eval ρ and final test ρ measures how much eval-overfitting occurred. The final test set is never touched until after the last iteration.

---

## 3. Pre-Populated Structures

Structure prediction is computationally expensive and would consume a disproportionate share of the agent's compute budget if performed during the iterative loop. To eliminate this bottleneck, we pre-compute predicted structures for the CR9114–HA1 complex and provide them as given inputs in the agent's workspace.

### 3.1 Structures to Generate

- **Mature CR9114 + HA1:** The fully somatically mutated heavy chain (genotype 1111111111111111) paired with the constant light chain, in complex with HA1. This is the primary reference structure for the agent. No experimentally determined structure exists for this exact complex.
- **Germline CR9114 + HA1:** The fully germline-reverted heavy chain (genotype 0000000000000000) paired with the constant light chain, in complex with HA1. This serves as a second structural reference point, particularly useful for understanding how accumulated mutations reshape the binding interface.

### 3.2 Prediction Method

Structures will be predicted using Boltz-2 via the autobio toolkit. Boltz-2 handles antibody-antigen complexes and produces confidence scores (pLDDT, PAE) that can be used to assess prediction quality. If Boltz-2 predictions show poor interface confidence, Chai-1 will be used as a fallback. Predictions will be validated by visual inspection of the paratope-epitope interface and comparison with published crystal structures of related CR9114 complexes (PDB: 4FQI for CR9114-H3, 4FQY for CR9114-H1 stem).

### 3.3 Agent Usage

The agent receives these structures as PDB files in its workspace. To evaluate a mutant variant, the agent's pipeline will: (1) introduce the relevant mutations into the closest reference structure via in silico mutagenesis, (2) minimize or relax the mutant structure to resolve steric clashes, and (3) score the binding energy of the resulting complex. The choice of reference structure, mutagenesis approach, minimization protocol, and scoring function are all degrees of freedom that the agent optimizes.

---

## 4. Autobio Toolkit Integration

### 4.1 Available Tools

The agent has access to a curated subset of the autobio toolkit. Tools are selected based on relevance to the affinity prediction task and compatibility with the per-iteration credit budget. All tools are invoked via the autobio CLI and run in Docker containers with standardized input/output contracts.

| Tool               | Category      | Credit Rate      | Description                                                                                           |
|--------------------|---------------|------------------|-------------------------------------------------------------------------------------------------------|
| rosetta-score      | Scoring       | 1 / structure    | Evaluate structure energy using Rosetta ref2015 or other score functions. Returns total score, per-residue breakdown, and interface energy. |
| evoef2             | Scoring       | 1 / structure    | Fast physics-based scoring. Complementary to Rosetta with different energy decomposition.              |
| stab-ddg           | Scoring       | 3 / mutation     | ML-based ΔΔG prediction using ProteinMPNN architecture. Predicts binding energy change from mutation. Requires GPU but fast inference. |
| esm-2              | Embeddings    | 2 / sequence     | 650M-parameter protein language model embeddings. Useful for building lightweight ML regressors on top of structural scores. |
| esm-1b             | Embeddings    | 1 / sequence     | Smaller ESM model. Cheaper embeddings, potentially sufficient for this task.                           |
| rosetta-minimize   | Minimization  | 5 / structure    | Gradient-based energy minimization in Rosetta. Fast, resolves minor clashes.                          |
| omm-amber-minimize | Minimization  | 8 / structure    | OpenMM energy minimization with Amber force field. More thorough than Rosetta minimize.               |
| omm-amber-relax    | Relaxation    | 15 / structure   | Short OpenMM relaxation protocol. Allows limited backbone flexibility.                                 |

### 4.2 Excluded Tools

The following autobio tools are excluded from the agent's available set. They are either too computationally expensive for iterative use or not relevant to the affinity prediction task.

| Tool                          | Reason for Exclusion                                                                    |
|-------------------------------|----------------------------------------------------------------------------------------|
| flex-ddg                      | Rosetta ensemble backrub sampling. Hours per mutation — incompatible with credit budget. |
| rosetta-relax                 | Full Rosetta FastRelax protocol. Significantly slower than rosetta-minimize.             |
| omm-md-simulate               | Full molecular dynamics trajectories. Expensive and unnecessary when minimize/relax suffice. |
| boltz-1, boltz-2              | Structure prediction. Pre-populated structures eliminate the need for runtime prediction. |
| chai-1, esmfold, openfold3    | Structure prediction. Same rationale.                                                    |
| rfdiffusion3, proteina-*      | Structure design tools. Not relevant to scoring existing variants.                       |
| proteinmpnn, ligandmpnn       | Inverse folding. Not directly relevant to affinity scoring.                              |

### 4.3 Tool Invocation

All tools are invoked via the autobio CLI with JSON-formatted output for programmatic parsing:

    autobio run <tool-name> --config <config.json> --format json

The agent can query tool specifications and usage instructions at any time:

    autobio info <tool-name> --format json

Each tool execution produces a standardized output in the workspace directory, including raw tool output, standardized results, logs, and execution metadata. The harness captures these outputs and populates the results cache.

---

## 5. Credit System

### 5.1 Design Rationale

The credit system serves as a forcing function for computational efficiency. Without it, the agent could spend unlimited time running expensive tools on all variants, making iterations impractically long. Credits force the agent to make strategic choices: run a cheap tool broadly, an expensive tool selectively, or find creative combinations that balance coverage and depth.

### 5.2 Budget

- **Per-iteration budget: 4,000 credits.** This is enough to run a single cheap scoring tool (cost: 1 credit/structure) across the entire ~3,000-variant training set plus the ~500-variant eval set, with a modest surplus for experimentation. Alternatively, it funds rosetta-minimize on ~500 structures, or omm-amber-minimize on ~350. The budget is deliberately tight to drive creative pipeline design.
- **Total iteration cap: 15 iterations.** This gives the agent a maximum of 60,000 credits over the full experiment. Given caching, the effective budget grows over time as prior results accumulate.
- **Eval set cost:** Running the agent's pipeline on the 500-variant eval set costs credits like any other tool use. The agent must budget for this within each iteration's allocation.

### 5.3 Cost Accounting

Credits are computed as: credit_rate × n_structures (or n_sequences / n_mutations as appropriate for the tool). Multi-step pipelines accumulate costs at each stage. For example, a pipeline of rosetta-minimize (5/structure) followed by rosetta-score (1/structure) applied to 300 variants costs 300 × (5 + 1) = 1,800 credits.

Cached results are free. If the agent scored 200 variants with rosetta-score in iteration 1, those 200 scores are available at zero cost in all subsequent iterations. Only novel tool invocations on previously-unscored (tool, parameters, variant) combinations cost credits.

The harness tracks credit expenditure in real time and halts tool execution if the budget would be exceeded. The agent can query its remaining budget at any point.

---

## 6. Queryable Training Dataset API

### 6.1 Design Principles

The training dataset (~3,000 variants) is too large to hold in the agent's context window, so it is exposed through a query API. The API is designed to guide the agent toward productive exploration by providing both raw data retrieval and analytical summaries. All queries are free (no credit cost). A hard cap of 200 records per query prevents context-window overflow and forces the agent to think about which subset is most informative.

### 6.2 Query Operations

#### 6.2.1 Summary Queries

These return aggregate statistics about the training set without retrieving individual records.

    query summary

Returns total variant count, KD distribution statistics (min, max, median, quartiles, mean), mutation count distribution, and the number of cached tool results available.

    query affinity_distribution --bins 20

Returns a histogram of –log₁₀(KD) values, showing how training variants are distributed across the affinity range.

    query mutation_frequency

For each of the 16 mutation positions, returns the fraction of training variants carrying the somatic vs. germline residue, along with the mean KD for variants with/without each mutation. Provides a quick first-pass view of per-position effects.

#### 6.2.2 Filtered Data Retrieval

The primary query operation retrieves specific subsets of training data matching filter criteria.

    query variants \
      --antigen H1 \
      --n-mutations 3:5 \
      --mutations-include 10 \
      --mutations-exclude 4,9 \
      --kd-range 1e-10:1e-8 \
      --sort-by kd \
      --limit 100

Supported filter parameters:

| Parameter          | Type             | Description                                                                  |
|--------------------|------------------|-----------------------------------------------------------------------------|
| --n-mutations      | int or range     | Total somatic mutations present. Single value or min:max range.              |
| --mutations-include| list of ints     | Bit positions that MUST be somatic (mutated) in returned variants.           |
| --mutations-exclude| list of ints     | Bit positions that MUST be germline in returned variants.                    |
| --kd-range         | float:float      | Min:max KD in Molar. Filters to variants within this affinity range.        |
| --genotypes        | list of strings  | Specific 16-bit genotype strings to retrieve.                                |
| --hamming-from     | string           | Reference genotype for Hamming distance filtering.                           |
| --hamming-range    | int:int          | Min:max Hamming distance from the --hamming-from reference.                  |
| --sort-by          | string           | Sort order: kd (ascending), kd_desc, n_mutations, random.                   |
| --limit            | int              | Maximum records to return. Hard cap at 200.                                  |

Each returned record includes: genotype (16-bit string), KD value, –log₁₀(KD), mutation count, and a list of which positions are mutated.

#### 6.2.3 Comparative Queries

Specialized queries designed for probing mutational effects and epistasis.

    query neighbors --genotype 1111111111111111 --distance 1

Returns all training variants within the specified Hamming distance of a reference genotype, with their KD values. Essential for discovering which single mutations cause large affinity shifts.

    query mutation_impact --position 10

For a given mutation position, returns paired statistics: mean and median KD for variants carrying the somatic vs. germline residue at that position, stratified by total mutation count. Estimates the marginal effect of a single position.

    query epistasis --positions 4,10

Returns the 2×2 table of mean affinities for the four combinations of germline/somatic at two specified positions (both germline, position A only, position B only, both somatic). Probes pairwise epistatic interactions.

#### 6.2.4 Sequence Retrieval

    query sequences --genotypes 1111111111111111,0000000000000000 --format fasta

Returns the full amino acid sequences (heavy chain, light chain, antigen) for specified genotypes. The agent needs these when preparing inputs for structural tools. Sequences are long (121 + 111 + 547 = 779 residues per complex), so this query should be used selectively.

#### 6.2.5 Eval and Test Set Queries

    query eval_genotypes

Returns the list of genotypes and mutation counts in the eval set (no KD values). The agent needs this to know what it must predict.

    query test_exists

Confirms the final test set exists and returns its size. No genotypes or sequences are revealed until the experiment concludes.

---

## 7. Results Cache

### 7.1 Architecture

The results cache is a two-tier system: a full results store on disk (Tier 1) and a compressed context-window summary (Tier 2). Tier 1 stores every tool output ever computed; Tier 2 is what the agent actually sees at the start of each iteration.

### 7.2 Tier 1: Disk Store

Every tool execution produces a cache record stored as a JSON file on disk:

```json
{
  "run_id": "iter03_run07",
  "iteration": 3,
  "tool": "rosetta-score",
  "parameters": {"score_function": "ref2015"},
  "genotypes": ["1111111111111111", "1111111111111110", ...],
  "results": {
    "1111111111111111": {"total_score": -342.5, "interface_score": -18.3},
    ...
  },
  "credits_spent": 48,
  "wall_time_seconds": 120
}
```

The agent can query cached results at any time at zero credit cost:

    cache lookup --tool rosetta-score --params ref2015 --genotypes 1111111111111111

Returns cached scores if available, or indicates which genotypes have no cached results for the specified tool/parameter combination. The agent uses this to determine which variants need new tool runs versus which can be retrieved from cache.

### 7.3 Tier 2: Context-Window Summary

At the start of each iteration, the agent receives a compressed summary of all cached results and performance history. This summary is the agent's persistent memory across iterations. It consists of three sections, targeting a total length of approximately 1,500–2,000 tokens.

#### 7.3.1 Section A: Run Log

A compact table of all tool runs completed across all iterations:

```
COMPLETED RUNS (23 total | 47,200 credits spent of 60,000 lifetime)
════════════════════════════════════════════════════════════════
Iter  Tool               N_var  Params             Credits  Cache_key
1     rosetta-score       500    ref2015            500      rs_ref_500
1     stab-ddg            200    default            600      sd_def_200
2     omm-amber-min+rs    300    steps=500,ref2015  1800     omm500_rs_300
3     esm-2               1000   layer=-1           2000     esm2_1000
```

As this table grows, older entries may be collapsed to show only the tool, variant count, and cache key, preserving the detailed entries for the most recent 3–5 iterations.

#### 7.3.2 Section B: Performance Tracker

The cumulative record of pipeline performance across iterations:

```
PIPELINE PERFORMANCE HISTORY
════════════════════════════════════════════════════════════════
Iter  Pipeline                            Train ρ   Eval ρ   Top-50  Pair Acc
1     rosetta-score(ref2015)              0.23      0.19     24%     58%
2     omm-min(500)+rosetta-score          0.31      0.27     32%     63%
3     0.6*stabddg + 0.4*omm_min+rs       0.41      0.38     48%     71%
4     esm2-ridge + stabddg (ensemble)     0.52      0.44     56%     74%
```

Each row records the agent's best pipeline for that iteration and four metrics: Spearman ρ on train and eval, top-50 precision (*of the predicted top 50 binders, how many are actually in the true top 50*), and pairwise accuracy on single-mutation neighbor pairs in the eval set.

#### 7.3.3 Section C: Agent Insights

A free-form scratchpad maintained by the agent. At the end of each iteration, the agent is prompted to write or update its accumulated knowledge about the problem. This section has no imposed structure — the agent writes whatever it finds most useful for guiding its next iteration. Examples of what might appear:

```
AGENT INSIGHTS (updated after iteration 4)
════════════════════════════════════════════════════════════════
Position 10 (S75F) dominates H1 affinity. Rosetta interface_score
captures this well but underestimates the effect when positions 4
and 9 are both germline. StaB-ddG is better at these epistatic
cases. Current best strategy: use StaB-ddG as the primary scorer,
with Rosetta interface_score as a tiebreaker for close calls.

OpenMM minimization with 500 steps improves Rosetta scoring by
~0.05 rho over unminimized structures, but 2000 steps gives no
further improvement. 500 steps is the sweet spot for cost/benefit.

ESM-2 embeddings + ridge regression achieves rho ~0.45 on train
without any structural information. This provides a strong baseline
that structure-based methods should beat. Combining ESM-2 features
with structural scores may be the winning approach.
```

Section C has a soft limit of ~500 tokens. The agent is encouraged to be concise and to overwrite stale insights rather than appending indefinitely.

---

## 8. Iteration Loop

### 8.1 Iteration Flow

Each iteration follows a fixed six-phase sequence:

1. **Orient.** The agent receives: its iteration number, remaining per-iteration credits (4,000), remaining lifetime credits, the full Tier 2 context summary (Sections A, B, and C from the previous iteration), and the list of available tools with their credit costs.
2. **Hypothesize.** The agent articulates what it wants to try this iteration and why, based on prior results. On the first iteration, this is an initial strategy; on subsequent iterations, it should reference specific insights from the performance history and its own notes.
3. **Execute.** The agent queries the training dataset, checks the cache for existing results, and runs tools on new variants. It assembles a scoring pipeline — a defined procedure that takes a genotype and produces a predicted binding score. Credits are deducted in real time.
4. **Evaluate on train.** The agent computes its pipeline's predictions for all training variants (using cached results where available) and calculates Spearman ρ, top-k precision, and pairwise accuracy against the known training KD values. The agent has full visibility into individual training errors.
5. **Evaluate on eval.** The agent applies its pipeline to the eval set. The harness computes Spearman ρ, top-50 precision, and pairwise accuracy, and returns only these aggregate metrics — the agent never sees individual eval KD values. This is the signal the agent uses to judge generalization.
6. **Reflect.** The agent updates Section C (insights scratchpad) and the harness updates Sections A and B. The iteration advances.

### 8.2 Eval Metrics

Three metrics are reported for both train and eval evaluations:

| Metric            | Definition                                                                                 | What It Tests                                                                                          |
|-------------------|-------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------|
| Spearman ρ        | Rank correlation between predicted scores and true –log₁₀(KD)                             | Overall ranking accuracy across the full affinity range. Primary objective.                             |
| Top-50 Precision  | Fraction of the predicted 50 strongest binders that are in the true top 50                 | Ability to identify high-affinity variants from a pool. The practical use case for antibody optimization. |
| Pairwise Accuracy | Fraction of single-mutation neighbor pairs where the predicted higher-scorer is the true higher-affinity binder | Fine-grained local resolution. Can the pipeline distinguish between closely related variants?           |

The eval set contains approximately 50–75 single-mutation neighbor pairs (identified during dataset construction) for the pairwise accuracy calculation. Eval-set metrics are broken down by affinity quartile when possible to diagnose whether the pipeline performs better on strong vs. weak binders.

### 8.3 Termination

The experiment terminates when any of the following conditions are met: the agent exhausts its 15-iteration cap, the agent exhausts its lifetime credit budget (60,000 credits), or the agent declares convergence (eval ρ has not improved by more than 0.01 for three consecutive iterations). Upon termination, the best pipeline (by eval Spearman ρ) is applied to the final test set, and the results are reported.

---

## 9. Agent Capabilities and Constraints

### 9.1 What the Agent Can Do

- Invoke any available autobio tool on any subset of variants, subject to credit limits.
- Query the training dataset using any combination of filter parameters.
- Look up cached results from any prior iteration at zero cost.
- Write and execute custom Python scripts for data analysis, including training lightweight ML models (e.g., ridge regression, gradient boosting) on cached features.
- Define ensemble scoring pipelines that combine outputs from multiple tools with learned weights.
- Choose which reference structure (mature or germline) to use for mutagenesis.
- Specify tool parameters (score functions, minimization steps, force field variants).
- Query autobio tool documentation via `autobio info` for parameter guidance.

### 9.2 What the Agent Cannot Do

- Access individual KD values in the eval or final test sets.
- Modify the train/eval/test split.
- Use excluded tools (flex-ddg, rosetta-relax, MD simulation, structure prediction).
- Exceed the per-iteration or lifetime credit budget.
- Persist state between iterations other than through the Tier 2 summary and disk cache (no hidden state).

### 9.3 Agent Autonomy: Custom Scripts

The agent is explicitly permitted to write and execute Python scripts during its iterations. This is critical because the most powerful pipelines are likely to combine structural tool outputs (e.g., Rosetta interface scores, StaB-ddG ΔΔG predictions) with sequence features (e.g., ESM-2 embeddings) through a learned model. The agent can use standard Python scientific computing libraries (numpy, scipy, scikit-learn, pandas) to train regressors, compute feature combinations, optimize ensemble weights, or perform any other data analysis. Scripts are sandboxed to the agent's workspace and cannot access external resources or the eval/test set ground truth.

---

## 10. Agent Instruction Sheet (PROGRAM.md)

The following is the initial draft of the instruction sheet provided to the agent at the start of the experiment. This document defines the agent's task, constraints, and available interfaces.

---

### Your Task

You are a computational biology research agent. Your goal is to discover an optimized pipeline of structural biology tools that accurately predicts the relative binding affinity of antibody variants for an antigen, given only the amino acid sequences of the antibody heavy chain, light chain, and antigen.

You will work iteratively: each iteration, you propose a pipeline, evaluate it, and refine based on results. Your performance is measured by Spearman rank correlation (ρ) between your predicted binding scores and experimentally measured affinities on a held-out eval set.

### The System

**Antibody:** CR9114, a broadly neutralizing anti-influenza antibody. The heavy chain has 16 positions where somatic mutations were introduced during affinity maturation. Each variant is defined by a 16-bit genotype string indicating which positions carry the somatic (1) or germline (0) residue. The light chain is constant.

**Antigen:** HA1 (H1N1 influenza hemagglutinin). All affinity measurements are for CR9114 variants binding HA1.

**Affinities:** Experimentally measured KD values ranging from ~150 pM to ~100 nM (about 3 orders of magnitude). You do not need to predict exact KD values — a unitless score that preserves the rank ordering is sufficient.

**Structures:** Two pre-computed structures are available in your workspace: the mature CR9114 (all 16 somatic mutations) in complex with HA1, and the germline CR9114 (all germline residues) in complex with HA1. Both were predicted with Boltz-2. Use these as starting points for in silico mutagenesis.

### Available Tools

You have access to the following autobio tools. Each has a credit cost per structure/sequence/mutation processed:

rosetta-score (1 credit), evoef2 (1 credit), stab-ddg (3 credits), esm-2 (2 credits), esm-1b (1 credit), rosetta-minimize (5 credits), omm-amber-minimize (8 credits), omm-amber-relax (15 credits).

Run `autobio info <tool> --format json` for detailed usage instructions and parameter options. All tools run in containers and accept standardized JSON configs.

### Budget

You have 4,000 credits per iteration and 60,000 credits total across all iterations (max 15 iterations). Credits are spent when tools process structures. Cached results from prior iterations are free. Budget your credits carefully — consider running cheap tools broadly first, then expensive tools selectively on informative subsets.

### Interfaces

**Training data:** Use the query API to retrieve subsets of training data. You can filter by mutation count, specific positions, affinity range, Hamming distance, and more. Each query returns up to 200 records. Training KD values are visible to you. Use `query summary`, `query mutation_frequency`, and `query mutation_impact` to orient yourself before diving in.

**Cache:** Use `cache lookup` to check for existing results before running tools. Always check the cache first to avoid wasting credits on work that's already been done.

**Eval feedback:** After running your pipeline on the eval set, you receive Spearman ρ, top-50 precision, and pairwise accuracy. You do NOT see individual predictions or KD values for the eval set.

### Iteration Protocol

1. Orient: Review your context summary (run log, performance history, your prior insights).
2. Hypothesize: State what you want to try and why.
3. Execute: Query data, check cache, run tools, assemble your scoring pipeline.
4. Train evaluate: Compute metrics against known training KDs. Examine individual errors.
5. Eval evaluate: Submit predictions for the eval set. Receive aggregate metrics.
6. Reflect: Update your insights scratchpad. What worked? What failed? What should you try next?

### Strategic Guidance

- Start simple. Your first iteration should establish a baseline with a single cheap scoring tool (e.g., rosetta-score) on a moderate-sized subset. Understand the data before building complex pipelines.
- Use the query API to understand the landscape. Which mutation positions matter most? Where are the biggest affinity jumps? This knowledge should guide which variants you spend expensive tool credits on.
- Cache results aggressively in early iterations. Broad, cheap scoring builds a foundation you can refine later.
- Consider hybrid approaches. Structure-based scores capture physics; ESM embeddings capture evolutionary patterns. The best pipeline may combine both.
- Pay attention to pairwise accuracy. If your overall ρ is decent but pairwise accuracy is low, your pipeline struggles with fine-grained local predictions — consider more thorough minimization or different scoring.
- The affinity range is only ~3 orders of magnitude. Subtle signal matters. Noisy or poorly calibrated tools may be worse than simpler ones.

---

## 11. Implementation Plan

### 11.1 Pre-Experiment Setup

These steps are completed before the agent begins iterating.

1. **Structure prediction.** Run Boltz-2 via autobio to predict the mature CR9114–HA1 and germline CR9114–HA1 complexes. Validate interface quality. Estimate time: 1–2 hours per structure.
2. **Dataset construction.** Run the stratified downsampling procedure (Section 2.3) on the full H1 binding landscape. Output three partition files: train.csv, eval.csv, final_test.csv. Compute and store the set of informative single-mutation pairs for eval. Estimate time: minutes.
3. **Query API implementation.** Implement the queryable training dataset as a Python module exposing the operations defined in Section 6. The API reads from train.csv and returns formatted text suitable for injection into an LLM context. Estimate time: 1–2 days.
4. **Results cache implementation.** Implement the two-tier cache (Section 7): JSON-based disk store with cache lookup querying, plus the Tier 2 summary generator that produces the context-window summary from the disk store. Estimate time: 1–2 days.
5. **Credit accounting.** Implement the credit tracker as part of the harness. Wraps autobio tool invocations to deduct credits before execution and enforce budget limits. Estimate time: half a day.
6. **Iteration harness.** Build the outer loop that orchestrates agent iterations: injects the Tier 2 summary, presents tools and APIs, captures the agent's pipeline definition, runs eval-set evaluation, computes metrics, and updates the cache and summary. Estimate time: 2–3 days.
7. **Agent system prompt.** Finalize the PROGRAM.md instruction sheet (Section 10) and test it with a dry-run iteration to verify the agent understands its interfaces and constraints.

### 11.2 Execution

1. **Launch.** Start the agent with iteration 1. Monitor the first 2–3 iterations closely to verify the agent interacts correctly with all interfaces.
2. **Autonomous iteration.** Allow the agent to run iterations 4–15 (or until termination conditions are met) with minimal human intervention. Monitor credit consumption and wall-clock time.
3. **Final evaluation.** Upon termination, apply the best pipeline to the final test set. Report Spearman ρ, top-50 precision, and pairwise accuracy.

### 11.3 Post-Experiment Analysis

1. **Pipeline characterization.** Document the agent's best pipeline: which tools it uses, in what order, with what parameters, and how it combines their outputs. This is the primary scientific deliverable.
2. **Ablation study.** Evaluate simplified versions of the best pipeline (dropping tools, changing parameters) to understand which components contribute most.
3. **Agent trajectory analysis.** Review the agent's iteration history: how did its strategy evolve? Did it converge smoothly or exhibit exploratory behavior? Were its insights accurate? This informs future experiments.
4. **Generalization test.** Apply the discovered pipeline to CR6261 variants and/or CR9114 with H3/Flu-B antigens to assess whether the workflow generalizes beyond the training system.

---

## Appendix A: Detailed Affinity Landscape Statistics

The following statistics characterize the full CR9114 × H1 binding landscape (65,094 variants with measured affinities) and inform the downsampling strategy.

### A.1 Affinity Distribution

| Statistic                  | Value                          |
|----------------------------|--------------------------------|
| Total variants             | 65,094                         |
| KD minimum                 | 1.46 × 10⁻¹⁰ M (146 pM)      |
| KD maximum                 | 1.00 × 10⁻⁷ M (100 nM)       |
| KD median                  | 4.31 × 10⁻¹⁰ M (431 pM)      |
| –log₁₀(KD) range          | 7.0 to 9.84                   |
| Strong binders (< 1 nM)   | 44,801 (68.8%)                 |
| Moderate (1–100 nM)        | 18,618 (28.6%)                 |
| Weak (> 100 nM)            | 1,675 (2.6%)                   |

### A.2 Single-Mutation Effects

| Category                                | Count                                     |
|-----------------------------------------|-------------------------------------------|
| Pairs with > 10-fold KD change          | 28,474                                    |
| Pairs with > 100-fold KD change         | 2,258                                     |
| Most impactful position                 | 10 (S75F): 12,925 large-effect pairs      |
| Second most impactful                   | 4 (I52S): 5,365 large-effect pairs        |
| Third most impactful                    | 9 (K74I): 4,973 large-effect pairs        |
| Largest single effect observed          | 499-fold (position 10, 200 pM → 100 nM)  |

### A.3 Mean Affinity by Mutation Count

| Mutations       | N Variants | Mean –log₁₀(KD) | Interpretation                            |
|-----------------|------------|------------------|-------------------------------------------|
| 0 (germline)    | 1          | 8.42             | Moderate binder (3.8 nM)                  |
| 1–4             | 2,482      | 8.62–8.68        | Slightly improved over germline           |
| 5–8             | 36,353     | 8.77–9.04        | Good binders; bulk of the library         |
| 9–12            | 25,561     | 9.14–9.40        | Strong binders                            |
| 13–16 (mature)  | 698        | 9.47–9.59        | Strongest binders (~250–400 pM)           |
