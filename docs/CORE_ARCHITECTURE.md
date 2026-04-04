# Core Architecture: Experiment Harness

**Technical Design for Credit System, Data API, and Results Cache**

Version 0.1 | April 2026

---

## 1. Design Principles

Five principles guide the architecture of the experiment harness:

1. **Experiment-agnostic core.** The framework knows about variants, tools, credits, and
   caches — not about antibodies, genotypes, or binding affinities. All experiment-specific
   semantics live in configuration and pluggable components.
2. **Python-first, CLI-adapted.** Every capability is a Python API call. The Typer CLI is a
   thin translation layer that parses text commands and formats output for LLM or human
   consumption. No logic lives in the CLI layer.
3. **Human-readable storage.** CSV is the primary data format — readable by humans, efficiently
   queryable via DuckDB, and interoperable with pandas, R, and spreadsheet tools.
4. **Budget as a first-class constraint.** The credit system is not an afterthought bolted onto
   tool execution. It is a pre-condition: every tool invocation must pass budget validation
   before it runs.
5. **Build for today, design for tomorrow.** Each component has a clear extension point
   (registries, protocol classes, config-driven behavior) but ships only what the current
   pilot requires.

---

## 2. System Architecture

### 2.1 Component Overview

```
┌─────────────────────────────────────────────────────────────────────┐
│                         Experiment Config                           │
│                          (YAML + Pydantic)                          │
└──────────────┬──────────────────────────────────────────────────────┘
               │ configures
               ▼
┌──────────────────────────────────────────────────────────────────────┐
│                          Harness Core                                │
│                                                                      │
│  ┌──────────────┐  ┌───────────────┐  ┌──────────────────────────┐  │
│  │ Credit Ledger │  │ Variant Store │  │     Results Cache        │  │
│  │              │  │  (Data API)   │  │  (Tier 1 + Tier 2)       │  │
│  └──────┬───────┘  └───────┬───────┘  └────────────┬─────────────┘  │
│         │                  │                        │                │
│  ┌──────┴───────┐  ┌──────┴────────┐  ┌───────────┴─────────────┐  │
│  │ Tool Registry │  │  Evaluation   │  │   Context Summary       │  │
│  │              │  │   Registry    │  │   Generator             │  │
│  └──────────────┘  └───────────────┘  └─────────────────────────┘  │
│                                                                      │
└──────────────────────────────────────────────────────────────────────┘
               │
               ▼
┌──────────────────────────────────────────────────────────────────────┐
│                         CLI Layer (Typer)                            │
│                                                                      │
│   autoimmune query ...    autoimmune cache ...    autoimmune run ... │
│   autoimmune budget ...   autoimmune eval ...     autoimmune info ..│
└──────────────────────────────────────────────────────────────────────┘
```

### 2.2 Data Flow

A single tool invocation follows this path:

```
Agent request (Python call or CLI command)
  │
  ├─ 1. Cache check: has (tool, params, variant_ids) been computed?
  │     ├─ All cached → return results, zero credits
  │     └─ Some/none cached → continue with uncached subset
  │
  ├─ 2. Budget check: does remaining budget ≥ cost(tool, n_uncached)?
  │     ├─ Yes → proceed
  │     └─ No  → raise BudgetExhausted with details
  │
  ├─ 3. Execute: invoke the tool's registered executor on uncached variants
  │
  ├─ 4. Record: write results to cache (Tier 1), debit ledger
  │
  └─ 5. Return: merge cached + fresh results, return to caller
```

---

## 3. Configuration

### 3.1 Schema

All experiment-specific parameters live in a single YAML file, validated at load time by a
Pydantic model. The schema is designed to be forward-compatible with OmegaConf for future
features (defaults, interpolation, CLI overrides).

```yaml
# experiment.yaml — CR9114 × H1 pilot example

experiment:
  name: cr9114_h1_pilot
  description: >
    Autonomous pipeline discovery for CR9114 heavy-chain variant
    binding affinity prediction against H1 influenza hemagglutinin.

dataset:
  train: data/train.csv
  eval: data/eval.csv
  test: data/final_test.csv

  # Schema: defines which columns exist and their roles.
  # The framework uses these roles for generic operations.
  schema:
    variant_id: genotype          # column used as the unique variant identifier
    target: kd                    # column holding the ground-truth metric
    target_transform: neg_log10   # transformation applied before metric computation
    sequence_columns:             # columns containing sequence data (for tool input)
      - heavy_chain_seq
      - light_chain_seq
      - antigen_seq
    feature_columns:              # columns representing individual mutable features
      - pos_29                    # each is 0/1 in the pilot; could be categorical
      - pos_30                    # in future experiments
      - pos_31
      # ... (all 16 positions for the pilot)
    metadata_columns:             # other columns available for filtering
      - n_mutations
      - region                    # e.g., CDR1, FR3, etc.

  query:
    record_limit: 200             # hard cap on records per query response

credits:
  per_iteration: 4000
  lifetime: 60000
  max_iterations: 15

tools:
  - name: rosetta-score
    credit_cost: 1
    cost_unit: structure          # credits = credit_cost × n_structures
    executor: autobio             # executor backend identifier
    default_params:
      score_function: ref2015

  - name: evoef2
    credit_cost: 1
    cost_unit: structure
    executor: autobio

  - name: stab-ddg
    credit_cost: 3
    cost_unit: mutation
    executor: autobio

  - name: esm-2
    credit_cost: 2
    cost_unit: sequence
    executor: autobio
    default_params:
      layer: -1

  - name: esm-1b
    credit_cost: 1
    cost_unit: sequence
    executor: autobio

  - name: rosetta-minimize
    credit_cost: 5
    cost_unit: structure
    executor: autobio

  - name: omm-amber-minimize
    credit_cost: 8
    cost_unit: structure
    executor: autobio

  - name: omm-amber-relax
    credit_cost: 15
    cost_unit: structure
    executor: autobio

evaluation:
  metrics:
    - spearman_rho
    - top_k_precision
    - pairwise_accuracy
  top_k: 50                       # k value for top-k precision
  pair_distance: 1                # Hamming distance for pairwise accuracy pairs

summary:
  token_budget: 2000              # total Tier 2 budget
  run_log_budget: 600             # Section A allocation
  performance_budget: 400         # Section B allocation
  insights_budget: 500            # Section C allocation
  # ~500 tokens reserved for headers, separators, and breathing room

workspace:
  cache_dir: workspace/cache
  structures_dir: workspace/structures
  scripts_dir: workspace/scripts
  insights_file: workspace/insights.md
```

### 3.2 Pydantic Models

```python
# src/autoimmune/config.py

class DatasetSchema(BaseModel):
    variant_id: str
    target: str
    target_transform: str | None = None
    sequence_columns: list[str] = []
    feature_columns: list[str] = []
    metadata_columns: list[str] = []

class DatasetConfig(BaseModel):
    train: Path
    eval: Path
    test: Path
    schema_: DatasetSchema = Field(alias="schema")
    query: QueryConfig = QueryConfig()

class ToolConfig(BaseModel):
    name: str
    credit_cost: int | float
    cost_unit: CostUnit              # StrEnum: structure, sequence, mutation
    executor: str
    default_params: dict[str, Any] = {}

class CreditConfig(BaseModel):
    per_iteration: int
    lifetime: int
    max_iterations: int

class EvalConfig(BaseModel):
    metrics: list[str]
    top_k: int = 50
    pair_distance: int = 1

class SummaryConfig(BaseModel):
    token_budget: int = 2000
    run_log_budget: int = 600
    performance_budget: int = 400
    insights_budget: int = 500

class ExperimentConfig(BaseModel):
    experiment: ExperimentMeta
    dataset: DatasetConfig
    credits: CreditConfig
    tools: list[ToolConfig]
    evaluation: EvalConfig
    summary: SummaryConfig = SummaryConfig()
    workspace: WorkspaceConfig = WorkspaceConfig()
```

Validation rules enforced at load time:

- All referenced CSV paths must exist.
- `schema.variant_id` and `schema.target` must be columns in the train CSV.
- Tool names must be unique.
- `per_iteration × max_iterations` should be ≥ `lifetime` (warn if not — it means the
  iteration cap is unreachable).
- All metric names in `evaluation.metrics` must be registered in the metric registry.

---

## 4. Credit System

### 4.1 Concepts

The credit system has three entities:

| Entity         | Description                                                       |
|----------------|-------------------------------------------------------------------|
| **CreditLedger** | The stateful tracker. Holds the budget, records transactions, enforces limits. |
| **CostUnit**     | An enum (`structure`, `sequence`, `mutation`) defining what the tool charges per. |
| **LedgerEntry**  | An immutable record of a single debit: tool, params, variant count, cost, timestamp, iteration. |

### 4.2 CreditLedger

```python
class CreditLedger:
    """Tracks credit expenditure against per-iteration and lifetime budgets."""

    def __init__(self, config: CreditConfig) -> None: ...

    # --- State ---
    @property
    def iteration(self) -> int: ...

    @property
    def iteration_remaining(self) -> int: ...

    @property
    def lifetime_remaining(self) -> int: ...

    @property
    def lifetime_spent(self) -> int: ...

    # --- Operations ---
    def estimate_cost(self, tool: ToolConfig, n_items: int) -> int:
        """Compute cost without spending. Pure function."""

    def can_afford(self, tool: ToolConfig, n_items: int) -> bool:
        """Check whether both iteration and lifetime budgets allow this cost."""

    def debit(
        self,
        tool: ToolConfig,
        n_items: int,
        *,
        run_id: str,
        variant_ids: list[str],
        params: dict[str, Any],
    ) -> LedgerEntry:
        """Deduct credits. Raises BudgetExhausted if insufficient.

        Returns the LedgerEntry for the transaction.
        """

    def advance_iteration(self) -> None:
        """Move to the next iteration. Resets per-iteration counter."""

    # --- Queries ---
    def entries(
        self,
        *,
        iteration: int | None = None,
        tool: str | None = None,
    ) -> list[LedgerEntry]:
        """Filter the transaction log."""

    def summary(self) -> CreditSummary:
        """Snapshot of current budget state for display."""
```

### 4.3 Cost Computation

Cost is always:

```
cost = tool.credit_cost × n_items
```

where `n_items` is the count of the unit specified by `tool.cost_unit`. For multi-step
pipelines (e.g., minimize → score), costs accumulate across steps — the harness debits each
step independently.

The `cost_unit` enum exists so the framework can present meaningful messages ("this will cost
1,500 credits for 300 structures at 5 credits/structure") without understanding what a
"structure" is.

### 4.4 Budget Enforcement

Budget enforcement is **pre-execution**: the harness calls `can_afford()` before invoking any
tool. If the budget is insufficient, execution never starts and the agent receives an error
with its remaining budget and the cost that was attempted. This prevents partial runs that
waste credits.

Two budgets are checked independently:

- **Per-iteration budget** resets at each `advance_iteration()` call.
- **Lifetime budget** is monotonically decreasing and never resets.

A tool invocation is rejected if **either** budget would be exceeded.

### 4.5 Persistence

The ledger writes its transaction log to a CSV file (`workspace/cache/_ledger.csv`) after
every `debit()` call. On harness restart, the ledger reconstructs its state from this file.
Columns:

```
run_id, iteration, tool, params_hash, n_items, cost, timestamp, variant_ids_file
```

The `variant_ids_file` column points to a sidecar file listing the variant IDs processed in
that transaction (avoids unbounded row width in the ledger CSV).

### 4.6 Extension Points

- **Dynamic cost functions.** The current model is `cost × n_items`. A future version could
  allow tools to register a callable `cost_fn(params, n_items) -> int` for parameter-dependent
  pricing (e.g., minimization cost scaling with step count).
- **Cost warnings.** Configurable thresholds that warn (but don't block) when a single tool
  invocation exceeds N% of the iteration budget.

---

## 5. Queryable Training Data API (VariantStore)

### 5.1 Architecture

```
┌──────────────┐          ┌──────────────┐
│  CSV files   │◄─────────│ DuckDB       │
│  (on disk)   │  reads   │ (in-process) │
└──────────────┘          └──────┬───────┘
                                 │ query results
                                 ▼
                          ┌──────────────┐
                          │ VariantStore │
                          │ (Python API) │
                          └──────┬───────┘
                                 │
                          ┌──────┴───────┐
                          ▼              ▼
                    Python callers   CLI (Typer)
```

The `VariantStore` loads CSV files via DuckDB's native CSV reader. DuckDB runs in-process
(no server) and can query CSV files directly without import — the CSV files remain the
source of truth on disk.

### 5.2 VariantStore API

```python
class VariantStore:
    """Queryable interface to experiment training data."""

    def __init__(self, config: DatasetConfig) -> None:
        """Load train/eval/test CSVs. Eval and test are loaded with
        restricted column visibility (no target column)."""

    # --- Summary Queries ---

    def summary(self) -> DatasetSummary:
        """Aggregate stats: count, target distribution (min/max/median/
        quartiles/mean), feature-column distributions, cached result counts."""

    def target_distribution(self, *, bins: int = 20) -> Histogram:
        """Histogram of (transformed) target values."""

    def feature_frequency(self) -> list[FeatureStats]:
        """Per-feature-column: value counts, mean target for each
        value. Generalizes 'mutation_frequency' from the pilot."""

    # --- Filtered Retrieval ---

    def query_variants(
        self,
        *,
        filters: list[ColumnFilter] | None = None,
        sort_by: str | None = None,
        sort_desc: bool = False,
        limit: int | None = None,
    ) -> list[VariantRecord]:
        """Retrieve training variants matching filter criteria.

        Limit is capped at config.query.record_limit regardless of
        the value passed here.
        """

    # --- Comparative Queries ---

    def neighbors(
        self,
        reference_id: str,
        *,
        distance: int = 1,
        feature_columns: list[str] | None = None,
    ) -> list[VariantRecord]:
        """Variants within Hamming distance on feature columns.

        Uses the feature_columns from config.schema by default.
        Generalizes the 'neighbors' query — works on any set of
        categorical feature columns, not just binary genotypes.
        """

    def feature_impact(
        self,
        feature_column: str,
        *,
        stratify_by: str | None = None,
    ) -> FeatureImpactReport:
        """Marginal effect of a single feature on the target.

        Generalizes 'mutation_impact': for each distinct value of
        feature_column, reports mean/median target, optionally
        stratified by another column (e.g., mutation count).
        """

    def feature_interaction(
        self,
        columns: list[str],
    ) -> InteractionTable:
        """N-way interaction table of mean target for all value
        combinations across the specified feature columns.

        Generalizes 'epistasis': works for any 2+ feature columns
        with any number of distinct values, not just binary pairs.
        """

    # --- Sequence Retrieval ---

    def get_sequences(
        self,
        variant_ids: list[str],
        *,
        format: SequenceFormat = SequenceFormat.DICT,
    ) -> list[SequenceRecord]:
        """Retrieve full sequences for specified variants.

        format controls output: DICT (Python), FASTA (text), or
        CSV. Only returns columns listed in schema.sequence_columns.
        """

    # --- Eval/Test Accessors ---

    def eval_variant_ids(self) -> list[str]:
        """Return variant IDs (and metadata, but NOT target values)
        for the eval set."""

    def test_exists(self) -> TestSetInfo:
        """Confirm test set exists, return its size. No data revealed."""
```

### 5.3 ColumnFilter

Filters are the generic mechanism for variant retrieval. Each filter targets a single column
and applies a typed predicate.

```python
@dataclass
class ColumnFilter:
    column: str
    op: FilterOp          # eq, ne, gt, ge, lt, le, between, in_, contains
    value: Any            # scalar or tuple (for between) or list (for in_)
```

The `query_variants` method translates filters into a DuckDB SQL WHERE clause. This means
the framework doesn't need to know about antibody-specific concepts like "mutation count" —
if it's a column in the CSV, it's filterable.

Example: the pilot's `--n-mutations 3:5 --kd-range 1e-10:1e-8 --mutations-include 10` becomes:

```python
store.query_variants(
    filters=[
        ColumnFilter("n_mutations", FilterOp.BETWEEN, (3, 5)),
        ColumnFilter("kd", FilterOp.BETWEEN, (1e-10, 1e-8)),
        ColumnFilter("pos_75", FilterOp.EQ, 1),  # "position 10" is the pos_75 column
    ],
    sort_by="kd",
    limit=100,
)
```

### 5.4 DuckDB Integration

DuckDB is used exclusively as a **query engine**, not as a persistent database. The pattern:

```python
import duckdb

conn = duckdb.connect()  # in-memory, no file
conn.execute("CREATE VIEW train AS SELECT * FROM read_csv_auto(?)", [str(csv_path)])
# All subsequent queries run against this view.
```

This gives us SQL's full filtering, aggregation, and sorting power while the CSV files remain
the sole persistent representation. DuckDB's CSV reader handles type inference, but the
schema config provides explicit column roles so the framework doesn't depend on inference.

### 5.5 Record Limit Enforcement

The hard cap (`config.query.record_limit`, default 200) is enforced at the `VariantStore`
level, not the CLI level. Even Python callers are subject to it. This is intentional: the
limit exists to prevent context-window overflow when results are injected into an LLM prompt,
and that concern applies regardless of the caller. Callers that need bulk access (e.g., the
evaluation system computing metrics across the full training set) use a separate internal
method (`_query_uncapped`) that is not exposed through the CLI.

### 5.6 Partition Isolation

The `VariantStore` enforces data partitioning rules:

- **Train**: full access (all columns, including target).
- **Eval**: variant IDs and metadata columns only. Target column is never loaded from the eval
  CSV into DuckDB views accessible through public methods.
- **Test**: only existence and size are queryable until the harness explicitly unlocks it at
  experiment termination.

This is enforced at the DuckDB view level: the eval view simply omits the target column.

---

## 6. Results Cache

### 6.1 Two-Tier Design

```
Tier 1 (Disk Store)                      Tier 2 (Context Summary)
┌─────────────────────┐                  ┌──────────────────────────┐
│ CSV files per        │  ──generates──▶ │ Compressed text summary  │
│ (tool, params_hash)  │                 │ injected into agent      │
│                      │                 │ context at each iteration│
│ + _runs.csv log      │                 │                          │
│ + _ledger.csv        │                 │ Sections:                │
│                      │                 │  A: Run log              │
└─────────────────────┘                  │  B: Performance tracker  │
                                         │  C: Agent insights       │
                                         └──────────────────────────┘
```

### 6.2 Tier 1: Disk Store

#### Storage Layout

```
workspace/cache/
├── _runs.csv                          # run log (one row per tool invocation)
├── _ledger.csv                        # credit transactions
├── _performance.csv                   # eval metrics per iteration
├── rosetta_score__ref2015.csv         # results: one file per (tool, params_hash)
├── rosetta_score__beta_nov16.csv
├── stab_ddg__default.csv
├── esm2__layer_neg1.csv
├── omm_amber_min__steps500.csv
└── ...
```

Each results CSV has columns:

```
variant_id, <tool-specific result columns...>, run_id, iteration
```

The tool-specific columns depend on the tool (e.g., `total_score, interface_score,
per_residue_json` for rosetta-score). The framework does not interpret these columns — it
stores and retrieves them opaquely. Tools declare their output schema when registering.

The `_runs.csv` log has columns:

```
run_id, iteration, tool, params_hash, params_json, n_variants, credits_spent,
wall_time_seconds, timestamp
```

#### Cache Key

The cache key for a result is the tuple `(tool_name, params_hash, variant_id)`. The
`params_hash` is a deterministic hash of the canonical JSON serialization of the tool's
parameters (sorted keys, no whitespace). Two invocations of the same tool with identical
parameters on the same variant are guaranteed to hit cache.

#### ResultsCache API

```python
class ResultsCache:
    """Tier 1 disk store for tool execution results."""

    def __init__(self, cache_dir: Path, tools: list[ToolConfig]) -> None: ...

    def lookup(
        self,
        tool: str,
        params: dict[str, Any],
        variant_ids: list[str],
    ) -> CacheLookupResult:
        """Check cache for existing results.

        Returns a CacheLookupResult with:
          - hits: dict[variant_id, dict] — cached results
          - misses: list[variant_id] — IDs with no cached result
        """

    def store(
        self,
        tool: str,
        params: dict[str, Any],
        results: dict[str, dict[str, Any]],
        *,
        run_id: str,
        iteration: int,
        wall_time: float,
    ) -> None:
        """Write results to the appropriate CSV file and append to _runs.csv."""

    def get_results(
        self,
        tool: str,
        params: dict[str, Any],
        variant_ids: list[str] | None = None,
    ) -> dict[str, dict[str, Any]]:
        """Retrieve cached results. Returns empty dict for misses."""

    def run_log(self, *, iteration: int | None = None) -> list[RunRecord]:
        """Read the _runs.csv log, optionally filtered by iteration."""

    def record_performance(
        self,
        iteration: int,
        pipeline_description: str,
        metrics: dict[str, float],
    ) -> None:
        """Append a row to _performance.csv."""

    def performance_history(self) -> list[PerformanceRecord]:
        """Read the full _performance.csv."""

    def cached_variant_ids(self, tool: str, params: dict[str, Any]) -> set[str]:
        """Return the set of variant IDs with cached results for this tool+params."""

    def cached_tools(self) -> list[CachedToolSummary]:
        """Summary of what's in cache: tool, params, variant count per file."""
```

#### DuckDB for Cache Queries

Like the VariantStore, the cache uses DuckDB to query its CSV files. The `lookup` method runs:

```sql
SELECT * FROM read_csv_auto('workspace/cache/rosetta_score__ref2015.csv')
WHERE variant_id IN (?, ?, ...)
```

This is efficient even for large cache files because DuckDB pushes predicates down into the
CSV scan.

### 6.3 Tier 2: Context Summary Generator

The `ContextSummaryGenerator` reads from Tier 1 and produces a token-budgeted text summary
for injection into the agent's context window at the start of each iteration.

```python
class ContextSummaryGenerator:
    """Generates the compressed context-window summary from Tier 1 data."""

    def __init__(
        self,
        cache: ResultsCache,
        ledger: CreditLedger,
        config: SummaryConfig,
    ) -> None: ...

    def generate(self, iteration: int) -> str:
        """Produce the full Tier 2 summary for the given iteration.

        Combines Sections A, B, and C within token budgets.
        """

    def section_a_run_log(self, iteration: int) -> str:
        """Compact table of all tool runs. Recent iterations get full
        detail; older iterations are collapsed to one line each."""

    def section_b_performance(self) -> str:
        """Pipeline performance history table."""

    def section_c_insights(self) -> str:
        """Read the agent's insights scratchpad from disk."""

    def update_insights(self, text: str) -> None:
        """Write updated insights to the scratchpad file."""
```

#### Token Budget Management

Each section has an allocated token budget from the config. The generator uses a simple
token estimation heuristic (words × 1.3, since most tokenizers average ~1.3 tokens per
English word for technical text). If a section exceeds its budget:

- **Section A (run log):** Collapse oldest iterations to single-line summaries first. If still
  over budget, drop the oldest entries entirely (they're still on disk in `_runs.csv`).
- **Section B (performance):** Keep all rows — this table grows slowly (one row per iteration,
  max 15). If somehow over budget, drop the least-informative columns.
- **Section C (insights):** Truncate to budget with a `[truncated — see workspace/insights.md
  for full text]` marker. The agent is encouraged to self-edit for conciseness.

#### Insights Scratchpad

The agent's insights (Section C) are stored as a plain Markdown file at the path specified
in `config.workspace.insights_file`. The harness:

1. Reads this file at summary generation time and includes it in Tier 2.
2. Provides `update_insights()` for the iteration harness to persist the agent's updated
   insights after the reflect phase.
3. Does **not** version the insights file across iterations — the agent is expected to
   overwrite stale content. The full iteration history is recoverable from `_runs.csv` and
   `_performance.csv`.

---

## 7. Tool Registry

### 7.1 ToolSpec and Registry

```python
@dataclass
class ToolSpec:
    """A registered tool's specification."""
    name: str
    credit_cost: int | float
    cost_unit: CostUnit
    executor: ToolExecutor                 # Protocol class
    default_params: dict[str, Any]
    output_columns: list[str]              # columns this tool produces in cache CSVs
    description: str = ""

class ToolExecutor(Protocol):
    """Protocol for tool execution backends."""

    def execute(
        self,
        tool_name: str,
        params: dict[str, Any],
        variant_ids: list[str],
        sequences: dict[str, SequenceRecord],
        workspace: Path,
    ) -> dict[str, dict[str, Any]]:
        """Run the tool and return results keyed by variant_id."""
        ...

    def validate_params(self, tool_name: str, params: dict[str, Any]) -> None:
        """Raise ValueError if params are invalid for this tool."""
        ...

class ToolRegistry:
    """Registry of available tools for an experiment."""

    def __init__(self) -> None: ...

    def register(self, spec: ToolSpec) -> None: ...

    def get(self, name: str) -> ToolSpec: ...

    def list_tools(self) -> list[ToolSpec]: ...

    def run(
        self,
        tool_name: str,
        variant_ids: list[str],
        *,
        params: dict[str, Any] | None = None,
        cache: ResultsCache,
        ledger: CreditLedger,
        store: VariantStore,
    ) -> ToolRunResult:
        """High-level tool execution with cache check and budget enforcement.

        1. Merge params with tool defaults.
        2. Check cache for existing results.
        3. Validate budget for uncached variants.
        4. Fetch sequences from store for uncached variants.
        5. Execute tool on uncached variants.
        6. Store results in cache, debit ledger.
        7. Return merged (cached + fresh) results.
        """
```

### 7.2 Executor Backends

The pilot uses `autobio` as the executor backend. The `AutobioExecutor` implements the
`ToolExecutor` protocol by shelling out to `autobio run <tool> --config <config.json>
--format json`.

Future backends might include:
- `WebAPIExecutor` — calls external REST APIs, maps responses to the result schema.
- `PythonCallableExecutor` — wraps a Python function directly (useful for lightweight tools
  like feature extractors or custom scoring functions).
- `MockExecutor` — returns synthetic results for testing and development.

The registry doesn't care which backend a tool uses — it only calls the `ToolExecutor`
protocol.

---

## 8. Evaluation System

### 8.1 Metric Protocol

```python
class Metric(Protocol):
    """Protocol for evaluation metrics."""

    @property
    def name(self) -> str: ...

    def compute(
        self,
        predictions: np.ndarray,
        targets: np.ndarray,
        *,
        variant_ids: np.ndarray | None = None,
        store: VariantStore | None = None,
        **kwargs: Any,
    ) -> MetricResult:
        """Compute the metric.

        predictions and targets are aligned arrays of floats.
        variant_ids and store are provided for metrics that need
        access to variant metadata (e.g., pairwise accuracy needs
        to identify neighbor pairs from feature columns).
        """
        ...
```

### 8.2 Built-in Metrics

Three metrics ship with the framework, matching the pilot requirements:

| Metric             | Implementation                                                    |
|--------------------|-------------------------------------------------------------------|
| `spearman_rho`     | `scipy.stats.spearmanr(predictions, targets).statistic`           |
| `top_k_precision`  | Fraction of predicted top-k that overlap with true top-k.         |
| `pairwise_accuracy`| Among pairs of variants differing by Hamming distance `pair_distance` on feature columns, fraction where predicted ordering matches true ordering. |

### 8.3 MetricRegistry

```python
class MetricRegistry:
    """Registry of available evaluation metrics."""

    def __init__(self) -> None:
        self._metrics: dict[str, Metric] = {}

    def register(self, metric: Metric) -> None: ...

    def compute_all(
        self,
        metric_names: list[str],
        predictions: np.ndarray,
        targets: np.ndarray,
        **kwargs: Any,
    ) -> dict[str, MetricResult]: ...
```

The three built-in metrics are auto-registered at import time. New metrics are added by
implementing the `Metric` protocol and calling `registry.register()`.

### 8.4 Evaluation Flow

The harness evaluates in two modes:

- **Train evaluation.** The agent triggers this. Full visibility: individual predictions,
  errors, and per-variant diagnostics are returned. Uses `_query_uncapped` internally to
  access all training targets.
- **Eval evaluation.** The harness triggers this after the agent submits predictions. Only
  aggregate metrics are returned — no individual predictions, no target values. Results are
  recorded in `_performance.csv`.

---

## 9. CLI Layer

### 9.1 Structure

The CLI is a Typer application with subcommand groups mirroring the Python API:

```
autoimmune
├── query                         # VariantStore operations
│   ├── summary
│   ├── distribution
│   ├── features                  # feature_frequency
│   ├── variants                  # filtered retrieval
│   ├── neighbors
│   ├── impact                    # feature_impact
│   ├── interaction               # feature_interaction
│   ├── sequences
│   ├── eval-ids                  # eval_variant_ids
│   └── test-info                 # test_exists
│
├── cache                         # ResultsCache operations
│   ├── lookup
│   ├── runs                      # run log
│   ├── tools                     # cached_tools summary
│   └── performance               # performance history
│
├── budget                        # CreditLedger queries
│   ├── status                    # current budget state
│   ├── estimate                  # estimate cost for a tool invocation
│   └── history                   # transaction log
│
├── run                           # ToolRegistry execution
│   └── <tool-name>               # dynamic subcommand per registered tool
│
├── eval                          # Evaluation
│   ├── train                     # compute metrics on training set
│   └── submit                    # submit eval-set predictions
│
└── info                          # Informational
    ├── tools                     # list available tools with costs
    ├── schema                    # show dataset schema
    └── config                    # show active experiment config
```

### 9.2 Output Formatting

CLI output is formatted as plain text tables by default, designed for LLM consumption:
compact, no ANSI colors (unless a `--rich` flag is passed for human terminals). Numeric
values use consistent formatting (scientific notation for KD, 4 decimal places for metrics).

An optional `--json` flag on all commands outputs machine-readable JSON for programmatic
consumers.

### 9.3 Config Resolution

The CLI resolves the experiment config from (in priority order):

1. `--config <path>` flag
2. `AUTOIMMUNE_CONFIG` environment variable
3. `experiment.yaml` in the current working directory

---

## 10. Module Layout

```
src/autoimmune/
├── __init__.py
├── config.py                    # Pydantic config models, YAML loader
├── types.py                     # Shared types: CostUnit, FilterOp, SequenceFormat, etc.
│
├── credits/
│   ├── __init__.py
│   ├── ledger.py                # CreditLedger
│   └── models.py                # LedgerEntry, CreditSummary, BudgetExhausted
│
├── data/
│   ├── __init__.py
│   ├── store.py                 # VariantStore
│   ├── filters.py               # ColumnFilter, filter-to-SQL translation
│   └── models.py                # VariantRecord, DatasetSummary, Histogram, etc.
│
├── cache/
│   ├── __init__.py
│   ├── store.py                 # ResultsCache (Tier 1)
│   ├── summary.py               # ContextSummaryGenerator (Tier 2)
│   └── models.py                # RunRecord, CacheLookupResult, PerformanceRecord
│
├── tools/
│   ├── __init__.py
│   ├── registry.py              # ToolRegistry, ToolSpec
│   ├── executor.py              # ToolExecutor protocol, ToolRunResult
│   └── autobio.py               # AutobioExecutor (pilot backend)
│
├── evaluation/
│   ├── __init__.py
│   ├── registry.py              # MetricRegistry
│   ├── metrics.py               # Built-in metrics: spearman, top_k, pairwise
│   └── models.py                # MetricResult
│
└── cli/
    ├── __init__.py
    ├── app.py                   # Typer app root, config resolution
    ├── query.py                 # query subcommands
    ├── cache_cmd.py             # cache subcommands (cache.py conflicts with package)
    ├── budget.py                # budget subcommands
    ├── run.py                   # run subcommands
    ├── eval.py                  # eval subcommands
    ├── info.py                  # info subcommands
    └── formatting.py            # Output formatting utilities
```

---

## 11. Dependency Inventory

These are the external dependencies required by the core framework:

| Package     | Role                                        | Notes                          |
|-------------|---------------------------------------------|--------------------------------|
| `pydantic`  | Config validation, data models              | v2+                            |
| `pyyaml`    | YAML config parsing                         | pydantic-settings YAML source  |
| `duckdb`    | Query engine for CSV data and cache files   | In-process, zero config        |
| `typer`     | CLI framework                               | With `rich` extra for help     |
| `rich`      | Terminal formatting, progress bars          | Optional for human output      |
| `numpy`     | Array operations for metrics                |                                |
| `scipy`     | `spearmanr` for rank correlation            | Only stats subpackage used     |

Dev-only additions to existing `[dev]` extras: `duckdb` type stubs if available, plus any
test fixtures.

---

## 12. Open Questions

Items to resolve during implementation:

1. **Params hash stability.** The cache key includes a hash of tool parameters. We need a
   canonical serialization that is stable across Python versions and dict ordering. Proposal:
   `json.dumps(params, sort_keys=True, separators=(",", ":"))` piped through SHA-256,
   truncated to 12 hex chars.

2. **DuckDB connection lifecycle.** Should the `VariantStore` hold a persistent DuckDB
   connection for the process lifetime, or create a fresh connection per query? A persistent
   connection with views is simpler and avoids repeated CSV scans. Proposal: single connection
   per `VariantStore` instance, views created at init.

3. **Concurrent tool execution.** The pilot runs tools sequentially, but future experiments
   may benefit from parallel tool runs (e.g., running rosetta-score and evoef2 simultaneously
   on different variant subsets). The credit system is synchronous; parallel execution would
   require either pre-reserving credits or a lock. Proposal: defer to a future version; keep
   synchronous execution for now.

4. **Token estimation accuracy.** The 1.3 words-per-token heuristic is approximate. For
   tighter budget control, we could integrate `tiktoken` (for OpenAI models) or a similar
   tokenizer. Proposal: use the heuristic for now; add configurable tokenizer support later.

5. **Cache invalidation.** Tool results are assumed deterministic (same tool + params +
   variant = same result). If a tool's implementation changes across versions, cached results
   may be stale. Proposal: include a `tool_version` field in cache records; the harness warns
   (but doesn't invalidate) when the current tool version differs from cached records.
