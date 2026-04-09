# SAAINT-DB Data Processing

Data preparation pipeline for the SAAINT-DB antibody-antigen binding affinity experiment.
Produces quality-filtered, deduplicated, antigen-split datasets suitable for an autonomous
agent experiment predicting relative binding affinity from structural/interface features.

Run once: `python prepare.py`

---

## 1. Raw Data

Two TSV files in `data/`, both from the SAAINT-DB (Structural Antibody-Antigen Interaction
Database):

### `saaintdb_affinity_all.tsv` (6,158 rows)

Affinity measurements linked to PDB structures.

| Column | Type | Description |
|--------|------|-------------|
| `PDB_ID` | str | PDB accession (e.g., `5zxv`) |
| `PMID` | str | PubMed ID of source publication |
| `DOI` | str | Digital Object Identifier |
| `Model_index` | int | Model index in structure (almost always 0) |
| `Asym_ID_type` | str | Chain ID convention (`auth_asym_id` or `label_asym_id`) |
| `H_chain_ID` | str | Heavy chain ID in PDB |
| `L_chain_ID` | str | Light chain ID (`N.A.` for VHH/nanobodies) |
| `Ag_chain_ID(s)` | str | Antigen chain(s), semicolon-separated for multi-chain |
| `Ag_type(s)` | str | Antigen type (`protein`, `peptide`) |
| `Affinity_KD(nM)` | str | Dissociation constant in nanomolar; `N.A.` if unavailable |
| `Affinity_method` | str | Measurement method (`SPR`, `BLI`, `ITC`, `ELISA`, etc.) |
| `Affinity_temp(K)` | str | Temperature in Kelvin; `N.A.` if not recorded |
| `Affinity_notes` | str | Free-text notes |

**Key statistics:**
- 2,707 rows (44%) have affinity data; 3,451 are `N.A.`
- SPR: 1,265 entries; BLI: 1,188 entries (89% of rows with methods)
- Remaining methods: ITC (139), ELISA (62), grating (22), fluorescence (19), etc.
- 69 entries have inequality KD values (`<0.001`, `<10`, `>100000`, etc.)
- 12 entries have asterisk-annotated KD values (e.g., `0.4*`)

### `saaintdb_20260326_all.tsv` (21,627 rows)

Full structural metadata for all antibody-antigen complexes in SAAINT-DB.

| Column | Type | Description |
|--------|------|-------------|
| `PDB_ID` | str | PDB accession |
| `Title` | str | PDB title |
| `Mutation(s)` | str | Reported mutations |
| `Classification` | str | PDB classification |
| `Deposit_date` | str | PDB deposit date |
| `Release_date` | str | PDB release date |
| `Method` | str | Experimental method (x-ray diffraction, etc.) |
| `Resolution` | str | Resolution in angstroms; `N.A.` for NMR/other |
| `R_free` | str | R-free value |
| `R_work` | str | R-work value |
| `PMID` | str | PubMed ID |
| `DOI` | str | DOI |
| `Model_index` | int | Model index |
| `Asym_ID_type` | str | Chain ID convention |
| `Ab_type` | str | Antibody format (`FabH:FabL`, `VH:VL`, `VHH`, `scFv`, `VH`) |
| `H_subgroup` | str | Heavy chain germline subgroup (e.g., `IGHV3`) |
| `L_subgroup` | str | Light chain germline subgroup |
| `H_chain_ID` | str | Heavy chain ID |
| `L_chain_ID` | str | Light chain ID |
| `H_fas_seq` | str | Heavy chain FASTA sequence |
| `L_fas_seq` | str | Light chain FASTA sequence |
| `H_filled_pdb_seq` | str | Heavy chain PDB sequence (gap-filled) |
| `L_filled_pdb_seq` | str | Light chain PDB sequence (gap-filled) |
| `H_mean_radius` | float | Heavy chain mean radius |
| `L_mean_radius` | float | Light chain mean radius |
| `H_fas_seq_len` | int | Heavy chain sequence length |
| `L_fas_seq_len` | int | Light chain sequence length |
| `H_pdb_seq_len` | int | Heavy chain PDB sequence length |
| `L_pdb_seq_len` | int | Light chain PDB sequence length |
| `H_filled_seq_len` | int | Heavy chain filled sequence length |
| `L_filled_seq_len` | int | Light chain filled sequence length |
| `HL_inf_res_num` | int | H-L interface residue count |
| `H_mol_name` | str | Heavy chain molecule name |
| `L_mol_name` | str | Light chain molecule name |
| `H_species` | str | Heavy chain species |
| `L_species` | str | Light chain species |
| `Ag_chain_ID(s)` | str | Antigen chain(s) |
| `Ag_type(s)` | str | Antigen type(s) |
| `Ag_mol_name(s)` | str | Antigen molecule name(s) |
| `Ag_species` | str | Antigen species |
| `Ab_ag_inf_res_num` | int | Ab-Ag interface residue count |
| `CDR_inf_res_num` | int | CDR interface residue count |
| `CDR_inf_res_ratio` | float | CDR / total interface ratio |

---

## 2. Filtering Criteria

### Method filter: SPR and BLI only

Surface Plasmon Resonance (SPR) and Biolayer Interferometry (BLI) are the gold-standard
methods for measuring antibody-antigen binding kinetics. They provide real-time, label-free
KD measurements with quantitative accuracy. Other methods in the dataset are either
semi-quantitative (ELISA, FACS), measure different physical properties (ITC measures
thermodynamics, not kinetics), or have insufficient representation (<20 entries each).

### Resolution filter: <= 3.5 angstroms

The experiment predicts binding affinity from structural/interface features. At resolutions
worse than 3.5A, sidechain conformations (critical for interface contacts) become unreliable.
The 3.5A cutoff retains 97% of crystal structures while excluding the worst-quality data.

### Numeric KD only

Entries with `N.A.`, inequality prefixes (`<`, `>`), or other non-numeric KD values are
dropped. Inequality values (69 entries) represent assay detection limits, not true
measurements. Asterisk-annotated values (e.g., `0.4*`, 12 entries) are treated as valid
after stripping the asterisk.

### Post-filter dataset

| Metric | Value |
|--------|-------|
| After merge (Model_index=0) | 6,098 |
| SPR/BLI only | 2,453 |
| Numeric KD (inequalities dropped) | 2,321 (90 dropped) |
| Resolution <= 3.5 A | 1,745 |
| Unique PDB IDs | 849 |
| KD range | 0.0008 -- 12,900 nM |
| neg_log10_KD (pKD) range | 4.89 to 12.10 |
| Ab types | FabH:FabL (891), VH:VL (408), VHH (291), scFv (105), VH (45), other (5) |

Note: Both source files contain multi-model duplicates (same PDB_ID + chain IDs at
different `Model_index` values). Filtering to `Model_index == 0` before the merge
eliminates these (70 duplicates in the affinity file, 503 in the structural file).

---

## 3. Deduplication

### Within-PDB crystal copies

Many PDB entries contain 2-6 copies of the same antibody-antigen complex in the asymmetric
unit. All copies share identical sequences, identical KD, and differ only in chain letter
assignments. These are not independent data points.

**Strategy**: Group by `(PDB_ID, H_fas_seq, L_fas_seq)`. Within each group, keep the first
entry (sorted by H_chain_ID for determinism). For VHH/nanobodies where `L_fas_seq = "N.A."`,
group by `(PDB_ID, H_fas_seq)` only.

**Result**: 1,745 -> 923 unique complexes. 822 crystal copies removed.

### Antigen name normalization

Antigen names contain semicolon-separated components for multi-chain antigens. Homomultimeric
antigens produce apparent "distinct" entries that are actually the same target:

- `"spike glycoprotein;spike glycoprotein"` -> `"spike glycoprotein"` (homomultimer collapse)
- `"ricin a chain;ricin a chain"` -> `"ricin a chain"`
- `"guanine nucleotide-binding protein g(i) subunit alpha-1;..."` -> kept as-is (heteromultimer)

**Algorithm**: Split on `;`, deduplicate identical components, sort unique components
alphabetically, rejoin with `;`.

### Cross-PDB deduplication

The same antibody can be solved in multiple crystal forms (different PDB entries) bound to
the same antigen. These provide no additional information for training.

**Strategy**: Group by `(H_fas_seq, L_fas_seq, normalized_ag_name)`. Within each group,
keep the entry with the best (lowest) resolution.

Note: the same antibody bound to *different* antigens represents genuinely different
data points and is retained.

**Result**: 923 -> 802 unique complexes (121 cross-PDB duplicates removed), from 752
unique PDB IDs.

---

## 4. SARS-CoV-2 Downsampling

### The problem

After deduplication, SARS-CoV-2 accounts for 42.5% of the dataset (341 of 802 entries):
331 spike entries and 10 non-spike entries (nucleoprotein, 3CLpro).

If left unaddressed, the dataset would be dominated by a single antigen target, biasing
both training and evaluation. However, antibodies to different spike epitopes have
structurally distinct interfaces and should not be treated as redundant.

Only CoV-2 spike entries are clustered and downsampled. Non-spike CoV-2 entries (10 total)
are kept as-is in the "everything else" partition. CoV-2 spike entries without downloadable
structures are dropped (0 in the current run).

### Approach: interface contact clustering

Rather than arbitrary capping, we cluster CoV-2 spike antibodies by their epitope footprint
on the antigen surface, then downsample within each epitope cluster.

#### Step 1: Compute interface contacts

For each CoV-2 spike complex:
1. Parse the mmCIF structure with BioPython's `MMCIFParser`
2. Identify antibody chains (H_chain_ID, L_chain_ID) and antigen chains (Ag_chain_ID(s))
3. Compute all antigen residues with any heavy atom within 4.5A of any antibody heavy atom
4. Record the set of antigen residue positions (using author residue numbering)

SARS-CoV-2 spike structures overwhelmingly use numbering consistent with the full-length
spike protein (UniProt P0DTC2, 1273 residues). This means contact residue positions are
directly comparable across structures without explicit sequence alignment.

**Sanity check**: If a structure's contact residues fall outside the expected range (1-1273),
it may use a renumbered construct. Flag these for review; assign to an "unassigned" cluster.

**Result**: 324/331 CoV-2 spike complexes successfully parsed (97.9%). 7 failures (missing
chains, parse errors) assigned to the "unassigned" cluster.

#### Step 2: Jaccard-based clustering

1. Build a binary contact matrix: rows = complexes, columns = all unique residue positions
2. Compute pairwise Jaccard distance: `J(A,B) = 1 - |A intersect B| / |A union B|`
3. Hierarchical clustering with average linkage
4. Cut dendrogram at distance threshold 0.65

Average linkage balances between single-linkage (which chains distant outliers) and
complete-linkage (which over-splits). A Jaccard distance cutoff of 0.65 means two antibodies
sharing >= 35% of their contact residues are grouped together. This is intentionally loose --
the goal is broad epitope classes, not fine-grained epitope bins.

**Expected clusters** (~5-10 groups):
- RBD class 1/2 (ACE2-competing, top of RBD)
- RBD class 3 (side of RBD)
- RBD class 4 (cryptic, base of RBD)
- NTD (N-terminal domain)
- S2 stem helix
- S2 fusion peptide
- Miscellaneous / cross-domain

If clustering produces <4 groups, decrease the cutoff by 0.05 increments. If >12, increase.
Log cluster counts and sizes for manual review.

**Result**: The cutoff auto-tuned from the initial 0.65 up to 0.90, producing 11 clusters.
The upward adjustment indicates that CoV-2 spike epitope footprints overlap more than
anticipated -- most antibodies contact overlapping RBD regions. Cluster sizes ranged from
1 (singletons) to 165 (dominant RBD epitope group), with three large clusters (165, 70, 54)
accounting for most entries.

#### Step 3: Downsample within clusters

Cap each cluster at 12 representatives. Selection priorities:
1. **Ab type diversity**: Include at least one of each antibody format present in the cluster
2. **Best resolution**: Among remaining candidates, prefer lower resolution values
3. **KD range**: Take entries at the extremes of the affinity range first, then fill from
   the middle

**Target**: Reduce SARS-CoV-2 from ~380 PDBs to ~60-80 total entries.

Entries that failed interface analysis (unparseable structures, missing chains) are assigned
to an "unassigned" cluster, capped at 5 entries.

**Result**: CoV-2 spike reduced from 331 -> 73 entries (7 parse failures, 5 kept as
unassigned). Three large clusters (165, 70, 54 entries) were each capped at 12. Eight
smaller clusters (1-15 entries) were kept in full or nearly in full.

---

## 5. Antigen-Based Train/Eval/Test Splitting

### Rationale

The experiment tests whether structural/interface features generalize to predict binding
affinity for unseen antigen targets. To evaluate this, all entries for a given antigen must
be in the same split. If antibodies to the same antigen appear in both train and test, the
model could exploit antigen-specific features rather than learning generalizable structural
principles.

### Antigen group construction

Each entry is assigned a canonical antigen group from `(normalized_species, normalized_mol_name)`:

- **Species normalization**: Take the first semicolon-component, lowercase, strip strain
  identifiers (e.g., `"plasmodium falciparum (isolate 3d7)"` -> `"plasmodium falciparum"`)
- **Molecule normalization**: Use the homomultimer-collapsed name from Stage 2
- **SARS-CoV-2 override**: All CoV-2 entries (post-downsampling) form a single group
  `"sars-cov-2_spike"` regardless of specific spike domain
- **Unknown antigens**: Entries with `N.A.` species or molecule name form individual
  groups keyed by PDB_ID

### Greedy splitting algorithm

Target ratios: 70% train, 15% eval, 15% test.

1. Collect all antigen groups with entry counts
2. Sort by entry count (descending) -- largest groups assigned first
3. For each group, assign to the split furthest below its target ratio
4. Constraints: eval >= 50 entries, test >= 50 entries
5. Post-hoc: verify KD distribution balance; swap smallest groups if a split is missing
   an entire affinity quartile

**Assertions**:
- No antigen group spans multiple splits
- Eval and test each contain >= 50 entries
- Realized ratios are within 5 percentage points of targets

### Result

After CoV-2 downsampling, the dataset has 544 entries across 267 antigen groups.

| Split | Entries | % | Ag Groups | pKD (min/med/max) | Res median |
|-------|---------|------|-----------|-------------------|------------|
| Train | 381 | 70.0 | 186 | 5.38 / 8.68 / 12.10 | 2.60 A |
| Eval | 82 | 15.1 | 41 | 4.89 / 8.79 / 11.56 | 2.75 A |
| Test | 81 | 14.9 | 40 | 5.71 / 8.09 / 11.77 | 2.71 A |

Ab type distribution across splits:

| Ab type | Train | Eval | Test |
|---------|-------|------|------|
| FabH:FabL | 208 | 46 | 47 |
| VH:VL | 61 | 10 | 22 |
| VHH | 80 | 13 | 9 |
| scFv | 17 | 12 | 3 |
| VH | 12 | 1 | 0 |

Method distribution: SPR and BLI are represented in all splits. Train has
SPR=228/BLI=153; eval SPR=46/BLI=36; test BLI=43/SPR=38.

---

## 6. Structure Download

All unique PDB IDs in the final dataset have their mmCIF files downloaded from RCSB:

```
https://files.rcsb.org/download/{pdb_id}.cif
```

**Design choices**:
- mmCIF over PDB format: large spike+antibody complexes can exceed PDB format's 99,999
  atom limit and 62 chain ID limit
- Downloaded to `structures/{pdb_id}.cif`
- Idempotent: existing files are skipped on re-run
- Failures are logged; entries with missing structures are flagged but not dropped
  (the structure may be fetchable manually)

**Result**: 752/752 structures downloaded successfully (100%). Total download size ~3-5 GB.

---

## 7. Output Format

### Split files: `splits/train.csv`, `splits/eval.csv`, `splits/test.csv`

| Column | Type | Description |
|--------|------|-------------|
| `PDB_ID` | str | PDB accession |
| `H_chain_ID` | str | Heavy chain ID in structure |
| `L_chain_ID` | str | Light chain ID (`N.A.` for VHH) |
| `Ag_chain_ID(s)` | str | Antigen chain(s), semicolon-separated |
| `Ab_type` | str | Antibody format |
| `Ag_species` | str | Antigen species (original, un-normalized) |
| `Ag_mol_name(s)` | str | Antigen molecule name(s) (original) |
| `Resolution` | float | Structure resolution in angstroms |
| `KD_nM` | float | Dissociation constant in nanomolar |
| `neg_log10_KD` | float | -log10(KD in M); higher = tighter binding (pKD) |
| `Affinity_method` | str | `SPR` or `BLI` |
| `antigen_group` | str | Canonical antigen group for splitting |

### Antigen map: `splits/antigen_split_map.csv`

| Column | Type | Description |
|--------|------|-------------|
| `antigen_group` | str | Canonical antigen group |
| `split` | str | `train`, `eval`, or `test` |
| `n_entries` | int | Number of entries in this group |

---

## 8. Prediction Target

**Global ranking via pKD = -log10(KD in M)**. The raw KD values are in nanomolar; the
script converts to molar before taking -log10, producing the standard pKD scale used in
drug discovery and consistent with what structure-based binding energy predictors output.
On this scale, higher values = tighter binding (e.g., pKD 9 = 1 nM, pKD 7 = 100 nM).

The autonomous agent predicts a score for each antibody-antigen complex from
structural/interface features. Evaluation uses Spearman rho (rank correlation) between
predicted scores and measured pKD across all complexes in the eval/test set.

This differs from the pilot experiment, which ranked variants of a single antibody against
a single antigen. Here, the ranking is across diverse antibodies and antigens, testing
whether structural features capture binding strength in a target-agnostic manner.

---

## 9. Script Structure

`prepare.py` follows the pattern established by `experiments/pilot/prepare.py`:

```
constants block (SEED, cutoffs, sizes)
    |
main() {
    Stage 1: merge and filter
    Stage 2: deduplicate
    Stage 3: download structures
    Stage 4: interface analysis (CoV-2)
    Stage 5: CoV-2 clustering and downsampling
    Stage 6: antigen-based splitting
    Stage 7: write output and summary
}
    |
helper functions below main()
```

### Dependencies

Script-level only (not added to `pyproject.toml`):

| Package | Use |
|---------|-----|
| `pandas` | Data loading, merging, filtering |
| `numpy` | Numerics, RNG |
| `biopython` | mmCIF parsing (`MMCIFParser`), neighbor search (`NeighborSearch`) |
| `scipy` | Hierarchical clustering (`linkage`, `fcluster`), Jaccard distance (`pdist`) |
| `requests` | RCSB structure downloads |

### Constants

```python
SEED = 42
RESOLUTION_CUTOFF = 3.5           # angstroms
VALID_METHODS = {"SPR", "BLI"}
CONTACT_DIST_CUTOFF = 4.5         # angstroms, for interface contacts
COV2_SPECIES = "severe acute respiratory syndrome coronavirus 2"
CLUSTER_LINKAGE = "average"
JACCARD_DISTANCE_CUTOFF = 0.65    # initial cutoff for epitope clustering (auto-tuned)
COV2_CAP_PER_CLUSTER = 12         # max entries per epitope cluster
COV2_FAILED_CAP = 5               # max unassigned (parse failure) entries to keep
SPLIT_RATIOS = (0.70, 0.15, 0.15)
MIN_EVAL_ENTRIES = 50
MIN_TEST_ENTRIES = 50
```

---

## 10. Verification

### In-script assertions

- No duplicate `(PDB_ID, H_chain_ID, L_chain_ID)` tuples after deduplication
- No antigen group appears in multiple splits
- Eval >= 50 entries and test >= 50 entries
- No NaN values in `neg_log10_KD`
- All PDB IDs in final output have corresponding structure files

### Distribution checks (printed to stdout)

- KD quartiles per split
- Ab type counts per split
- Resolution distribution per split
- SARS-CoV-2 fraction per split
- Total antigen group count per split

### Reproducibility

Running the script twice with the same inputs produces identical outputs. Achieved via:
- Fixed random seed (`np.random.default_rng(SEED)`)
- Deterministic sort order in all groupby/dedup operations
- Deterministic greedy splitting (sorted antigen groups, tie-breaking by name)

### All assertions passed

All in-script assertions passed on the production run (2026-04-09):
- Zero duplicate keys after dedup
- Zero antigen groups spanning multiple splits
- Eval (82) and test (81) both exceed the 50-entry minimum
- Zero NaN values in neg_log10_KD
- Zero missing structure files

---

## 11. Implementation Notes

### Multi-model duplicates

Both source TSV files contain duplicate `(PDB_ID, H_chain_ID, L_chain_ID)` rows caused
by different `Model_index` values (70 in the affinity file, 503 in the structural file).
These come from NMR ensembles or alternate biological assemblies. The script filters to
`Model_index == 0` before merging to avoid creating spurious many-to-many join rows.

### Data partitioning for CoV-2 downsampling

The dataset is split into two mutually exclusive partitions before clustering:

1. **CoV-2 spike with structure** -- sent through interface analysis, clustering, and
   downsampling
2. **Everything else** (non-CoV-2 + CoV-2 non-spike) -- kept as-is

After downsampling, the two partitions are concatenated back together. CoV-2 spike
entries without downloadable structures are dropped entirely (0 in the current run).

### Clustering cutoff auto-tuning

The Jaccard distance cutoff auto-tuned upward from the initial 0.65 to 0.90 to achieve
the target range of 4-12 clusters. This indicates most CoV-2 spike antibodies in the
dataset target overlapping epitopes on the RBD, which is consistent with the literature --
the pandemic-era structural biology effort was heavily focused on RBD-targeting
neutralizing antibodies. A higher cutoff (0.90 = requiring only 10% contact overlap to
be grouped) was needed to separate these into distinct clusters.

---

## 12. Pipeline Summary

```
Raw data (6,158 affinity + 21,627 structural rows)
  |
  v
Stage 1: Merge + quality filters ................ -> 1,745 entries
  |  (Model_index=0, SPR/BLI, numeric KD, resolution <= 3.5 A)
  v
Stage 2: Deduplicate ............................ -> 802 entries
  |  (within-PDB crystal copies, cross-PDB same-antibody)
  v
Stage 3: Download structures .................... 752 mmCIF files (100% success)
  |
  v
Stage 4: Interface analysis (CoV-2 spike) ....... 324/331 parsed (97.9%)
  |
  v
Stage 5: Epitope clustering + downsampling ...... CoV-2: 331 -> 73
  |  (11 clusters at Jaccard cutoff 0.90)
  v
Stage 6: Antigen-based splitting ................ 267 antigen groups
  |  (greedy assignment, 70/15/15 target)
  v
Stage 7: Output ................................. 544 total entries
                                                   train=381, eval=82, test=81
```
