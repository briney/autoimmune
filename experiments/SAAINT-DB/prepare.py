"""One-time data preparation for the SAAINT-DB experiment.

Produces quality-filtered, deduplicated, antigen-split datasets from the
SAAINT-DB antibody-antigen structural database.  Each split entry has a
high-resolution (<= 3.5 A) crystal structure and an SPR- or BLI-measured
binding affinity.

Run once:  python prepare.py
"""

from __future__ import annotations

import re
import time
import warnings
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import requests
    from Bio.PDB import MMCIFParser
    from Bio.PDB.NeighborSearch import NeighborSearch
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import pdist
except ImportError as e:
    raise SystemExit(
        f"Missing dependency: {e}\nInstall with: pip install biopython scipy requests"
    ) from e

# ── Constants ────────────────────────────────────────────────────────────

SEED = 42

RESOLUTION_CUTOFF = 3.5  # angstroms
VALID_METHODS = {"SPR", "BLI"}

CONTACT_DIST_CUTOFF = 4.5  # angstroms, heavy-atom interface contacts
COV2_SPECIES = "severe acute respiratory syndrome coronavirus 2"
CLUSTER_LINKAGE = "average"
JACCARD_DISTANCE_CUTOFF = 0.65  # for epitope clustering
COV2_CAP_PER_CLUSTER = 12
COV2_FAILED_CAP = 5

SPLIT_RATIOS = (0.70, 0.15, 0.15)  # train, eval, test
MIN_EVAL_ENTRIES = 50
MIN_TEST_ENTRIES = 50

RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download/{}.cif"

OUTPUT_COLUMNS = [
    "PDB_ID",
    "H_chain_ID",
    "L_chain_ID",
    "Ag_chain_ID(s)",
    "Ab_type",
    "Ag_species",
    "Ag_mol_name(s)",
    "Resolution",
    "KD_nM",
    "neg_log10_KD",
    "Affinity_method",
    "antigen_group",
]

# Columns selected from each source file before merging (to avoid conflicts)
_AFF_COLS = [
    "PDB_ID",
    "Model_index",
    "H_chain_ID",
    "L_chain_ID",
    "Ag_chain_ID(s)",
    "Affinity_KD(nM)",
    "Affinity_method",
]
_STRUCT_COLS = [
    "PDB_ID",
    "Model_index",
    "H_chain_ID",
    "L_chain_ID",
    "Ab_type",
    "H_fas_seq",
    "L_fas_seq",
    "Ag_mol_name(s)",
    "Ag_species",
    "Resolution",
]


def main() -> None:
    root = Path(__file__).parent
    data_dir = root / "data"
    structures_dir = root / "structures"
    splits_dir = root / "splits"
    structures_dir.mkdir(exist_ok=True)
    splits_dir.mkdir(exist_ok=True)

    rng = np.random.default_rng(SEED)

    # ── Stage 1: Merge and filter ────────────────────────────────────────

    print("Stage 1: Merge and filter")
    df = _merge_and_filter(data_dir)

    # ── Stage 2: Deduplicate ─────────────────────────────────────────────

    print("\nStage 2: Deduplicate")
    n_before = len(df)
    df = _dedup_within_pdb(df)
    print(f"  Within-PDB crystal copies: {n_before} -> {len(df)}")

    df["normalized_ag_name"] = df["Ag_mol_name(s)"].apply(_normalize_ag_name)

    n_before = len(df)
    df = _dedup_cross_pdb(df)
    print(f"  Cross-PDB same-antibody: {n_before} -> {len(df)}")

    dup_mask = df.duplicated(subset=["PDB_ID", "H_chain_ID", "L_chain_ID"], keep=False)
    assert not dup_mask.any(), "Duplicate (PDB_ID, H_chain, L_chain) keys remain"
    print(f"  Unique PDB IDs: {df['PDB_ID'].nunique()}")

    # ── Stage 3: Download structures ─────────────────────────────────────

    print("\nStage 3: Download structures")
    pdb_ids = sorted(df["PDB_ID"].unique())
    download_ok = _download_structures(pdb_ids, structures_dir)
    df["_structure_ok"] = df["PDB_ID"].map(download_ok)
    n_avail = df["_structure_ok"].sum()
    print(f"  Structures available: {n_avail}/{len(df)}")

    # ── Stage 4: Interface analysis (CoV-2 spike) ────────────────────────

    print("\nStage 4: Interface analysis (CoV-2)")
    cov2_mask = df["Ag_species"].str.lower().str.contains(COV2_SPECIES, na=False)
    spike_mask = df["normalized_ag_name"].str.lower().str.contains("spike", na=False)
    cov2_spike_all = cov2_mask & spike_mask

    # Partition into mutually exclusive groups:
    #   1. CoV-2 spike with structure -> cluster and downsample
    #   2. CoV-2 spike without structure -> dropped (can't cluster)
    #   3. Everything else -> kept as-is
    cov2_spike_df = df[cov2_spike_all & df["_structure_ok"]].copy()
    rest_df = df[~cov2_spike_all].copy()
    n_nostruct = (cov2_spike_all & ~df["_structure_ok"]).sum()
    n_nonspike = (cov2_mask & ~spike_mask).sum()

    print(
        f"  CoV-2 spike (with structure): {len(cov2_spike_df)}, "
        f"non-spike CoV-2: {n_nonspike}, "
        f"spike no-structure (dropped): {n_nostruct}, "
        f"non-CoV-2: {len(rest_df) - n_nonspike}"
    )

    contact_sets = _compute_cov2_interfaces(cov2_spike_df, structures_dir)

    # ── Stage 5: CoV-2 epitope clustering and downsampling ───────────────

    print("\nStage 5: CoV-2 epitope clustering")
    cov2_down = _cluster_and_downsample(cov2_spike_df, contact_sets, rng)

    # Recombine: rest (non-CoV-2 + CoV-2 non-spike) + downsampled spike
    df = pd.concat([rest_df, cov2_down], ignore_index=True)
    df = df.sort_values(["PDB_ID", "H_chain_ID"]).reset_index(drop=True)
    print(f"  Dataset after downsampling: {len(df)} entries")

    # ── Stage 6: Antigen-based splitting ─────────────────────────────────

    print("\nStage 6: Antigen-based splitting")
    df = _make_antigen_groups(df)
    n_groups = df["antigen_group"].nunique()
    print(f"  Antigen groups: {n_groups}")

    train, eval_, test, ag_map = _split_by_antigen(df)

    # Assertions
    all_groups = (
        set(train["antigen_group"]) | set(eval_["antigen_group"]) | set(test["antigen_group"])
    )
    for g in all_groups:
        splits_in = []
        for name, part in [("train", train), ("eval", eval_), ("test", test)]:
            if g in part["antigen_group"].values:
                splits_in.append(name)
        assert len(splits_in) == 1, f"Antigen group '{g}' in multiple splits: {splits_in}"

    assert len(eval_) >= MIN_EVAL_ENTRIES, f"Eval has {len(eval_)} < {MIN_EVAL_ENTRIES}"
    assert len(test) >= MIN_TEST_ENTRIES, f"Test has {len(test)} < {MIN_TEST_ENTRIES}"
    assert not train["neg_log10_KD"].isna().any(), "NaN in train neg_log10_KD"
    assert not eval_["neg_log10_KD"].isna().any(), "NaN in eval neg_log10_KD"
    assert not test["neg_log10_KD"].isna().any(), "NaN in test neg_log10_KD"

    # Structure availability check
    all_pdb_ids = pd.concat([train, eval_, test])["PDB_ID"].unique()
    missing = [p for p in all_pdb_ids if not (structures_dir / f"{p}.cif").exists()]
    if missing:
        print(f"  WARNING: {len(missing)} PDB structures missing: {missing[:5]}...")

    # ── Stage 7: Write output ────────────────────────────────────────────

    print("\nStage 7: Write output")
    for name, part in [("train", train), ("eval", eval_), ("test", test)]:
        out_path = splits_dir / f"{name}.csv"
        part[OUTPUT_COLUMNS].to_csv(out_path, index=False)

    ag_map.to_csv(splits_dir / "antigen_split_map.csv", index=False)

    _print_summary(train, eval_, test, splits_dir)


# ── Stage 1 ──────────────────────────────────────────────────────────────


def _merge_and_filter(data_dir: Path) -> pd.DataFrame:
    """Load both TSVs, inner-join, and apply quality filters."""
    aff = pd.read_csv(
        data_dir / "saaintdb_affinity_all.tsv",
        sep="\t",
        keep_default_na=False,
        usecols=_AFF_COLS,
    )
    struct = pd.read_csv(
        data_dir / "saaintdb_20260326_all.tsv",
        sep="\t",
        keep_default_na=False,
        usecols=_STRUCT_COLS,
    )
    print(f"  Loaded: {len(aff)} affinity, {len(struct)} structural rows")

    # Filter to Model_index 0 to avoid multi-model duplicates
    aff = aff[aff["Model_index"].astype(int) == 0]
    struct = struct[struct["Model_index"].astype(int) == 0]

    # Inner join
    join_keys = ["PDB_ID", "H_chain_ID", "L_chain_ID"]
    df = aff.drop(columns=["Model_index"]).merge(
        struct.drop(columns=["Model_index"]),
        on=join_keys,
        how="inner",
    )
    n_merged = len(df)
    print(f"  After merge (Model_index=0): {n_merged}")

    # Method filter
    df = df[df["Affinity_method"].isin(VALID_METHODS)]
    print(f"  SPR/BLI only: {len(df)}")

    # KD filter: drop N.A. and inequality values
    df = df[df["Affinity_KD(nM)"] != "N.A."]
    inequality_mask = df["Affinity_KD(nM)"].str.match(r"^[<>]")
    n_ineq = inequality_mask.sum()
    df = df[~inequality_mask]
    df["KD_nM"] = df["Affinity_KD(nM)"].str.rstrip("*").astype(float)
    print(f"  Numeric KD: {len(df)} ({n_ineq} inequalities dropped)")

    # Resolution filter
    df = df[df["Resolution"] != "N.A."]
    df["Resolution"] = df["Resolution"].astype(float)
    df = df[df["Resolution"] <= RESOLUTION_CUTOFF]
    print(f"  Resolution <= {RESOLUTION_CUTOFF} A: {len(df)}")

    # Compute prediction target
    # Convert KD from nM to M, then compute -log10(KD_M).
    # This is the standard convention (pKD) and matches the scale that
    # structure-based binding energy predictors produce.
    df["neg_log10_KD"] = -np.log10(df["KD_nM"] * 1e-9)

    # Clean up
    df = df.drop(columns=["Affinity_KD(nM)"]).reset_index(drop=True)
    print(f"  Unique PDB IDs: {df['PDB_ID'].nunique()}")
    return df


# ── Stage 2 ──────────────────────────────────────────────────────────────


def _normalize_ag_name(name: str) -> str:
    """Collapse homomultimer antigen names and sort components.

    ``"spike glycoprotein;spike glycoprotein"`` becomes
    ``"spike glycoprotein"``.
    """
    parts = [p.strip() for p in name.split(";")]
    return ";".join(sorted(set(parts)))


def _dedup_within_pdb(df: pd.DataFrame) -> pd.DataFrame:
    """Remove crystal-copy duplicates within each PDB entry.

    Groups by (PDB_ID, H_fas_seq, L_fas_seq) and keeps the entry with the
    alphabetically-first H_chain_ID for determinism.
    """
    df = df.sort_values(["PDB_ID", "H_chain_ID", "L_chain_ID"]).reset_index(drop=True)
    return df.drop_duplicates(subset=["PDB_ID", "H_fas_seq", "L_fas_seq"], keep="first")


def _dedup_cross_pdb(df: pd.DataFrame) -> pd.DataFrame:
    """Remove same-antibody duplicates across PDB entries.

    If the same antibody (identical H+L sequence) was solved multiple times
    against the same antigen, keep only the highest-resolution structure.
    """
    df = df.sort_values(["Resolution", "PDB_ID"]).reset_index(drop=True)
    return df.drop_duplicates(subset=["H_fas_seq", "L_fas_seq", "normalized_ag_name"], keep="first")


# ── Stage 3 ──────────────────────────────────────────────────────────────


def _download_structures(pdb_ids: list[str], structures_dir: Path) -> dict[str, bool]:
    """Download mmCIF files from RCSB.  Skips files already on disk."""
    results: dict[str, bool] = {}
    to_download: list[str] = []

    for pdb_id in pdb_ids:
        if (structures_dir / f"{pdb_id}.cif").exists():
            results[pdb_id] = True
        else:
            to_download.append(pdb_id)

    cached = len(results)
    print(f"  {cached} cached, {len(to_download)} to download")

    for i, pdb_id in enumerate(to_download, 1):
        url = RCSB_DOWNLOAD_URL.format(pdb_id)
        cif_path = structures_dir / f"{pdb_id}.cif"
        try:
            resp = requests.get(url, timeout=60)
            resp.raise_for_status()
            cif_path.write_bytes(resp.content)
            results[pdb_id] = True
        except (requests.RequestException, OSError) as exc:
            print(f"    FAILED {pdb_id}: {exc}")
            results[pdb_id] = False

        if i % 100 == 0 or i == len(to_download):
            print(f"    {i}/{len(to_download)} downloaded")
        # Brief pause every 50 to be polite to RCSB
        if i % 50 == 0:
            time.sleep(0.5)

    ok = sum(v for v in results.values())
    print(f"  Total: {ok} available, {len(results) - ok} failed")
    return results


# ── Stage 4 ──────────────────────────────────────────────────────────────


def _compute_interface_contacts(
    cif_path: Path,
    h_chain: str,
    l_chain: str,
    ag_chains: list[str],
) -> set[int] | None:
    """Compute antigen residue positions within contact distance of antibody.

    Returns the set of antigen author-residue sequence numbers that have at
    least one heavy atom within ``CONTACT_DIST_CUTOFF`` of an antibody heavy
    atom, or ``None`` on failure.
    """
    parser = MMCIFParser(QUIET=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            structure = parser.get_structure("s", str(cif_path))
        except Exception:
            return None

    model = structure[0]
    available_chains = {c.get_id() for c in model}

    ab_chain_ids = {h_chain}
    if l_chain != "N.A.":
        ab_chain_ids.add(l_chain)
    ag_chain_ids = set(ag_chains)

    # Check that required chains exist
    if not ab_chain_ids.issubset(available_chains):
        return None
    if not ag_chain_ids.intersection(available_chains):
        return None

    # Collect heavy atoms
    ab_atoms: list = []
    ag_atoms: list = []
    for chain in model:
        cid = chain.get_id()
        target = ab_atoms if cid in ab_chain_ids else (ag_atoms if cid in ag_chain_ids else None)
        if target is None:
            continue
        for residue in chain:
            if residue.get_id()[0] == "W":  # skip water
                continue
            for atom in residue:
                if atom.element in ("H", "D"):
                    continue
                target.append(atom)

    if not ab_atoms or not ag_atoms:
        return None

    # Neighbor search: find Ag residues contacted by Ab
    ns = NeighborSearch(ag_atoms)
    contact_residues: set[int] = set()
    for atom in ab_atoms:
        for res in ns.search(atom.coord, CONTACT_DIST_CUTOFF, level="R"):
            seq_id = res.get_id()[1]
            contact_residues.add(seq_id)

    return contact_residues if contact_residues else None


def _compute_cov2_interfaces(
    cov2_df: pd.DataFrame, structures_dir: Path
) -> dict[int, set[int] | None]:
    """Run interface analysis on all CoV-2 spike entries."""
    contact_sets: dict[int, set[int] | None] = {}
    n_ok = 0
    total = len(cov2_df)

    for i, (idx, row) in enumerate(cov2_df.iterrows(), 1):
        cif_path = structures_dir / f"{row['PDB_ID']}.cif"
        ag_chains = [c.strip() for c in row["Ag_chain_ID(s)"].split(";")]
        contacts = _compute_interface_contacts(
            cif_path, row["H_chain_ID"], row["L_chain_ID"], ag_chains
        )
        contact_sets[idx] = contacts
        if contacts is not None:
            n_ok += 1

        if i % 50 == 0 or i == total:
            print(f"  {i}/{total} parsed ({n_ok} with contacts)")

    fail_rate = (total - n_ok) / total * 100 if total else 0
    print(f"  Interface contacts: {n_ok}/{total} succeeded ({fail_rate:.1f}% failed)")
    if fail_rate > 10:
        print("  WARNING: >10% of CoV-2 structures failed interface analysis")
    return contact_sets


# ── Stage 5 ──────────────────────────────────────────────────────────────


def _cluster_and_downsample(
    cov2_df: pd.DataFrame,
    contact_sets: dict[int, set[int] | None],
    rng: np.random.Generator,
) -> pd.DataFrame:
    """Cluster CoV-2 entries by epitope footprint and downsample."""
    valid_idxs = [k for k, v in contact_sets.items() if v is not None]
    failed_idxs = [k for k, v in contact_sets.items() if v is None]

    if len(valid_idxs) < 2:
        print("  WARNING: <2 valid contact sets — returning all CoV-2 entries")
        return cov2_df

    # Build binary contact matrix
    all_positions = sorted(set().union(*(contact_sets[k] for k in valid_idxs)))
    pos_to_col = {p: i for i, p in enumerate(all_positions)}

    matrix = np.zeros((len(valid_idxs), len(all_positions)), dtype=bool)
    for row_i, idx in enumerate(valid_idxs):
        for pos in contact_sets[idx]:
            matrix[row_i, pos_to_col[pos]] = True

    # Jaccard distance + hierarchical clustering
    dist_vec = pdist(matrix, metric="jaccard")
    linkage_matrix = linkage(dist_vec, method=CLUSTER_LINKAGE)

    # Auto-tune cutoff for reasonable cluster count
    cutoff = JACCARD_DISTANCE_CUTOFF
    labels = fcluster(linkage_matrix, t=cutoff, criterion="distance")
    n_clusters = len(set(labels))

    while n_clusters < 4 and cutoff > 0.30:
        cutoff -= 0.05
        labels = fcluster(linkage_matrix, t=cutoff, criterion="distance")
        n_clusters = len(set(labels))
    while n_clusters > 12 and cutoff < 0.90:
        cutoff += 0.05
        labels = fcluster(linkage_matrix, t=cutoff, criterion="distance")
        n_clusters = len(set(labels))

    print(f"  Jaccard cutoff: {cutoff:.2f} -> {n_clusters} clusters")

    idx_to_label = dict(zip(valid_idxs, labels.tolist(), strict=True))

    # Downsample each cluster
    selected: list[pd.DataFrame] = []
    for cl in sorted(set(labels)):
        cl_idxs = [i for i, lab in idx_to_label.items() if lab == cl]
        cl_df = cov2_df.loc[cl_idxs]
        kept = _select_representatives(cl_df, COV2_CAP_PER_CLUSTER, rng)
        selected.append(kept)
        print(f"    Cluster {cl}: {len(cl_df)} -> {len(kept)}")

    # Include a few failed entries
    if failed_idxs:
        failed_df = cov2_df.loc[failed_idxs].sort_values("Resolution")
        n_keep = min(len(failed_df), COV2_FAILED_CAP)
        selected.append(failed_df.head(n_keep))
        print(f"    Unassigned (parse failures): {len(failed_idxs)} -> {n_keep}")

    result = pd.concat(selected, ignore_index=True)
    print(f"  CoV-2 spike: {len(cov2_df)} -> {len(result)} after downsampling")
    return result


def _select_representatives(df: pd.DataFrame, cap: int, rng: np.random.Generator) -> pd.DataFrame:
    """Select up to *cap* diverse representatives from *df*.

    Priority: (1) Ab-type diversity, (2) KD-range extremes, (3) best
    resolution.
    """
    if len(df) <= cap:
        return df

    selected_idxs: list = []
    remaining = set(df.index)

    # 1) One representative per Ab type (best resolution of that type)
    for ab_type in sorted(df["Ab_type"].unique()):
        if len(selected_idxs) >= cap:
            break
        candidates = df.loc[list(remaining)]
        typed = candidates[candidates["Ab_type"] == ab_type]
        if typed.empty:
            continue
        pick = typed.sort_values(["Resolution", "PDB_ID"]).index[0]
        selected_idxs.append(pick)
        remaining.discard(pick)

    # 2) KD extremes
    if len(selected_idxs) < cap and remaining:
        cand = df.loc[list(remaining)]
        best_kd = cand["KD_nM"].idxmin()
        selected_idxs.append(best_kd)
        remaining.discard(best_kd)

    if len(selected_idxs) < cap and remaining:
        cand = df.loc[list(remaining)]
        worst_kd = cand["KD_nM"].idxmax()
        selected_idxs.append(worst_kd)
        remaining.discard(worst_kd)

    # 3) Fill by resolution
    if len(selected_idxs) < cap and remaining:
        cand = df.loc[list(remaining)].sort_values(["Resolution", "PDB_ID"])
        n_need = cap - len(selected_idxs)
        selected_idxs.extend(cand.index[:n_need].tolist())

    return df.loc[selected_idxs]


# ── Stage 6 ──────────────────────────────────────────────────────────────


def _make_antigen_groups(df: pd.DataFrame) -> pd.DataFrame:
    """Assign a canonical ``antigen_group`` label to each entry."""
    groups: list[str] = []
    for _, row in df.iterrows():
        ag_species: str = row["Ag_species"]
        ag_mol: str = row.get("normalized_ag_name", row["Ag_mol_name(s)"])

        # CoV-2 override
        if COV2_SPECIES in ag_species.lower():
            groups.append("sars-cov-2_spike")
            continue

        # Unknown antigen
        if ag_species == "N.A." or ag_mol == "N.A.":
            groups.append(f"unknown__{row['PDB_ID']}")
            continue

        # Normalise species: first ;-component, lowercase, strip strain info
        sp = ag_species.split(";")[0].strip().lower()
        sp = sp.split("(")[0].strip()
        sp = re.sub(r"[^a-z0-9]+", "-", sp).strip("-")

        # Normalise molecule name
        mol = _normalize_ag_name(ag_mol).lower()
        mol = re.sub(r"[^a-z0-9;]+", "-", mol).strip("-")

        groups.append(f"{sp}/{mol}")

    df = df.copy()
    df["antigen_group"] = groups
    return df


def _split_by_antigen(
    df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Greedy antigen-based split targeting 70/15/15."""
    group_counts = df.groupby("antigen_group").size().sort_values(ascending=False).to_dict()

    targets = {"train": SPLIT_RATIOS[0], "eval": SPLIT_RATIOS[1], "test": SPLIT_RATIOS[2]}
    counts: dict[str, int] = {"train": 0, "eval": 0, "test": 0}
    assignments: dict[str, str] = {}

    for group, n in group_counts.items():
        total_after = sum(counts.values()) + n

        # Pick split that minimises total squared deviation from targets
        best_split = "train"
        best_err = float("inf")
        for split in ("train", "eval", "test"):
            trial = counts.copy()
            trial[split] += n
            err = sum((trial[s] / total_after - targets[s]) ** 2 for s in ("train", "eval", "test"))
            if err < best_err:
                best_err = err
                best_split = split

        assignments[group] = best_split
        counts[best_split] += n

    # Apply
    df = df.copy()
    df["_split"] = df["antigen_group"].map(assignments)
    train = df[df["_split"] == "train"].drop(columns=["_split"])
    eval_ = df[df["_split"] == "eval"].drop(columns=["_split"])
    test = df[df["_split"] == "test"].drop(columns=["_split"])

    ag_map = pd.DataFrame(
        [
            {"antigen_group": g, "split": s, "n_entries": group_counts[g]}
            for g, s in sorted(assignments.items())
        ]
    )

    for name, part in [("train", train), ("eval", eval_), ("test", test)]:
        pct = len(part) / len(df) * 100
        print(
            f"  {name:>5}: {len(part)} entries ({pct:.1f}%), "
            f"{part['antigen_group'].nunique()} antigen groups"
        )

    return train, eval_, test, ag_map


# ── Stage 7 ──────────────────────────────────────────────────────────────


def _print_summary(
    train: pd.DataFrame,
    eval_: pd.DataFrame,
    test: pd.DataFrame,
    splits_dir: Path,
) -> None:
    """Print a summary banner."""
    total = len(train) + len(eval_) + len(test)
    bar = "=" * 65

    print(f"\n{bar}")
    print("  SAAINT-DB Split Summary")
    print(bar)
    print(f"  Total entries: {total}")
    print()

    for name, part in [("Train", train), ("Eval", eval_), ("Test", test)]:
        pct = len(part) / total * 100
        kd_lo = part["neg_log10_KD"].min()
        kd_hi = part["neg_log10_KD"].max()
        kd_med = part["neg_log10_KD"].median()
        res_med = part["Resolution"].median()
        n_groups = part["antigen_group"].nunique()

        print(
            f"  {name:>5}: {len(part):>4} entries ({pct:5.1f}%)  "
            f"neg_log10_KD {kd_lo:+.2f} / {kd_med:+.2f} / {kd_hi:+.2f}  "
            f"res_med {res_med:.2f} A  "
            f"{n_groups} Ag groups"
        )

        ab_counts = Counter(part["Ab_type"])
        ab_str = ", ".join(f"{t}={c}" for t, c in ab_counts.most_common())
        print(f"         Ab types: {ab_str}")

        method_counts = Counter(part["Affinity_method"])
        m_str = ", ".join(f"{m}={c}" for m, c in method_counts.most_common())
        print(f"         Methods:  {m_str}")
        print()

    print(bar)
    print(f"  Splits written to {splits_dir}/")
    print(bar)


if __name__ == "__main__":
    main()
