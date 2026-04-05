"""One-time data splitting for CR9114 PoC.

Produces stratified train/eval/test splits from the full H1 binding dataset.
Restricts to a narrow mutation-count window (6-8 mutations) to prevent the
agent from using mutation count as a proxy for affinity.

Run once:  python prepare.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

SEED = 42
N_TRAIN = 500
N_EVAL = 100
N_TEST = 100

# Mutation-count window: restrict to variants with 6-8 somatic mutations.
# This reduces Spearman(n_mut, h1_mean) from ~0.37 (full dataset) to ~0.14,
# forcing the agent to rely on sequence/structure signal rather than mutation count.
MUT_COUNT_MIN = 6
MUT_COUNT_MAX = 8

# High-impact pairs: Hamming-1 neighbors with >=10-fold KD difference
HIGH_IMPACT_DELTA = 1.0
N_HIGH_IMPACT_EVAL = 15
N_HIGH_IMPACT_TEST = 15

# Fraction of high-impact pairs that should be "reverse" (more mutations but
# LOWER affinity) — breaks the mutation-count heuristic
REVERSE_FRACTION = 0.25


def main() -> None:
    root = Path(__file__).parent
    splits_dir = root / "splits"
    splits_dir.mkdir(exist_ok=True)

    # ── Load and filter ──────────────────────────────────────────────────

    df = pd.read_csv(root / "data" / "cr9114_h1_binding_data.csv", dtype={"genotype": str})
    n_raw = len(df)
    df = df.dropna(subset=["h1_mean"])
    df["n_mut"] = df["genotype"].apply(lambda g: g.count("1"))

    pool_df = df[(df["n_mut"] >= MUT_COUNT_MIN) & (df["n_mut"] <= MUT_COUNT_MAX)].copy()
    rho_full, _ = spearmanr(df["n_mut"], df["h1_mean"])
    rho_window, _ = spearmanr(pool_df["n_mut"], pool_df["h1_mean"])
    print(f"Loaded {n_raw} variants, {len(df)} with H1 data")
    print(f"Mutation window {MUT_COUNT_MIN}-{MUT_COUNT_MAX}: {len(pool_df)} variants")
    print(f"  Spearman(n_mut, h1_mean): full={rho_full:.3f}, window={rho_window:.3f}")

    rng = np.random.default_rng(SEED)

    # ── Phase 1: stratified random sampling ──────────────────────────────

    pool_df["quartile"] = pd.qcut(pool_df["h1_mean"], q=4, labels=False)

    train_idx: list[int] = []
    eval_idx: list[int] = []
    test_idx: list[int] = []

    per_q_train = N_TRAIN // 4  # 125
    per_q_eval = N_EVAL // 4  # 25
    per_q_test = N_TEST // 4  # 25

    for q in range(4):
        idxs = pool_df.index[pool_df["quartile"] == q].to_numpy().copy()
        rng.shuffle(idxs)
        a = per_q_train
        b = a + per_q_eval
        c = b + per_q_test
        train_idx.extend(idxs[:a])
        eval_idx.extend(idxs[a:b])
        test_idx.extend(idxs[b:c])

    train = pool_df.loc[train_idx, ["genotype", "h1_mean"]].reset_index(drop=True)
    eval_ = pool_df.loc[eval_idx, ["genotype", "h1_mean"]].reset_index(drop=True)
    test = pool_df.loc[test_idx, ["genotype", "h1_mean"]].reset_index(drop=True)

    all_genos = set(train["genotype"]) | set(eval_["genotype"]) | set(test["genotype"])
    assert len(all_genos) == N_TRAIN + N_EVAL + N_TEST, "Duplicate genotypes across splits"

    # ── Phase 2: high-impact pair injection ──────────────────────────────

    high_impact = _find_high_impact_pairs(pool_df, min_delta=HIGH_IMPACT_DELTA)
    normal = [p for p in high_impact if not p[4]]
    reverse = [p for p in high_impact if p[4]]
    print(f"High-impact pairs in window: {len(high_impact)} total "
          f"({len(normal)} normal, {len(reverse)} reverse)")

    # Filter to pairs where neither member is already in a split
    def available(pairs: list[_Pair]) -> list[_Pair]:
        return [p for p in pairs if p[0] not in all_genos and p[1] not in all_genos]

    normal_avail = available(normal)
    reverse_avail = available(reverse)
    rng.shuffle(normal_avail)
    rng.shuffle(reverse_avail)

    # Select for eval: fill reverse quota first, then normal
    n_reverse_eval = round(N_HIGH_IMPACT_EVAL * REVERSE_FRACTION)
    eval_reverse = _select_diverse_pairs(reverse_avail, n_reverse_eval)
    used = _genotypes_in(eval_reverse)
    eval_normal = _select_diverse_pairs(
        [p for p in normal_avail if p[0] not in used and p[1] not in used],
        N_HIGH_IMPACT_EVAL - len(eval_reverse),
    )
    eval_extra = eval_reverse + eval_normal
    used |= _genotypes_in(eval_normal)

    # Select for test: same approach
    n_reverse_test = round(N_HIGH_IMPACT_TEST * REVERSE_FRACTION)
    test_reverse = _select_diverse_pairs(
        [p for p in reverse_avail if p[0] not in used and p[1] not in used],
        n_reverse_test,
    )
    used |= _genotypes_in(test_reverse)
    test_normal = _select_diverse_pairs(
        [p for p in normal_avail if p[0] not in used and p[1] not in used],
        N_HIGH_IMPACT_TEST - len(test_reverse),
    )
    test_extra = test_reverse + test_normal

    # Inject pairs into splits
    train, eval_ = _inject_pairs(train, eval_, eval_extra, pool_df, rng)
    train, test = _inject_pairs(train, test, test_extra, pool_df, rng)

    n_eval_rev = len(eval_reverse)
    n_test_rev = len(test_reverse)
    print(f"Injected {len(eval_extra)} pairs into train/eval ({n_eval_rev} reverse)")
    print(f"Injected {len(test_extra)} pairs into train/test ({n_test_rev} reverse)")

    # Re-check uniqueness
    all_genos = set(train["genotype"]) | set(eval_["genotype"]) | set(test["genotype"])
    assert len(all_genos) == len(train) + len(eval_) + len(test), (
        "Duplicate genotypes across splits after injection"
    )

    # ── Verification ─────────────────────────────────────────────────────

    for name, part in [("train", train), ("eval", eval_), ("test", test)]:
        for pos in range(16):
            states = {g[pos] for g in part["genotype"]}
            if states != {"0", "1"}:
                missing = "germline (0)" if "0" not in states else "somatic (1)"
                print(f"  WARNING: {name} missing {missing} at position {pos + 1}")

    eval_pairs = _find_cross_pairs(train, eval_)
    test_pairs = _find_cross_pairs(train, test)

    # ── Write outputs ────────────────────────────────────────────────────

    train.to_csv(splits_dir / "train.csv", index=False)
    eval_[["genotype"]].to_csv(splits_dir / "eval_genotypes.csv", index=False)
    eval_.to_csv(splits_dir / "eval_truth.csv", index=False)
    test[["genotype"]].to_csv(splits_dir / "test_genotypes.csv", index=False)
    test.to_csv(splits_dir / "test_truth.csv", index=False)
    pd.DataFrame(eval_pairs, columns=["genotype_a", "genotype_b"]).to_csv(
        splits_dir / "eval_pairs.csv", index=False
    )
    pd.DataFrame(test_pairs, columns=["genotype_a", "genotype_b"]).to_csv(
        splits_dir / "test_pairs.csv", index=False
    )

    # ── Summary ──────────────────────────────────────────────────────────

    print(f"\n{'=' * 65}")
    print("  Split summary")
    print(f"{'=' * 65}")
    for name, part in [("Train", train), ("Eval", eval_), ("Test", test)]:
        part_mut = part["genotype"].apply(lambda g: g.count("1"))
        rho_part, _ = spearmanr(part_mut, part["h1_mean"])
        lo, hi = part["h1_mean"].min(), part["h1_mean"].max()
        print(f"  {name:>5}: {len(part):>4} variants  "
              f"(h1_mean {lo:.2f}–{hi:.2f}, rho(n_mut,h1)={rho_part:+.3f})")
    print(f"  Eval cross-boundary Hamming-1 pairs: {len(eval_pairs)}")
    print(f"  Test cross-boundary Hamming-1 pairs: {len(test_pairs)}")
    print(f"{'=' * 65}")
    print(f"\nSplits written to {splits_dir}/")


# ── Type alias for pair tuples ───────────────────────────────────────────

# (genotype_a, genotype_b, bit_position, delta, is_reverse)
type _Pair = tuple[str, str, int, float, bool]


# ── Helper functions ─────────────────────────────────────────────────────


def _find_high_impact_pairs(df: pd.DataFrame, min_delta: float) -> list[_Pair]:
    """Find all Hamming-1 pairs in df with large affinity difference.

    Each pair is tagged as "reverse" if the more-mutated member has lower
    affinity (breaks the mutation-count → affinity heuristic).
    """
    h1 = df.set_index("genotype")["h1_mean"]
    nmut = df.set_index("genotype")["n_mut"]
    geno_set = set(h1.index)
    pairs: list[_Pair] = []
    for g in geno_set:
        for i in range(16):
            flipped = g[:i] + ("1" if g[i] == "0" else "0") + g[i + 1:]
            if flipped in geno_set and g < flipped:
                delta = abs(h1[g] - h1[flipped])
                if delta >= min_delta:
                    # Determine if "reverse": more mutations but worse binding
                    if nmut[g] != nmut[flipped]:
                        more = g if nmut[g] > nmut[flipped] else flipped
                        is_reverse = h1[more] < h1[g if more == flipped else flipped]
                    else:
                        is_reverse = False  # same mutation count — not applicable
                    pairs.append((g, flipped, i + 1, delta, is_reverse))
    pairs.sort(key=lambda p: -p[3])
    return pairs


def _select_diverse_pairs(pairs: list[_Pair], n: int) -> list[_Pair]:
    """Select up to n pairs, preferring diversity across mutation positions."""
    selected: list[_Pair] = []
    used_positions: set[int] = set()
    used_genos: set[str] = set()

    # First pass: one pair per position
    for p in pairs:
        if len(selected) >= n:
            break
        if p[2] not in used_positions and p[0] not in used_genos and p[1] not in used_genos:
            selected.append(p)
            used_positions.add(p[2])
            used_genos.update((p[0], p[1]))

    # Second pass: fill remaining slots regardless of position
    for p in pairs:
        if len(selected) >= n:
            break
        if p[0] not in used_genos and p[1] not in used_genos:
            selected.append(p)
            used_genos.update((p[0], p[1]))

    return selected


def _genotypes_in(pairs: list[_Pair]) -> set[str]:
    """Collect all genotypes appearing in a list of pairs."""
    return {g for p in pairs for g in (p[0], p[1])}


def _inject_pairs(
    train: pd.DataFrame,
    other: pd.DataFrame,
    pairs: list[_Pair],
    source_df: pd.DataFrame,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Inject pair members: one into train, the other into other partition."""
    for ga, gb, *_ in pairs:
        row_a = source_df.loc[source_df["genotype"] == ga, ["genotype", "h1_mean"]].iloc[0]
        row_b = source_df.loc[source_df["genotype"] == gb, ["genotype", "h1_mean"]].iloc[0]
        if rng.random() < 0.5:
            train = pd.concat([train, row_a.to_frame().T], ignore_index=True)
            other = pd.concat([other, row_b.to_frame().T], ignore_index=True)
        else:
            train = pd.concat([train, row_b.to_frame().T], ignore_index=True)
            other = pd.concat([other, row_a.to_frame().T], ignore_index=True)
    return train, other


def _find_cross_pairs(
    train_df: pd.DataFrame, other_df: pd.DataFrame
) -> list[tuple[str, str]]:
    """Find all Hamming-1 neighbor pairs spanning train and the other partition."""
    train_set = set(train_df["genotype"])
    pairs: list[tuple[str, str]] = []
    for g in other_df["genotype"]:
        for i in range(16):
            flipped = g[:i] + ("1" if g[i] == "0" else "0") + g[i + 1:]
            if flipped in train_set:
                pairs.append((flipped, g))
    return pairs


if __name__ == "__main__":
    main()
