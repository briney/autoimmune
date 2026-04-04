"""One-time data splitting for CR9114 PoC.

Produces stratified train/eval/test splits from the full H1 binding dataset.
Run once:  python prepare.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

SEED = 42
N_TRAIN = 500
N_EVAL = 100
N_TEST = 100


def main() -> None:
    root = Path(__file__).parent
    splits_dir = root / "splits"
    splits_dir.mkdir(exist_ok=True)

    # Load and clean
    df = pd.read_csv(root / "data" / "cr9114_h1_binding_data.csv", dtype={"genotype": str})
    n_raw = len(df)
    df = df.dropna(subset=["h1_mean"])
    print(f"Loaded {n_raw} variants, {len(df)} with H1 binding data")

    rng = np.random.default_rng(SEED)

    # Stratify by affinity quartiles so each split covers the full dynamic range
    df["quartile"] = pd.qcut(df["h1_mean"], q=4, labels=False)

    train_idx: list[int] = []
    eval_idx: list[int] = []
    test_idx: list[int] = []

    per_q_train = N_TRAIN // 4  # 125
    per_q_eval = N_EVAL // 4  # 25
    per_q_test = N_TEST // 4  # 25

    for q in range(4):
        pool = df.index[df["quartile"] == q].to_numpy().copy()
        rng.shuffle(pool)
        a = per_q_train
        b = a + per_q_eval
        c = b + per_q_test
        train_idx.extend(pool[:a])
        eval_idx.extend(pool[a:b])
        test_idx.extend(pool[b:c])

    train = df.loc[train_idx, ["genotype", "h1_mean"]].reset_index(drop=True)
    eval_ = df.loc[eval_idx, ["genotype", "h1_mean"]].reset_index(drop=True)
    test = df.loc[test_idx, ["genotype", "h1_mean"]].reset_index(drop=True)

    # Sanity: no genotype appears in multiple splits
    all_genos = set(train["genotype"]) | set(eval_["genotype"]) | set(test["genotype"])
    assert len(all_genos) == N_TRAIN + N_EVAL + N_TEST, "Duplicate genotypes across splits"

    # Verify all 16 positions appear in both states in each split
    for name, part in [("train", train), ("eval", eval_), ("test", test)]:
        for pos in range(16):
            states = {g[pos] for g in part["genotype"]}
            if states != {"0", "1"}:
                missing = "germline (0)" if "0" not in states else "somatic (1)"
                print(f"  WARNING: {name} missing {missing} at position {pos + 1}")

    # Find cross-boundary Hamming-1 pairs for pairwise accuracy
    eval_pairs = _find_cross_pairs(train, eval_)
    test_pairs = _find_cross_pairs(train, test)

    # Write splits
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

    # Summary
    print(f"\n{'=' * 55}")
    print("  Split summary")
    print(f"{'=' * 55}")
    for name, part in [("Train", train), ("Eval", eval_), ("Test", test)]:
        lo, hi = part["h1_mean"].min(), part["h1_mean"].max()
        print(f"  {name:>5}: {len(part):>4} variants  (h1_mean {lo:.2f} – {hi:.2f})")
    print(f"  Eval cross-boundary Hamming-1 pairs: {len(eval_pairs)}")
    print(f"  Test cross-boundary Hamming-1 pairs: {len(test_pairs)}")
    print(f"{'=' * 55}")
    print(f"\nSplits written to {splits_dir}/")


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
