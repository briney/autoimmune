"""Evaluate agent predictions against held-out ground truth.

Usage:
    python evaluate.py predictions.csv [--partition eval|test]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Evaluate binding affinity predictions against held-out truth",
    )
    parser.add_argument("predictions", type=Path, help="CSV with genotype,predicted_score columns")
    parser.add_argument(
        "--partition",
        choices=["eval", "test"],
        default="eval",
        help="Which partition to evaluate against (default: eval)",
    )
    args = parser.parse_args()

    splits_dir = Path(__file__).parent / "splits"

    # Load predictions (force genotype to string to handle binary-like values)
    preds = pd.read_csv(args.predictions, dtype={"genotype": str})
    if "genotype" not in preds.columns or "predicted_score" not in preds.columns:
        print(
            "ERROR: predictions CSV must have 'genotype' and 'predicted_score' columns",
            file=sys.stderr,
        )
        sys.exit(1)
    preds["genotype"] = preds["genotype"].str.zfill(16)
    preds = preds.drop_duplicates(subset="genotype", keep="first")

    # Load truth for the requested partition
    truth = pd.read_csv(splits_dir / f"{args.partition}_truth.csv", dtype={"genotype": str})
    merged = truth.merge(preds[["genotype", "predicted_score"]], on="genotype", how="inner")

    n_total = len(truth)
    n_matched = len(merged)
    if n_matched < 2:
        print(
            f"ERROR: only {n_matched} predictions matched {args.partition} genotypes",
            file=sys.stderr,
        )
        sys.exit(1)
    if n_matched < n_total:
        print(
            f"WARNING: predictions cover {n_matched}/{n_total} {args.partition} variants",
            file=sys.stderr,
        )

    # -- Spearman rho --
    rho, pval = spearmanr(merged["predicted_score"], merged["h1_mean"])

    # -- Top-k precision (k = 10% of partition, minimum 1) --
    k = max(1, n_matched // 10)
    top_pred = set(merged.nlargest(k, "predicted_score")["genotype"])
    top_true = set(merged.nlargest(k, "h1_mean")["genotype"])
    top_k_prec = len(top_pred & top_true) / k

    # -- Pairwise accuracy on cross-boundary Hamming-1 pairs --
    pairs_path = splits_dir / f"{args.partition}_pairs.csv"
    pair_str = "N/A (no pairs file)"
    if pairs_path.exists():
        pairs = pd.read_csv(pairs_path, dtype=str)
        # Truth comes from train labels + partition labels
        train_truth = pd.read_csv(splits_dir / "train.csv", dtype={"genotype": str})
        all_truth = (
            pd.concat([train_truth, truth], ignore_index=True)
            .drop_duplicates(subset="genotype")
            .set_index("genotype")["h1_mean"]
        )
        pred_scores = preds.set_index("genotype")["predicted_score"]

        correct, total = 0, 0
        for _, row in pairs.iterrows():
            ga, gb = row.iloc[0], row.iloc[1]
            has_preds = ga in pred_scores.index and gb in pred_scores.index
            has_truth = ga in all_truth.index and gb in all_truth.index
            if has_preds and has_truth:
                if (pred_scores[ga] > pred_scores[gb]) == (all_truth[ga] > all_truth[gb]):
                    correct += 1
                total += 1

        if total > 0:
            pair_str = f"{correct / total:.2f} ({correct}/{total} pairs)"
        else:
            pair_str = "N/A (no pairs with predictions for both members)"

    # -- Output --
    print(f"\n{'=' * 50}")
    print(f"  Evaluation ({args.partition}, {n_matched} variants)")
    print(f"{'=' * 50}")
    print(f"  Spearman rho:       {rho:+.4f}  (p = {pval:.2e})")
    print(f"  Top-{k} precision:   {top_k_prec:.2f}")
    print(f"  Pairwise accuracy:  {pair_str}")
    print(f"{'=' * 50}\n")


if __name__ == "__main__":
    main()
