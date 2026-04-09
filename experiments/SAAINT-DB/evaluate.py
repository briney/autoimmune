"""Evaluate agent predictions against held-out ground truth.

Usage:
    python evaluate.py predictions.csv [--partition eval|test]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import pearsonr, spearmanr


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Evaluate binding affinity predictions against held-out truth",
    )
    parser.add_argument(
        "predictions",
        type=Path,
        help="CSV with complex_id,predicted_score columns",
    )
    parser.add_argument(
        "--partition",
        choices=["eval", "test"],
        default="eval",
        help="Which partition to evaluate against (default: eval)",
    )
    args = parser.parse_args()

    splits_dir = Path(__file__).parent / "splits"

    # Load predictions
    preds = pd.read_csv(args.predictions)
    if "complex_id" not in preds.columns or "predicted_score" not in preds.columns:
        print(
            "ERROR: predictions CSV must have 'complex_id' and 'predicted_score' columns",
            file=sys.stderr,
        )
        sys.exit(1)
    preds = preds.drop_duplicates(subset="complex_id", keep="first")

    # Load truth for the requested partition
    truth = pd.read_csv(splits_dir / f"{args.partition}_truth.csv")
    merged = truth.merge(
        preds[["complex_id", "predicted_score"]], on="complex_id", how="inner"
    )

    n_total = len(truth)
    n_matched = len(merged)
    coverage = n_matched / n_total

    if n_matched < 2:
        print(
            f"ERROR: only {n_matched} predictions matched {args.partition} complexes",
            file=sys.stderr,
        )
        sys.exit(1)
    if n_matched < n_total:
        print(
            f"WARNING: predictions cover {n_matched}/{n_total} {args.partition} complexes",
            file=sys.stderr,
        )

    # -- Spearman rho --
    rho, rho_pval = spearmanr(merged["predicted_score"], merged["neg_log10_KD"])

    # -- Pearson r --
    r, r_pval = pearsonr(merged["predicted_score"], merged["neg_log10_KD"])

    # -- Top-k precision (k = 10% of matched, minimum 1) --
    k = max(1, n_matched // 10)
    top_pred = set(merged.nlargest(k, "predicted_score")["complex_id"])
    top_true = set(merged.nlargest(k, "neg_log10_KD")["complex_id"])
    top_k_prec = len(top_pred & top_true) / k

    # -- Output --
    print(f"\n{'=' * 56}")
    print(f"  Evaluation ({args.partition}, {n_matched}/{n_total} complexes)")
    print(f"{'=' * 56}")
    print(f"  Spearman rho:       {rho:+.4f}  (p = {rho_pval:.2e})")
    print(f"  Pearson r:          {r:+.4f}  (p = {r_pval:.2e})")
    print(f"  Top-{k} precision:   {top_k_prec:.2f}")
    print(f"  Coverage:           {coverage:.2f}  ({n_matched}/{n_total})")
    print(f"{'=' * 56}\n")


if __name__ == "__main__":
    main()
