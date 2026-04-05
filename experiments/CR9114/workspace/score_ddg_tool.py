from __future__ import annotations

import argparse
import json
import logging
import subprocess
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


LOGGER = logging.getLogger(__name__)
SUPPORTED_TOOLS = {"stabddg", "baddg"}


@dataclass(slots=True)
class MutationSite:
    """Single reversible mutation site from the CR9114 mutation key."""

    bit: int
    position: int
    germline_aa: str
    somatic_aa: str


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Score CR9114 genotypes with a direct ddG autobio tool "
            "(currently stabddg or baddg)."
        )
    )
    parser.add_argument(
        "--tool",
        choices=sorted(SUPPORTED_TOOLS),
        required=True,
        help="Which autobio ddG tool to run.",
    )
    parser.add_argument(
        "--input-csv",
        type=Path,
        required=True,
        help="CSV containing at least a genotype column.",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        required=True,
        help="Where to write the summarized ddG scores.",
    )
    parser.add_argument(
        "--experiment-root",
        type=Path,
        default=Path(__file__).resolve().parent.parent,
        help="Path to the CR9114 experiment directory.",
    )
    parser.add_argument(
        "--interface-chains",
        default="ABC_DE",
        help="Interface chain specification passed to the ddG tool.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Random seed passed through to the tool when supported.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional row limit for quick calibrations.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Reuse existing per-genotype JSON outputs when present.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity.",
    )
    return parser.parse_args()


def configure_logging(log_level: str) -> None:
    """Configure logging."""
    logging.basicConfig(level=getattr(logging, log_level), format="%(levelname)s %(message)s")


def load_mutation_key(key_path: Path) -> list[MutationSite]:
    """Load the CR9114 mutation key into a positional lookup list."""
    key = pd.read_csv(key_path).rename(
        columns={
            "binary_position": "bit",
            "heavy_chain_sequence_position_1_based": "position",
        }
    )
    sites = [
        MutationSite(
            bit=int(row["bit"]),
            position=int(row["position"]),
            germline_aa=str(row["germline_aa"]),
            somatic_aa=str(row["somatic_aa"]),
        )
        for _, row in key.sort_values("bit").iterrows()
    ]
    if len(sites) != 16:
        raise ValueError(f"Expected 16 mutation sites, found {len(sites)}")
    return sites


def genotype_to_mutations(genotype: str, sites: list[MutationSite]) -> list[str]:
    """Convert a CR9114 genotype into a list of mature-to-germline reversions."""
    normalized = genotype.strip().zfill(16)
    mutations: list[str] = []
    for bit_char, site in zip(normalized, sites, strict=True):
        if bit_char == "0":
            mutations.append(f"{site.somatic_aa}D{site.position}{site.germline_aa}")
    return mutations


def score_genotype(
    *,
    tool: str,
    genotype: str,
    mutations: list[str],
    structure_path: Path,
    interface_chains: str,
    seed: int,
    run_root: Path,
    experiment_root: Path,
    resume: bool,
) -> tuple[float, dict]:
    """Run one ddG tool invocation and return its parsed score."""
    run_dir = run_root / genotype
    score_path = run_dir / "score.json"
    if resume and score_path.exists():
        payload = json.loads(score_path.read_text(encoding="utf-8"))
    else:
        run_dir.mkdir(parents=True, exist_ok=True)
        config = {
            "structure_path": str(structure_path),
            "extra": {
                "mutations": mutations,
                "chains": interface_chains,
                "seed": seed,
            },
        }
        config_path = run_dir / "config.json"
        config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
        command = [
            "autobio",
            "run",
            tool,
            "--config",
            str(config_path),
            "--output-dir",
            str(run_dir / "workspace"),
            "--format",
            "json",
        ]
        completed = subprocess.run(
            command,
            check=True,
            cwd=experiment_root,
            capture_output=True,
            text=True,
        )
        payload = json.loads(completed.stdout)
        score_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    score = float(payload["scores"][0]["ddg"])
    return score, payload


def main() -> None:
    """Score a genotype CSV with the selected ddG tool."""
    args = parse_args()
    configure_logging(args.log_level)

    experiment_root = args.experiment_root.resolve()
    structure_path = experiment_root / "structures/cr9114_mature_h1.pdb"
    sites = load_mutation_key(experiment_root / "data/cr9114_mutation_key.csv")

    frame = pd.read_csv(args.input_csv, dtype={"genotype": str})
    frame["genotype"] = frame["genotype"].str.zfill(16)
    if args.limit is not None:
        frame = frame.head(args.limit).copy()

    run_root = args.output_csv.resolve().parent / f"{args.tool}_runs"
    run_root.mkdir(parents=True, exist_ok=True)

    results: list[dict] = []
    total = len(frame)
    for index, row in enumerate(frame.itertuples(index=False), start=1):
        genotype = str(row.genotype).zfill(16)
        mutations = genotype_to_mutations(genotype, sites)
        LOGGER.info("[%d/%d] scoring %s with %d reversions", index, total, genotype, len(mutations))
        ddg, payload = score_genotype(
            tool=args.tool,
            genotype=genotype,
            mutations=mutations,
            structure_path=structure_path,
            interface_chains=args.interface_chains,
            seed=args.seed,
            run_root=run_root,
            experiment_root=experiment_root,
            resume=args.resume,
        )
        result = {
            "genotype": genotype,
            "predicted_ddg": ddg,
            "predicted_score": -ddg,
            "n_reversions": len(mutations),
            "mutations": ";".join(mutations),
            "wall_time_seconds": payload["metadata"]["wall_time_seconds"],
        }
        if hasattr(row, "h1_mean"):
            result["h1_mean"] = float(row.h1_mean)
        results.append(result)

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(results).to_csv(args.output_csv, index=False)
    LOGGER.info("Wrote %d scored rows to %s", len(results), args.output_csv)


if __name__ == "__main__":
    main()
