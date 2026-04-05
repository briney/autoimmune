"""
Build ESM2 input: construct variant heavy chain sequences from genotype + mutation key.
Saves a JSON config for autobio esm2 run.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

MATURE_HC = (
    "QVQLVQSGAEVKKPGSSVKVSCKASGGTSNNYAISWVRQAPGQGLEWMGGISPIFGSTAYAQKFQGRVTIS"
    "ADIFSNTAYMELNSLTSEDTAVYFCARHGNYYYYSGMDVWGQGTTVTVSS"
)

# Mutation key: 0-indexed position, germline aa
mutation_key = pd.read_csv(ROOT / "data/cr9114_mutation_key.csv")
mutations = {
    row.binary_position: {
        "pos0": row.heavy_chain_sequence_position_1_based - 1,  # 0-indexed
        "germline": row.germline_aa,
        "somatic": row.somatic_aa,
    }
    for _, row in mutation_key.iterrows()
}


def genotype_to_sequence(genotype: str) -> str:
    seq = list(MATURE_HC)
    for bit_idx, aa_str in enumerate(genotype):
        bit = bit_idx + 1  # 1-indexed
        if aa_str == "0":  # revert to germline
            m = mutations[bit]
            assert seq[m["pos0"]] == m["somatic"], (
                f"Bit {bit}: expected {m['somatic']} at pos {m['pos0']+1}, got {seq[m['pos0']]}"
            )
            seq[m["pos0"]] = m["germline"]
    return "".join(seq)


def main() -> None:
    train = pd.read_csv(ROOT / "splits/train.csv", dtype={"genotype": str})
    eval_g = pd.read_csv(ROOT / "splits/eval_genotypes.csv", dtype={"genotype": str})

    all_genotypes = pd.concat([
        train[["genotype"]],
        eval_g[["genotype"]],
    ]).drop_duplicates()

    # Verify mature sequence is correct
    mature_seq = genotype_to_sequence("1" * 16)
    assert mature_seq == MATURE_HC, "Mature reconstruction mismatch!"
    print(f"Mature sequence verified ({len(MATURE_HC)} residues)")

    sequences = {}
    for g in all_genotypes.genotype:
        sequences[g] = genotype_to_sequence(g)

    # Add mature as reference
    sequences["mature"] = MATURE_HC

    print(f"Total sequences: {len(sequences)}")

    config = {
        "sequences": sequences,
        "pooling": "mean",
    }

    out_path = ROOT / "workspace/esm2_config.json"
    with open(out_path, "w") as f:
        json.dump(config, f, indent=2)
    print(f"Written: {out_path}")


if __name__ == "__main__":
    main()
