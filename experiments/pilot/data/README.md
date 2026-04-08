# Pilot Experiment Data

This directory contains binding affinity data for the pilot antibody-antigen system
used to benchmark flex-ddG-like structural scoring pipelines.

Files in this directory:

- `binding_data.csv`
  - Full combinatorial binding dataset for all antibody heavy-chain variants.
  - Contains `genotype`, three replicate measurements
    (`affinity_repa`, `affinity_repb`, `affinity_repc`), and `affinity_mean`.
  - `affinity_mean` is reported as −log₁₀(KD); higher values indicate tighter binding.

- `mutation_key.csv`
  - Heavy-chain mutation key.
  - Maps each binary genotype position to the corresponding heavy-chain residue
    position and germline/somatic amino-acid pair.

- `sequences.fasta`
  - FASTA file containing the mature antibody heavy chain, the antibody light
    chain, and the antigen amino-acid sequence.
