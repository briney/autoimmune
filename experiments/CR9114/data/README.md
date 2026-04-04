# CR9114 / H1 Subset

This directory contains a reduced CR9114/H1-specific subset derived from the
Phillips et al. influenza antibody binding landscape dataset:

- Phillips AM et al. "Binding affinity landscapes constrain the evolution of
  broadly neutralizing anti-influenza antibodies" eLife 2021, PMID `34491198`,
  DOI `10.7554/eLife.71393`

The original deposited assets and the broader reconstruction workflow are stored
in `data/manuscript_binding_landscapes/`. Those files were built from:

- eLife source data CSVs
- eLife supplementary plasmid/sequence files
- the companion GitHub repository `klawrence26/bnab-landscapes`
- raw sequencing deposition reported by the paper: NCBI BioProject `PRJNA741613`

Files in this directory:

- `cr9114_h1_binding_data.csv`
  - Reduced H1-only table derived from the CR9114 source-data CSV.
  - Contains `genotype`, the three H1 replicate measurements
    (`h1_repa`, `h1_repb`, `h1_repc`), and `h1_mean`.
  - `h1_mean` is carried through exactly as reported in the source data
    (`-log10(KD)`), not converted to molar KD.

- `cr9114_mutation_key.csv`
  - CR9114-only heavy-chain mutation key.
  - Maps each binary genotype position to the corresponding heavy-chain residue
    position and germline/somatic amino-acid pair.

- `cr9114_h1_sequences.fasta`
  - FASTA file containing the mature CR9114 heavy chain, the CR9114 light
    chain, and the H1 hemagglutinin amino-acid sequence.
  - The H1 sequence is trimmed before the common Avi/His-tag tail, matching the
    trimming described in `data/manuscript_binding_landscapes/README.md`.
  - The deposited N-terminal secretion leader is retained.
