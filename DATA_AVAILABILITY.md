# Data Availability

This repository is a curated subset of the full analysis project. It includes all code, config, small metadata/sample-sheet files, final processed outputs, and result figures/tables — but excludes bulk raw sequencing data and large regenerable reference downloads, to keep the repository small and avoid redistributing data that is already public elsewhere.

## Included

- `data/raw/supplied/ASV_table_Jumpers_Laggards.GKAQUA.csv` — the original supplied 16S ASV abundance table (2 pooled samples). This is the only primary input not available from a public database, so it is bundled directly.
- `data/raw/sra/PRJNA875278/{metadata.tsv,samples.tsv,samples_corrected.tsv,runinfo_cleaned.csv}` — sample sheets and run metadata for the RNA-seq samples (sample IDs, sex, tissue, run accessions). These are small text files, not sequencing reads.
- `data/processed/` — final CLR-transformed microbial profiles, merged host gene-expression matrices, and WGCNA module assignments. These are the direct inputs to every downstream figure/table in `results/`.
- `results/` — all final figures, tables, validation reports, and the official ranked shortlists.

## Excluded (and how to regenerate)

| Excluded data | Size | Source / regeneration |
|---|---|---|
| Raw RNA-seq reads (`*.fastq.gz`, both raw and trimmed, for 20 samples) | ~156 GB | Public on NCBI SRA under BioProject **[PRJNA875278](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA875278)**. Fetch with `pipeline/src/rnaseq/fetch_sra_metadata.sh` and `pipeline/src/rnaseq/download_sra_safe.sh`, using the run accessions in `data/raw/sra/PRJNA875278/metadata.tsv`. |
| Reference transcriptome/proteome and salmon index (`M_rosenbergii_*.fasta`, `salmon_index/`) | ~2.9 GB | *M. rosenbergii* TSA + mRNA records from NCBI, combined and indexed by `pipeline/src/network/fetch_prawn_proteins.py`. |
| UniProt Swiss-Prot reference (`uniprot_sprot.fasta`, `.dmnd`) | ~580 MB | Downloaded and DIAMOND-indexed by `pipeline/src/network/map_uniprot_sets.py`. |
| Intermediate salmon per-sample quant outputs, literature/API caches (`data/interim/`) | ~166 MB | Regenerated automatically by the corresponding Snakemake rules once raw reads are present; not needed to reproduce any figure or table, since `data/processed/` already contains their final, merged form. |

## Reproducing the RNA-seq layer from scratch

1. Fetch raw reads for BioProject PRJNA875278 (see table above).
2. Place them under `data/raw/sra/PRJNA875278/` following the existing sample-sheet naming.
3. Run `snakemake -s pipeline/Snakefile --profile pipeline/profiles/local` — the workflow will fetch/build the reference and Swiss-Prot databases, quantify with salmon, and merge counts.

Note: the results shipped in this repository were originally produced via `pipeline/src/rnaseq/process_all_samples.sh` run as a serial driver over the raw downloads (documented in `pipeline/Snakefile`), not a single `snakemake` invocation end-to-end — minor nondeterminism in tool versions/thread scheduling could cause small numerical differences from a from-scratch rerun.

## Credentials

`pipeline/src/asv/03c_metacyc_mapping.py` requires BioCyc account credentials, read from the environment variables `BIOCYC_USER` and `BIOCYC_PASSWORD`. No credentials are stored in this repository; register for a free BioCyc account to run this step.
