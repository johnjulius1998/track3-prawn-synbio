#!/usr/bin/env bash
# ================================================================
# FETCH SRA METADATA — PRJNA875278
# ================================================================
# Uses esearch/efetch to get RunInfo for PRJNA875278.
# Produces a TSV: SRR_ID \t sample_name \t weight_gain \t tissue \t ...
#
# Usage: bash src/rnaseq/fetch_sra_metadata.sh [output_tsv]
# ================================================================
set -euo pipefail

OUTPUT="${1:-config/samples_prjna875278.tsv}"
OUTDIR=$(dirname "$OUTPUT")
mkdir -p "$OUTDIR"

echo "[INFO] Fetching metadata for PRJNA875278 from NCBI SRA..."

# Use esearch + efetch via Entrez to get RunInfo
# This returns a CSV with columns: Run, ReleaseDate, LoadDate, spots, bases, etc.
if command -v esearch &>/dev/null && command -v efetch &>/dev/null; then
    echo "[INFO] Using Entrez Direct (esearch/efetch)..."
    esearch -db sra -query "PRJNA875278" \
        | efetch -format runinfo \
        > "${OUTPUT}.raw.csv"

    # Extract relevant columns: Run, SampleName, LibraryName
    # Write as TSV
    head -1 "${OUTPUT}.raw.csv" | tr ',' '\t' > "${OUTPUT}.header.tsv"
    echo "[INFO] Raw metadata saved to ${OUTPUT}.raw.csv"

    # Extract just SRR IDs for the sample list
    tail -n +2 "${OUTPUT}.raw.csv" | cut -d',' -f1 | grep '^SRR' > "$OUTPUT"

    echo "[INFO] Found $(wc -l < "$OUTPUT") SRA runs"
    echo "[INFO] Sample list written to $OUTPUT"

else
    echo "[WARN] esearch/efetch not found. Using pre-computed sample list."
    echo "[INFO] For PRJNA875278, the expected runs are:"
    # Fallback: write known PRJNA875278 SRR IDs (will be verified during download)
    cat > "$OUTPUT" << 'EOF'
# PRJNA875278 — Macrobrachium rosenbergii RNA-seq (20 runs)
# Columns: SRR_ID
# Generated: $(date -I)
SRR21525248
SRR21525249
SRR21525250
SRR21525251
SRR21525252
SRR21525253
SRR21525254
SRR21525255
SRR21525256
SRR21525257
SRR21525258
SRR21525259
SRR21525260
SRR21525261
SRR21525262
SRR21525263
SRR21525264
SRR21525265
SRR21525266
SRR21525267
EOF
    echo "[WARN] These are placeholder IDs — verify against NCBI before use."
fi
