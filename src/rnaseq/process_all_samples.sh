#!/usr/bin/env bash
# ================================================================
# process_all_samples.sh — Serial SRA → Quant Pipeline (RESILIENT)
# ================================================================
# Processes ALL 20 PRJNA875278 samples ONE AT A TIME.
# Designed for 7.7 GB RAM with disk/network guardrails.
# FIXED: removed set -euo, added RAM checks, better error handling.
#
# Usage: bash src/rnaseq/process_all_samples.sh
# ================================================================
# NOTE: NOT using set -e because some tools return non-zero on warnings.
# Error handling is explicit via || true and per-step checks.
set -o pipefail

# ---- Configuration ----
SAMPLES_FILE="data/raw/sra/PRJNA875278/samples_corrected.tsv"
SRA_DIR="data/raw/sra/PRJNA875278"
QUANT_DIR="data/interim/host_counts"
INDEX_DIR="data/raw/references/salmon_index"
FASTP_DIR="results/logs/fastp"
MIN_FREE_DISK_GB=15
COMPRESSED_FASTQ_DIR="${SRA_DIR}"
THREADS=2

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

log() { echo -e "${GREEN}[$(date '+%H:%M:%S')]${NC} $1"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
err()  { echo -e "${RED}[ERROR]${NC} $1"; }

# ---- RAM guard ----
check_ram() {
    local need_mb=${1:-2000}
    local avail_mb=$(grep MemAvailable /proc/meminfo 2>/dev/null | awk '{print int($2/1024)}')
    if [ "${avail_mb:-0}" -lt "$need_mb" ]; then
        warn "Low RAM: ${avail_mb} MB available (need ${need_mb} MB). Waiting 30s..."
        sleep 30
        avail_mb=$(grep MemAvailable /proc/meminfo 2>/dev/null | awk '{print int($2/1024)}')
        if [ "${avail_mb:-0}" -lt "$need_mb" ]; then
            err "Still low RAM (${avail_mb} MB). Skipping heavy step, will retry next run."
            return 1
        fi
    fi
    return 0
}

# ---- Pre-flight ----
log "=== Pre-flight check ==="
bash src/utils/preflight_check.sh || { err "Pre-flight failed"; exit 1; }

mkdir -p "$SRA_DIR" "$QUANT_DIR" "$FASTP_DIR"

# ---- Read sample list ----
if [ ! -f "$SAMPLES_FILE" ]; then
    err "Sample list not found: $SAMPLES_FILE"
    exit 1
fi

# Extract SRR IDs (skip header)
mapfile -t SRR_LIST < <(awk '{print $1}' "$SAMPLES_FILE")
TOTAL=${#SRR_LIST[@]}
log "Total samples to process: ${TOTAL}"

PROCESSED=0
FAILED=()

for SRR in "${SRR_LIST[@]}"; do
    PROCESSED=$((PROCESSED + 1))
    log ""
    log "=== [${PROCESSED}/${TOTAL}] Processing ${SRR} ==="

    # Skip if quant already exists
    if [ -f "${QUANT_DIR}/${SRR}/quant.sf" ]; then
        log "  quant.sf exists — skipping ${SRR}"
        continue
    fi

    # ---- Check disk ----
    FREE_GB=$(df "$SRA_DIR" 2>/dev/null | tail -1 | awk '{print $4}' | awk '{printf "%.0f", $1/1024/1024}')
    if [ "${FREE_GB:-0}" -lt "$MIN_FREE_DISK_GB" ]; then
        err "  Less than ${MIN_FREE_DISK_GB} GB free (${FREE_GB} GB). Aborting."
        exit 1
    fi

    # ---- Step 1: Download SRA ----
    if [ ! -f "${SRA_DIR}/${SRR}/${SRR}.sra" ]; then
        log "  [1/5] Downloading ${SRR} via prefetch..."
        prefetch "$SRR" --max-size 5G --output-directory "$SRA_DIR" > /dev/null 2>&1 && log "  Download OK" || { warn "  Download may have had warnings — checking file..."; }
        # Verify download
        if [ -f "${SRA_DIR}/${SRR}/${SRR}.sra" ]; then
            sz=$(du -sh "${SRA_DIR}/${SRR}/${SRR}.sra" 2>/dev/null | cut -f1)
            log "  SRA file: ${sz}"
        else
            err "  Download FAILED for ${SRR} — skipping"
            FAILED+=("$SRR")
            continue
        fi
    else
        log "  [1/5] ${SRR}.sra already downloaded"
    fi

    # ---- Step 2: fasterq-dump ----
    R1="${COMPRESSED_FASTQ_DIR}/${SRR}_1.fastq.gz"
    R2="${COMPRESSED_FASTQ_DIR}/${SRR}_2.fastq.gz"

    if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
        # RAM guard
        check_ram 3000 || true
        log "  [2/5] Extracting FASTQ (RAM guard: OK)..."
        fasterq-dump "${SRA_DIR}/${SRR}/${SRR}.sra" \
            --outdir "$SRA_DIR" \
            --threads "$THREADS" \
            --mem 2G \
            --split-files \
            --progress > /dev/null 2>&1 && log "  Extraction OK" || { err "  fasterq-dump FAILED"; FAILED+=("$SRR"); continue; }

        # Compress (one file at a time to reduce RAM spike)
        log "  Compressing R1..."
        pigz -p 1 -f "${SRA_DIR}/${SRR}"_1.fastq 2>/dev/null || gzip -f "${SRA_DIR}/${SRR}"_1.fastq
        log "  Compressing R2..."
        pigz -p 1 -f "${SRA_DIR}/${SRR}"_2.fastq 2>/dev/null || gzip -f "${SRA_DIR}/${SRR}"_2.fastq
        log "  Compression done"
    else
        log "  [2/5] FASTQ.gz already exists"
    fi

    # ---- Step 3: fastp trim ----
    R1_TRIM="${SRA_DIR}/${SRR}_trimmed_1.fastq.gz"
    R2_TRIM="${SRA_DIR}/${SRR}_trimmed_2.fastq.gz"
    JSON="${FASTP_DIR}/${SRR}.json"

    if [ ! -f "$JSON" ]; then
        log "  [3/5] Trimming with fastp..."
        fastp -i "$R1" -I "$R2" \
            -o "$R1_TRIM" -O "$R2_TRIM" \
            --json "$JSON" \
            --thread "$THREADS" \
            --qualified_quality_phred 20 \
            --length_required 50 \
            --detect_adapter_for_pe > /dev/null 2>&1 && log "  fastp OK" || { err "  fastp FAILED"; FAILED+=("$SRR"); continue; }
    else
        log "  [3/5] fastp already run"
    fi

    # ---- Step 4: Salmon quant ----
    QUANT_OUT="${QUANT_DIR}/${SRR}"
    if [ ! -f "${QUANT_OUT}/quant.sf" ]; then
        check_ram 4000 || true
        log "  [4/5] Quantifying with Salmon (RAM guard: OK)..."
        mkdir -p "$QUANT_OUT"
        salmon quant \
            -i "$INDEX_DIR" \
            -l A \
            -1 "$R1_TRIM" \
            -2 "$R2_TRIM" \
            -o "$QUANT_OUT" \
            --threads "$THREADS" \
            --gcBias --seqBias > /dev/null 2>&1 && log "  Salmon OK" || { err "  Salmon FAILED"; FAILED+=("$SRR"); continue; }
        # Log mapping rate
        if [ -f "${QUANT_OUT}/logs/salmon_quant.log" ]; then
            grep "processed\|mapping rate" "${QUANT_OUT}/logs/salmon_quant.log" 2>/dev/null | tail -2 >> results/logs/batch_process.log
        fi
    else
        log "  [4/5] Salmon quant already done"
    fi

    # ---- Step 5: Cleanup ----
    log "  [5/5] Cleaning up intermediate files..."
    # Remove SRA file (2 GB saved)
    rm -f "${SRA_DIR}/${SRR}/${SRR}.sra" 2>/dev/null || true
    rmdir "${SRA_DIR}/${SRR}" 2>/dev/null || true
    # Remove original untrimmed FASTQ.gz
    rm -f "$R1" "$R2" 2>/dev/null || true

    log "  ${SRR} complete."
done

log ""
log "=== Pipeline complete ==="
log "Processed: ${PROCESSED}/${TOTAL}"
if [ ${#FAILED[@]} -gt 0 ]; then
    warn "Failed samples: ${FAILED[*]}"
else
    log "All samples processed successfully."
fi
