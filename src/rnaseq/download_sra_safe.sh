#!/usr/bin/env bash
# ================================================================
# SAFE SRA DOWNLOAD — Track 3 Host-Microbe Integration
# ================================================================
# Downloads PRJNA875278 ONE sample at a time with size checks.
# Designed for WSL/Ubuntu with limited RAM (7.7 GB) and disk.
#
# Usage: bash src/rnaseq/download_sra_safe.sh <sample_list.tsv> <output_dir>
# ================================================================
set -euo pipefail

SAMPLE_LIST="${1:-config/samples_prjna875278.tsv}"
OUTDIR="${2:-data/raw/sra/PRJNA875278}"

# ---- Guardrail Constants ----
MAX_RETRIES=3
MIN_FREE_DISK_GB=10
PAUSE_BETWEEN_SEC=5

RED='\033[0;31m'
YELLOW='\033[1;33m'
GREEN='\033[0;32m'
NC='\033[0m'

log_info()  { echo -e "${GREEN}[INFO]${NC}  $(date '+%H:%M:%S') $1"; }
log_warn()  { echo -e "${YELLOW}[WARN]${NC}  $(date '+%H:%M:%S') $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $(date '+%H:%M:%S') $1"; }

# --- Pre-flight disk check ---
check_disk() {
    local free_kb
    free_kb=$(df "$OUTDIR" 2>/dev/null | tail -1 | awk '{print $4}' || echo 0)
    local free_gb=$((free_kb / 1024 / 1024))
    if [ "$free_gb" -lt "$MIN_FREE_DISK_GB" ]; then
        log_error "Less than ${MIN_FREE_DISK_GB} GB free (${free_gb} GB). Aborting."
        exit 1
    fi
    log_info "Disk OK: ${free_gb} GB free"
}

# --- Download ONE run ---
download_one_run() {
    local srr="$1"
    local outdir="$2"

    log_info "Downloading ${srr} ..."

    # prefetch first (small metadata), then fasterq-dump
    for attempt in $(seq 1 $MAX_RETRIES); do
        log_info "  Attempt ${attempt}/${MAX_RETRIES}"

        if prefetch --max-size 50G "$srr" --output-directory "$outdir/.prefetch_cache" 2>&1 | tail -3; then
            log_info "  prefetch OK for ${srr}"
            break
        else
            log_warn "  prefetch failed for ${srr} (attempt ${attempt})"
            if [ "$attempt" -eq "$MAX_RETRIES" ]; then
                log_error "Giving up on ${srr} after ${MAX_RETRIES} attempts"
                return 1
            fi
            sleep $((attempt * 30))  # exponential backoff
        fi
    done

    # fasterq-dump (this is the heavy part)
    log_info "  fasterq-dump ${srr} ..."
    if fasterq-dump "$srr" \
        --outdir "$outdir" \
        --threads 2 \
        --mem 2G \
        --split-files \
        --progress 2>&1 | tail -5; then
        log_info "  fasterq-dump OK for ${srr}"

        # Compress to save space
        log_info "  Compressing ${srr} FASTQ files ..."
        pigz -p 2 "$outdir/${srr}"_*.fastq 2>/dev/null || gzip "$outdir/${srr}"_*.fastq
        log_info "  ${srr} complete."
    else
        log_error "fasterq-dump failed for ${srr}"
        return 1
    fi
}

# --- Main ---
mkdir -p "$OUTDIR"
mkdir -p "$OUTDIR/.prefetch_cache"
mkdir -p "$OUTDIR/logs"

log_info "=== Safe SRA Download ==="
log_info "Sample list: ${SAMPLE_LIST}"
log_info "Output dir:  ${OUTDIR}"
log_info "Max concurrent downloads: 1 (serial only)"
echo ""

check_disk

# Read sample list
if [ ! -f "$SAMPLE_LIST" ]; then
    log_error "Sample list not found: ${SAMPLE_LIST}"
    log_info "Run metadata fetch first: bash src/rnaseq/fetch_sra_metadata.sh"
    exit 1
fi

# Count total
total=$(grep -c '^SRR' "$SAMPLE_LIST" 2>/dev/null || echo "0")
log_info "Total runs to download: ${total}"
echo ""

current=0
failed_runs=""

while IFS=$'\t' read -r srr _; do
    # Skip header/comments
    [[ "$srr" =~ ^#.* ]] && continue
    [[ "$srr" =~ ^SRR ]] || continue
    [ -z "$srr" ] && continue

    current=$((current + 1))
    log_info "=== [${current}/${total}] ${srr} ==="

    # Skip if already downloaded
    if compgen -G "$OUTDIR/${srr}*.fastq.gz" > /dev/null 2>&1; then
        log_info "  Already exists, skipping."
        continue
    fi

    # Check disk before EACH download
    check_disk

    # Download ONE sample
    if ! download_one_run "$srr" "$OUTDIR"; then
        failed_runs="${failed_runs} ${srr}"
        log_warn "  ${srr} failed, continuing with next sample."
    fi

    # Pause to let filesystem sync
    log_info "  Pausing ${PAUSE_BETWEEN_SEC}s ..."
    sleep "$PAUSE_BETWEEN_SEC"

done < "$SAMPLE_LIST"

echo ""
log_info "=== Download Complete ==="
log_info "Downloaded: $((current - $(echo "$failed_runs" | wc -w))) / ${total}"

if [ -n "$failed_runs" ]; then
    log_warn "Failed runs:${failed_runs}"
    log_info "Re-run this script to retry failed downloads."
fi
