#!/usr/bin/env bash
# ================================================================
# PRE-FLIGHT SYSTEM CHECK — Track 3 Host-Microbe Integration
# ================================================================
# Runs BEFORE any pipeline step to verify the system has enough
# resources. Aborts with a clear message if constraints aren't met.
#
# Usage: bash src/utils/preflight_check.sh [--skip-disk]
# ================================================================
set -euo pipefail

# ---- Configurable thresholds ----
MIN_FREE_RAM_MB=2000       # 2 GB free RAM minimum
MIN_FREE_DISK_GB=10        # 10 GB free disk minimum
WARN_SWAP_MB=4000          # Warn if swap < 4 GB
MAX_LOAD=8                 # Warn if system load avg > this

RED='\033[0;31m'
YELLOW='\033[1;33m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

PASS=0
WARN=0
FAIL=0

check_pass()  { echo -e "  ${GREEN}[PASS]${NC} $1"; PASS=$((PASS+1)); }
check_warn()  { echo -e "  ${YELLOW}[WARN]${NC} $1"; WARN=$((WARN+1)); }
check_fail()  { echo -e "  ${RED}[FAIL]${NC} $1"; FAIL=$((FAIL+1)); }

echo "========================================"
echo " Track 3 Pre-Flight System Check"
echo " $(date '+%Y-%m-%d %H:%M:%S')"
echo "========================================"
echo ""

# --- 1. RAM ---
echo "--- Memory ---"
total_ram_kb=$(grep MemTotal /proc/meminfo 2>/dev/null | awk '{print $2}' || echo 0)
avail_ram_kb=$(grep MemAvailable /proc/meminfo 2>/dev/null | awk '{print $2}' || echo 0)
total_ram_mb=$((total_ram_kb / 1024))
avail_ram_mb=$((avail_ram_kb / 1024))

echo "  Total RAM:    ${total_ram_mb} MB"
echo "  Available:    ${avail_ram_mb} MB"

if [ "$avail_ram_mb" -lt "$MIN_FREE_RAM_MB" ]; then
    check_fail "Available RAM (${avail_ram_mb} MB) < minimum (${MIN_FREE_RAM_MB} MB)"
    echo "  → ACTION: Close other applications (browser, IDE tabs) and retry."
else
    check_pass "Available RAM (${avail_ram_mb} MB) >= minimum (${MIN_FREE_RAM_MB} MB)"
fi

# --- 2. Swap ---
echo ""
echo "--- Swap ---"
swap_total_kb=$(grep SwapTotal /proc/meminfo 2>/dev/null | awk '{print $2}' || echo 0)
swap_total_mb=$((swap_total_kb / 1024))

echo "  Swap total:   ${swap_total_mb} MB"

if [ "$swap_total_mb" -lt "$WARN_SWAP_MB" ]; then
    check_warn "Swap (${swap_total_mb} MB) is small (< ${WARN_SWAP_MB} MB). OOM risk is HIGH."
    echo "  → If the system crashes, increase WSL swap via .wslconfig:"
    echo "    [wsl2]"
    echo "    memory=8GB"
    echo "    swap=8GB"
else
    check_pass "Swap (${swap_total_mb} MB) >= ${WARN_SWAP_MB} MB"
fi

# --- 3. Disk Space ---
if [ "${1:-}" != "--skip-disk" ]; then
    echo ""
    echo "--- Disk Space ---"
    available_kb=$(df . 2>/dev/null | tail -1 | awk '{print $4}')
    available_gb=$((available_kb / 1024 / 1024))

    echo "  Free on $(pwd): ${available_gb} GB"

    if [ "$available_gb" -lt "$MIN_FREE_DISK_GB" ]; then
        check_fail "Free disk (${available_gb} GB) < minimum (${MIN_FREE_DISK_GB} GB)"
    else
        check_pass "Free disk (${available_gb} GB) >= ${MIN_FREE_DISK_GB} GB"
    fi
fi

# --- 4. System Load ---
echo ""
echo "--- System Load ---"
load=$(uptime 2>/dev/null | awk -F'load average:' '{print $2}' | awk '{print $1}' | tr -d ',' || echo "unknown")
echo "  Load average: ${load}"

if [ "$load" != "unknown" ]; then
    load_int=$(echo "$load" | cut -d. -f1)
    if [ "$load_int" -gt "$MAX_LOAD" ] 2>/dev/null; then
        check_warn "System load (${load}) is high. Pipeline may run slowly."
    else
        check_pass "System load (${load}) is acceptable"
    fi
fi

# --- 5. Required Tools ---
echo ""
echo "--- Required Tools ---"
for tool in conda mamba python R fastp salmon; do
    if command -v "$tool" &>/dev/null; then
        check_pass "$tool is available"
    else
        check_warn "$tool not found in PATH (may be in conda env)"
    fi
done

# --- SUMMARY ---
echo ""
echo "========================================"
echo " SUMMARY: ${GREEN}${PASS} passed${NC}, ${YELLOW}${WARN} warnings${NC}, ${RED}${FAIL} failed${NC}"
echo "========================================"

if [ "$FAIL" -gt 0 ]; then
    echo ""
    echo -e "${RED}Cannot proceed: ${FAIL} critical check(s) failed.${NC}"
    echo "Fix the issues above and re-run this check."
    exit 1
elif [ "$WARN" -gt 0 ]; then
    echo ""
    echo -e "${YELLOW}Proceeding with ${WARN} warning(s).${NC}"
    echo "The pipeline will use conservative resource limits."
    exit 0
else
    echo ""
    echo -e "${GREEN}All checks passed. System is ready.${NC}"
    exit 0
fi
