#!/bin/bash
# =============================================================================
# negishi_scratch_keepalive.sh
#
# Purpose: Prevent Purdue Negishi scratch file purging by "touching" all files
#          in specified directories. Negishi purges files not accessed in ~14
#          days. Running this script (e.g. via cron) resets the atime/mtime on
#          every file, keeping them alive.
#
# Usage:
#   bash negishi_scratch_keepalive.sh          # normal run
#   bash negishi_scratch_keepalive.sh --dry-run  # preview only, no changes
#
# Scheduling (add to crontab with: crontab -e):
#   0 9 * * 1,4  /path/to/negishi_scratch_keepalive.sh >> ~/scratch_keepalive.log 2>&1
#   (runs every Monday and Thursday at 9am — safely within the 14-day window)
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# CONFIGURATION — edit these paths to match your analysis input directories
# ---------------------------------------------------------------------------
PATHS_TO_KEEP=(
    # Primary scratch working area
    "/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/Quantiles_MB0to1_Apr12/ROOT"
    "/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/Resolution_MB11to21_Jan26/ROOT/"
    "/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/Resolution_MB0to1_Apr27_charge_pub_v0/ROOT/"
    "/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/Resolution_MB0to1_Apr15_charge_pub_v0/ROOT/"
    # Add your specific analysis input directories below:
    # "/scratch/negishi/$USER/ePIC/ToF_material_scan/input"
    # "/scratch/negishi/$USER/CMS/data/Run3"
    # "/scratch/negishi/$USER/EIC/simulations"
)

# ---------------------------------------------------------------------------
# OPTIONS
# ---------------------------------------------------------------------------
DRY_RUN=false
LOG_PREFIX="[$(date '+%Y-%m-%d %H:%M:%S')]"

# Parse flags
for arg in "$@"; do
    case "$arg" in
        --dry-run|-n)
            DRY_RUN=true
            ;;
        --help|-h)
            grep '^#' "$0" | head -20 | sed 's/^# \?//'
            exit 0
            ;;
    esac
done

# ---------------------------------------------------------------------------
# FUNCTIONS
# ---------------------------------------------------------------------------

log() {
    echo "$LOG_PREFIX $*"
}

keepalive_path() {
    local target_path="$1"

    if [[ ! -e "$target_path" ]]; then
        log "WARNING: Path does not exist, skipping: $target_path"
        return 1
    fi

    # Count files before touching
    local file_count
    file_count=$(find "$target_path" -type f 2>/dev/null | wc -l)

    log "Processing: $target_path  ($file_count files found)"

    if $DRY_RUN; then
        log "  [DRY RUN] Would touch $file_count files — no changes made."
        return 0
    fi

    # Touch all files recursively (updates atime + mtime)
    # -type f ensures we only touch regular files (not symlinks or dirs)
    local touched=0
    local failed=0

    while IFS= read -r -d '' file; do
        if touch "$file" 2>/dev/null; then
            (( touched++ )) || true
        else
            log "  WARNING: Could not touch: $file"
            (( failed++ )) || true
        fi
    done < <(find "$target_path" -type f -print0 2>/dev/null)

    log "  Done — touched: $touched  |  failed: $failed"
}

# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

log "====== Negishi Scratch Keepalive ======"
$DRY_RUN && log "*** DRY RUN MODE — no files will be modified ***"
log "User: $USER"
log "Host: $(hostname)"
log "Paths to process: ${#PATHS_TO_KEEP[@]}"
echo ""

total_paths=0
skipped_paths=0

for path in "${PATHS_TO_KEEP[@]}"; do
    # Skip comment-like or empty entries
    [[ -z "$path" || "$path" == \#* ]] && continue

    keepalive_path "$path" && (( total_paths++ )) || (( skipped_paths++ )) || true
    echo ""
done

log "====== Summary ======"
log "Paths processed : $total_paths"
log "Paths skipped   : $skipped_paths"
$DRY_RUN && log "No files were modified (dry run)."
log "====== Done ======"
