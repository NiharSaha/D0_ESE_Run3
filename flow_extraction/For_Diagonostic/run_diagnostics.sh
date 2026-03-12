#!/bin/bash
# run_diagnostics.sh
# Runs save_mass_distributions for ALL 1% cent bins and ALL pT bins,
# captures the diagnostic output, then calls the Python parser for a summary.
#
# Usage:
#   ./run_diagnostics.sh [executable_path] [log_file]
# Defaults:
#   executable_path = ./save_mass_distributions
#   log_file        = diagnostic_log.txt

EXEC="${1:-./save_mass_distributions}"
LOGFILE="${2:-diagnostic_log.txt}"
PARSER="$(dirname "$0")/parse_diagnostics.py"

if [[ ! -x "$EXEC" ]]; then
    echo "ERROR: executable not found or not executable: $EXEC"
    echo "Build it first, e.g.:  g++ -o save_mass_distributions save_mass_distributions.C \$(root-config --cflags --libs)"
    exit 1
fi

echo "========================================================"
echo "  Running $EXEC (all cent/pT bins)"
echo "  Log -> $LOGFILE"
echo "========================================================"
"$EXEC" 2>&1 | tee "$LOGFILE"

echo ""
echo "========================================================"
echo "  Parsing diagnostic log: $LOGFILE"
echo "========================================================"
python3 "$PARSER" "$LOGFILE"
