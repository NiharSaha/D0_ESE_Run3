#!/bin/bash

# ---------------------------------------------------------------
# Compile first, then run save_mass.out for every (cent, pT) pair
# ---------------------------------------------------------------

SRCDIR="/home/saha115/D0_ESE/CMSSW_13_2_11/src/flow_extraction"
EXE="${SRCDIR}/save_mass.out"
LOGDIR="${SRCDIR}/diag_logs"
SUMMARY_LOG="${SRCDIR}/all_bins_summary.log"

# Number of 1% cent bins and pT bins (must match the macro)
N_CENT=90
N_PT=10

mkdir -p "$LOGDIR"

# ---------------------------------------------------------------
# Build
# ---------------------------------------------------------------
echo "==> Compiling save_mass_distributions.C ..."
cd "$SRCDIR"
g++ -O2 -o save_mass.out save_mass_distributions.C \
    $(root-config --cflags --libs) \
    -lMathMore 2>&1

if [ $? -ne 0 ]; then
    echo "ERROR: Compilation failed. Exiting."
    exit 1
fi
echo "==> Compilation successful."

# ---------------------------------------------------------------
# Run all (cent, pT) combinations
# ---------------------------------------------------------------
echo "" > "$SUMMARY_LOG"
echo "==========================================" >> "$SUMMARY_LOG"
echo "  Run summary: all (cent1pct, pT) bins"   >> "$SUMMARY_LOG"
echo "  $(date)"                                  >> "$SUMMARY_LOG"
echo "==========================================" >> "$SUMMARY_LOG"

TOTAL=$(( N_CENT * N_PT ))
DONE=0

for (( icen=0; icen<N_CENT; icen++ )); do
    for (( ipt=0; ipt<N_PT; ipt++ )); do
        LOGFILE="${LOGDIR}/log_cen${icen}_pt${ipt}.txt"
        echo -n "  Running: ./save_mass.out ${icen} ${ipt}  ->  ${LOGFILE} ... "

        "${EXE}" "${icen}" "${ipt}" > "${LOGFILE}" 2>&1
        STATUS=$?

        if [ $STATUS -ne 0 ]; then
            echo "FAILED (exit ${STATUS})"
            echo "FAILED  cen=${icen}  pt=${ipt}  exit=${STATUS}" >> "$SUMMARY_LOG"
        else
            echo "OK"
            # Extract the ratio lines and append to summary log
            echo "--- cen=${icen}  pt=${ipt} ---" >> "$SUMMARY_LOG"
            grep -E "(q2 cand max/min ratio|q3 cand max/min ratio|WARNING)" \
                "${LOGFILE}" >> "$SUMMARY_LOG"
        fi

        DONE=$(( DONE + 1 ))
        echo "  Progress: ${DONE} / ${TOTAL}"
    done
done

echo ""
echo "==> All jobs done. Summary saved to: ${SUMMARY_LOG}"
echo "==> Individual logs in:              ${LOGDIR}/"
echo ""
echo "==> To find max ratios run:"
echo "    python3 ${SRCDIR}/find_max_ratio.py ${SUMMARY_LOG}"
