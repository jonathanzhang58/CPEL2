#!/bin/bash
#
# submit_cpel_pipeline.sh
# Usage: ./submit_cpel_pipeline.sh <tile_width> <p_value>
#

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <tile_width> <p_value>"
  exit 1
fi

TILE_WIDTH=$1
P_VALUE=$2

# Paths to your two sbatch scripts
PARALLEL_SCRIPT="/dcs07/afeinber/data/personal/jzhan/CPEL2/parallelcutpromoters.sh"
ANALYZE_SCRIPT="/dcs07/afeinber/data/personal/jzhan/CPEL2/analyze.sh"

# Base output dir created by parallelcut.sh
OUTPUT_BASE="/dcs07/afeinber/data/personal/jzhan/CPEL2/output"
# This glob will expand at job‐start time to the exact directory
OUTPUT_PATTERN="${OUTPUT_BASE}/cpel_old_young_${TILE_WIDTH}_*"

echo "1) Submitting parallelcut (width=${TILE_WIDTH})..."
JOB1=$(sbatch --parsable "${PARALLEL_SCRIPT}" "${TILE_WIDTH}")
echo "   → parallelcut job submitted as ${JOB1}"

echo "2) Submitting analyze (p_val=${P_VALUE}) after completion of ${JOB1}..."
JOB2=$(sbatch --parsable \
  --dependency=afterok:"${JOB1}" \
  "${ANALYZE_SCRIPT}" "${OUTPUT_PATTERN}" "${P_VALUE}")
echo "   → analyze job submitted as ${JOB2} (afterok:${JOB1})"

echo "All set! First job: ${JOB1}, second job: ${JOB2}"
