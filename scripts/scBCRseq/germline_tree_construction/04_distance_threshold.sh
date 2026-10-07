#!/bin/bash
# =============================================================================
# SCRIPT 04 — Docker wrapper for distance threshold R script
# Usage: bash 04_distance_threshold.sh
# Run from: ${PROJECT_DIR}/scripts/
# =============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../project_config.sh"

echo "============================================="
echo " Step 04 — Distance Threshold Analysis"
echo "============================================="
echo " Started: $(date)"
echo ""

docker run --rm \
    -v "${PROJECT_DIR}:/project" \
    -v "${SCRIPT_DIR}:/scripts" \
    -e PROJECT_DIR="/project" \
    ${DOCKER_IMAGE} \
    Rscript /scripts/04_distance_threshold.R \
    2>&1 | tee "${PROJECT_DIR}/logs/04_threshold.log"

echo ""
echo "============================================="
echo " STEP 04 COMPLETE"
echo "============================================="
echo ""
echo " !! IMPORTANT — Do this before continuing:"
echo ""
echo " Open this PDF and find the VALLEY between two peaks:"
echo " ${PROJECT_DIR}/06_results/qc_plots/04_distance_to_nearest.pdf"
echo ""
echo " Left peak  = sequences in the SAME clone"
echo " Right peak = sequences in DIFFERENT clones"
echo " Valley     = your threshold value"
echo ""
echo " Then edit CLONE_DIST in:"
echo " ${PROJECT_DIR}/project_config.sh"
echo ""
echo " Then run: bash 05_define_clones.sh"
