#!/bin/bash
# =============================================================================
# SCRIPT 07 — Docker wrapper for annotation and visualization R script
# Usage: bash 07_annotate_visualize.sh
# Run from: ${PROJECT_DIR}/scripts/
# =============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../project_config.sh"

echo "============================================="
echo " Step 07 — Annotate Trees & Visualize"
echo "============================================="
echo " Started: $(date)"
echo ""

docker run --rm \
    -v "${PROJECT_DIR}:/project" \
    -v "${SCRIPT_DIR}:/scripts" \
    -e PROJECT_DIR="/project" \
    ${DOCKER_IMAGE} \
    Rscript /scripts/07_annotate_and_visualize.R \
    2>&1 | tee "${PROJECT_DIR}/logs/07_visualize.log"

echo ""
echo "============================================="
echo " STEP 07 COMPLETE"
echo "============================================="
echo " Results: ${PROJECT_DIR}/06_results/"
