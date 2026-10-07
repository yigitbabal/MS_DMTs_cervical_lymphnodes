#!/bin/bash
# =============================================================================
# SCRIPT 01b — Docker wrapper
# Usage: bash 01b_filter_by_seurat_barcodes.sh
# Run from: ${PROJECT_DIR}/scripts/
# =============================================================================
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../project_config.sh"

echo "============================================="
echo " Step 01b — Filter to Seurat Barcodes"
echo "============================================="

docker run --rm \
    -v "${PROJECT_DIR}:/project" \
    -v "${SCRIPT_DIR}:/scripts" \
    -e PROJECT_DIR="/project" \
    ${DOCKER_IMAGE} \
    Rscript /scripts/01b_filter_by_seurat_barcodes.R \
    2>&1 | tee "${PROJECT_DIR}/logs/01b_seurat_filter.log"

echo ""
echo " ✓ Done. Check logs/01b_seurat_filter.log"
echo " NEXT: bash 02_filter_heavy_chain.sh"
echo "       (update IN_FILE in 02_filter_heavy_chain.R first)"
