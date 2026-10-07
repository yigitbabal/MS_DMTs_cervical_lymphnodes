#!/bin/bash
# =============================================================================
# SCRIPT 02 — Docker wrapper to run filter R script
# =============================================================================
# Usage: bash 02_filter_heavy_chain.sh
# Must be run from ${PROJECT_DIR}/scripts/
# =============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${SCRIPT_DIR}/../project_config.sh"

if [ ! -f "${CONFIG}" ]; then
    echo "ERROR: project_config.sh not found at: ${CONFIG}"
    echo "Make sure you are running from: \${PROJECT_DIR}/scripts/"
    exit 1
fi

source "${CONFIG}"

echo "============================================="
echo " Step 02 — Filter Heavy Chains (Docker)"
echo "============================================="
echo " Project dir:  ${PROJECT_DIR}"
echo " Docker image: ${DOCKER_IMAGE}"
echo ""

docker run --rm \
    -v "${PROJECT_DIR}:/project" \
    -e PROJECT_DIR="/project" \
    ${DOCKER_IMAGE} \
    Rscript /project/scripts/02_filter_heavy_chain.R \
    2>&1 | tee "${PROJECT_DIR}/logs/02_filter.log"

echo ""
echo " ✓ Done. Check: ${PROJECT_DIR}/logs/02_filter.log"
echo " NEXT: bash 03_create_germlines.sh"
