#!/bin/bash
# =============================================================================
# SCRIPT 06 — Build Germline Trees
# Usage: bash 06_build_trees.sh
# Run from: ${PROJECT_DIR}/scripts/
# =============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../project_config.sh"

LOG="${PROJECT_DIR}/logs/06_trees.log"
MIN_CLONE_SIZE=2   # increase to 5 to focus on expanded clones only

echo "=============================================" | tee ${LOG}
echo " Step 06 — Build Germline Trees" | tee -a ${LOG}
echo "=============================================" | tee -a ${LOG}
echo " Started:         $(date)" | tee -a ${LOG}
echo " Min clone size:  ${MIN_CLONE_SIZE}" | tee -a ${LOG}
echo " Threads:         ${NPROC}" | tee -a ${LOG}
echo "" | tee -a ${LOG}

# Quick clone size check before running Docker
CLONE_FILE="${PROJECT_DIR}/04_clones/cloned_clone-pass.tsv"
if [ ! -f "${CLONE_FILE}" ]; then
    echo " ✗ ERROR: Clone file not found: ${CLONE_FILE}" | tee -a ${LOG}
    exit 1
fi

TOTAL_SEQS=$(tail -n +2 "${CLONE_FILE}" | wc -l)
echo " Total sequences in clone file: ${TOTAL_SEQS}" | tee -a ${LOG}
echo "" | tee -a ${LOG}

docker run --rm \
    -v "${PROJECT_DIR}:/project" \
    -v "${SCRIPT_DIR}:/scripts" \
    -e PROJECT_DIR="/project" \
    -e MIN_CLONE_SIZE="${MIN_CLONE_SIZE}" \
    -e NPROC="${NPROC}" \
    -e CLONE_DIST="${CLONE_DIST}" \
    ${DOCKER_IMAGE} \
    Rscript /scripts/06_build_trees.R \
    2>&1 | tee -a ${LOG}

echo "" | tee -a ${LOG}
echo "=============================================" | tee -a ${LOG}
echo " STEP 06 COMPLETE" | tee -a ${LOG}
echo " Finished: $(date)" | tee -a ${LOG}
echo " NEXT: bash 07_annotate_visualize.sh" | tee -a ${LOG}
echo "=============================================" | tee -a ${LOG}
