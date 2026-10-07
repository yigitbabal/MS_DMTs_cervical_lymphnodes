#!/bin/bash
# =============================================================================
# SCRIPT 05 — Define Clones with Immcantation
# =============================================================================
# Groups sequences into B cell clones based on:
#   - Same V gene (family level)
#   - Same J gene
#   - Same CDR3 length
#   - CDR3 sequence similarity within your threshold
#
# This is cross-sample cloning — detects the SAME clone across all 16 samples!
# Cell Ranger clonotypes are per-sample only and cannot do this.
#
# Usage: bash 05_define_clones.sh
# =============================================================================

set -e
source "$(dirname "$0")/../project_config.sh"

IN_FILE="${PROJECT_DIR}/03_germlines/germlined_germ-pass.tsv"
OUT_DIR="${PROJECT_DIR}/04_clones"
LOG="${PROJECT_DIR}/logs/05_define_clones.log"

echo "============================================="
echo " Step 05 — Define Clones"
echo "============================================="
echo " Started:    $(date)"       | tee ${LOG}
echo " Threshold:  ${CLONE_DIST}" | tee -a ${LOG}
echo ""

docker run --rm \
    -v "${PROJECT_DIR}:/project" \
    ${DOCKER_IMAGE} \
    bash -c "
        set -e
        echo 'Running DefineClones.py...'
        echo 'Threshold: ${CLONE_DIST}'
        echo ''
        
        DefineClones.py \
            -d /project/03_germlines/germlined_germ-pass.tsv \
            --act set \
            --model ham \
            --norm len \
            --dist ${CLONE_DIST} \
            --sym min \
            --outdir /project/04_clones/ \
            --outname cloned \
            --log /project/logs/05_define_clones_defineclones.log
        
        echo ''
        echo 'DefineClones complete.'
    " 2>&1 | tee -a ${LOG}

echo "" | tee -a ${LOG}

# Summary stats
if [ -f "${OUT_DIR}/cloned_clone-pass.tsv" ]; then
    TOTAL=$(tail -n +2 "${OUT_DIR}/cloned_clone-pass.tsv" | wc -l)
    
    # Count unique clones
    CLONE_COL=$(head -1 "${OUT_DIR}/cloned_clone-pass.tsv" | \
                tr '\t' '\n' | grep -n "^clone_id$" | cut -d: -f1)
    N_CLONES=$(tail -n +2 "${OUT_DIR}/cloned_clone-pass.tsv" | \
               cut -f${CLONE_COL} | sort -u | wc -l)
    
    echo " Total sequences:    ${TOTAL}"  | tee -a ${LOG}
    echo " Unique clones:      ${N_CLONES}" | tee -a ${LOG}
else
    echo " ✗ ERROR: cloned_clone-pass.tsv not found!" | tee -a ${LOG}
    exit 1
fi

echo "" | tee -a ${LOG}
echo "=============================================" | tee -a ${LOG}
echo " CLONE DEFINITION COMPLETE" | tee -a ${LOG}
echo " Finished: $(date)" | tee -a ${LOG}
echo " NEXT: Run bash 06_build_trees.sh" | tee -a ${LOG}
echo "=============================================" | tee -a ${LOG}
