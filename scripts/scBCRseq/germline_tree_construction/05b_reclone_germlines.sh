#!/bin/bash
# =============================================================================
# SCRIPT 05b — Re-run CreateGermlines with --cloned flag
# =============================================================================
# WHY THIS STEP EXISTS:
#   CreateGermlines.py in step 03 ran WITHOUT --cloned because clone IDs
#   didn't exist yet. This produces slightly different germlines per sequence
#   within a clone, which causes formatClones() to fail in dowser.
#
#   Now that DefineClones has assigned clone_id to every sequence, we re-run
#   CreateGermlines WITH --cloned. This builds one consensus germline per clone
#   and assigns it identically to every member — required for tree building.
#
# ORDER: must run AFTER 05_define_clones.sh and BEFORE 06_build_trees.sh
#
# Usage: bash 05b_reclone_germlines.sh
# Run from: ${PROJECT_DIR}/scripts/
# =============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../project_config.sh"

IN_FILE="${PROJECT_DIR}/04_clones/cloned_clone-pass.tsv"
OUT_DIR="${PROJECT_DIR}/04_clones"
LOG="${PROJECT_DIR}/logs/05b_reclone_germlines.log"

echo "=============================================" | tee ${LOG}
echo " Step 05b — CreateGermlines with --cloned" | tee -a ${LOG}
echo "=============================================" | tee -a ${LOG}
echo " Started: $(date)" | tee -a ${LOG}
echo "" | tee -a ${LOG}

if [ ! -f "${IN_FILE}" ]; then
    echo " ✗ ERROR: Input not found: ${IN_FILE}" | tee -a ${LOG}
    echo "   Run 05_define_clones.sh first" | tee -a ${LOG}
    exit 1
fi

N_IN=$(tail -n +2 "${IN_FILE}" | wc -l)
echo " Input:    ${IN_FILE}" | tee -a ${LOG}
echo " Sequences: ${N_IN}" | tee -a ${LOG}
echo "" | tee -a ${LOG}

docker run --rm \
    -v "${PROJECT_DIR}:/project" \
    ${DOCKER_IMAGE} \
    bash -c "
        set -e

        echo 'Running CreateGermlines.py with --cloned flag...'
        echo ''

        CreateGermlines.py \
            -d /project/04_clones/cloned_clone-pass.tsv \
            -r /usr/local/share/germlines/imgt/human/vdj/ \
            -g dmask \
            --cloned \
            --outdir /project/04_clones/ \
            --outname cloned_germlined \
            --failed \
            --log /project/logs/05b_createdb_internal.log

        echo ''
        echo 'CreateGermlines.py done. Output:'
        ls -lh /project/04_clones/cloned_germlined* 2>/dev/null || echo 'No output files found'
    " 2>&1 | tee -a ${LOG}

echo "" | tee -a ${LOG}

PASS_FILE="${OUT_DIR}/cloned_germlined_germ-pass.tsv"
FAIL_FILE="${OUT_DIR}/cloned_germlined_germ-fail.tsv"

if [ ! -f "${PASS_FILE}" ]; then
    echo " ✗ ERROR: Output not found: ${PASS_FILE}" | tee -a ${LOG}
    exit 1
fi

N_PASS=$(tail -n +2 "${PASS_FILE}" | wc -l)
echo " ✓ Germline pass: ${N_PASS} sequences" | tee -a ${LOG}

if [ -f "${FAIL_FILE}" ]; then
    N_FAIL=$(tail -n +2 "${FAIL_FILE}" | wc -l)
    echo " ✗ Germline fail: ${N_FAIL} sequences" | tee -a ${LOG}
fi

echo "" | tee -a ${LOG}
echo "=============================================" | tee -a ${LOG}
echo " STEP 05b COMPLETE" | tee -a ${LOG}
echo " Finished: $(date)" | tee -a ${LOG}
echo " Output: ${PASS_FILE}" | tee -a ${LOG}
echo " NEXT: bash 06_build_trees.sh" | tee -a ${LOG}
echo "=============================================" | tee -a ${LOG}
