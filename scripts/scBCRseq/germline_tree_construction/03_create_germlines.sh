#!/bin/bash
# =============================================================================
# SCRIPT 03 — IgBLAST + MakeDb + CreateGermlines (v3 — no heredocs)
# =============================================================================
# Usage: bash 03_create_germlines.sh
# Run from: ${PROJECT_DIR}/scripts/
# =============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../project_config.sh"

IGBLAST_DIR="${PROJECT_DIR}/03_germlines/igblast"
MAKEDB_DIR="${PROJECT_DIR}/03_germlines/makedb"
GERM_DIR="${PROJECT_DIR}/03_germlines"
LOG="${PROJECT_DIR}/logs/03_germlines.log"

mkdir -p "${IGBLAST_DIR}" "${MAKEDB_DIR}" "${PROJECT_DIR}/logs"

echo "=============================================" | tee ${LOG}
echo " Step 03 — IgBLAST + MakeDb + CreateGermlines" | tee -a ${LOG}
echo "=============================================" | tee -a ${LOG}
echo " Started: $(date)" | tee -a ${LOG}
echo "" | tee -a ${LOG}

# =============================================================================
# STEP 3A — Extract sequences to FASTA using pre-written Python helper
# =============================================================================
echo "─── 3A: Extracting sequences to FASTA ──────────" | tee -a ${LOG}

# The helper script is mounted from the scripts directory
docker run --rm \
    -v "${PROJECT_DIR}:/project" \
    -v "${SCRIPT_DIR}:/scripts" \
    ${DOCKER_IMAGE} \
    python3 /scripts/helper_extract_fasta.py \
    2>&1 | tee -a ${LOG}

# Verify FASTA was created
FASTA_FILE="${IGBLAST_DIR}/heavy_chain.fasta"
if [ ! -f "${FASTA_FILE}" ]; then
    echo " ✗ ERROR: FASTA file not created: ${FASTA_FILE}" | tee -a ${LOG}
    exit 1
fi

N_SEQS=$(grep -c "^>" "${FASTA_FILE}" || true)
echo "  ✓ FASTA created: ${N_SEQS} sequences" | tee -a ${LOG}
echo "" | tee -a ${LOG}

# =============================================================================
# STEP 3B — IgBLAST via AssignGenes.py
# =============================================================================
echo "─── 3B: Running IgBLAST (AssignGenes.py) ───────" | tee -a ${LOG}
echo "  Threads: ${NPROC}" | tee -a ${LOG}
echo "  Sequences: ${N_SEQS}" | tee -a ${LOG}
echo "  (Expect 10-30 min for ~20k sequences)" | tee -a ${LOG}
echo "" | tee -a ${LOG}

docker run --rm \
    -v "${PROJECT_DIR}:/project" \
    ${DOCKER_IMAGE} \
    bash -c "
        set -e

        echo '  Checking IMGT references...'
        ls /usr/local/share/germlines/imgt/human/vdj/

        echo ''
        echo '  Running AssignGenes.py...'
        AssignGenes.py igblast \
            -s /project/03_germlines/igblast/heavy_chain.fasta \
            -b /usr/local/share/igblast \
            --organism human \
            --loci ig \
            --format blast \
            --outdir /project/03_germlines/igblast/ \
            --nproc ${NPROC}

        echo ''
        echo '  AssignGenes.py complete. Files created:'
        ls -lh /project/03_germlines/igblast/
    " 2>&1 | tee -a ${LOG}

# Verify fmt7 output
FMT7=$(find "${IGBLAST_DIR}" -name "*.fmt7" 2>/dev/null | head -1)
if [ -z "${FMT7}" ]; then
    echo "" | tee -a ${LOG}
    echo " ✗ ERROR: No .fmt7 file found after IgBLAST" | tee -a ${LOG}
    echo "   Files in igblast dir:" | tee -a ${LOG}
    ls "${IGBLAST_DIR}" | tee -a ${LOG}
    exit 1
fi
echo "  ✓ IgBLAST output: ${FMT7}" | tee -a ${LOG}
echo "" | tee -a ${LOG}

# =============================================================================
# STEP 3C — MakeDb.py
# =============================================================================
echo "─── 3C: Running MakeDb.py ───────────────────────" | tee -a ${LOG}

# Get the fmt7 filename relative to project dir for Docker path
FMT7_REL=$(echo "${FMT7}" | sed "s|${PROJECT_DIR}|/project|")

docker run --rm \
    -v "${PROJECT_DIR}:/project" \
    ${DOCKER_IMAGE} \
    bash -c "
        set -e

        echo '  IgBLAST input: ${FMT7_REL}'

        MakeDb.py igblast \
            -i '${FMT7_REL}' \
            -s /project/03_germlines/igblast/heavy_chain.fasta \
            -r /usr/local/share/germlines/imgt/human/vdj/ \
            --regions default --format airr --extended \
            --failed \
            --outdir /project/03_germlines/makedb/ \
            --outname heavy_igblast

        echo ''
        echo '  MakeDb.py complete. Files created:'
        ls -lh /project/03_germlines/makedb/
    " 2>&1 | tee -a ${LOG}

# Verify MakeDb output
MAKEDB_PASS="${MAKEDB_DIR}/heavy_igblast_db-pass.tsv"
if [ ! -f "${MAKEDB_PASS}" ]; then
    echo "" | tee -a ${LOG}
    echo " ✗ ERROR: MakeDb db-pass.tsv not found" | tee -a ${LOG}
    echo "   Check MakeDb output above for errors" | tee -a ${LOG}
    exit 1
fi
N_MAKEDB=$(tail -n +2 "${MAKEDB_PASS}" | wc -l)
echo "  ✓ MakeDb pass: ${N_MAKEDB} sequences" | tee -a ${LOG}
echo "" | tee -a ${LOG}

# =============================================================================
# STEP 3D — Merge sample metadata using pre-written Python helper
# =============================================================================
echo "─── 3D: Merging sample metadata ─────────────────" | tee -a ${LOG}

docker run --rm \
    -v "${PROJECT_DIR}:/project" \
    -v "${SCRIPT_DIR}:/scripts" \
    ${DOCKER_IMAGE} \
    python3 /scripts/helper_merge_meta.py \
    2>&1 | tee -a ${LOG}

META_FILE="${MAKEDB_DIR}/heavy_igblast_with_meta.tsv"
if [ ! -f "${META_FILE}" ]; then
    echo " ✗ ERROR: Metadata-merged file not created" | tee -a ${LOG}
    exit 1
fi
echo "  ✓ Metadata merged" | tee -a ${LOG}
echo "" | tee -a ${LOG}

# =============================================================================
# STEP 3E — CreateGermlines.py
# =============================================================================
echo "─── 3E: Running CreateGermlines.py ─────────────" | tee -a ${LOG}

docker run --rm \
    -v "${PROJECT_DIR}:/project" \
    ${DOCKER_IMAGE} \
    bash -c "
        set -e

        IN='/project/03_germlines/makedb/heavy_igblast_with_meta.tsv'
        ROWS=\$(tail -n +2 \$IN | wc -l)

        echo '  Input file: '\$IN
        echo '  Input rows: '\$ROWS

        CreateGermlines.py \
            -d \$IN \
            -r /usr/local/share/germlines/imgt/human/vdj/ \
            -g dmask \
            --outdir /project/03_germlines/ \
            --outname germlined \
            --failed \
            --log /project/logs/03_createdb_internal.log

        echo ''
        echo '  CreateGermlines.py complete. Output:'
        ls -lh /project/03_germlines/germlined* 2>/dev/null || echo '  No germlined files found'
    " 2>&1 | tee -a ${LOG}

echo "" | tee -a ${LOG}

# =============================================================================
# FINAL CHECK
# =============================================================================
PASS_FILE="${GERM_DIR}/germlined_germ-pass.tsv"
FAIL_FILE="${GERM_DIR}/germlined_germ-fail.tsv"

if [ -f "${PASS_FILE}" ]; then
    N_PASS=$(tail -n +2 "${PASS_FILE}" | wc -l)
    echo " ✓ Germline pass: ${N_PASS} sequences" | tee -a ${LOG}
else
    echo " ✗ ERROR: germlined_germ-pass.tsv not found" | tee -a ${LOG}
    echo "   Full log: ${LOG}" | tee -a ${LOG}
    exit 1
fi

if [ -f "${FAIL_FILE}" ]; then
    N_FAIL=$(tail -n +2 "${FAIL_FILE}" | wc -l)
    echo " ✗ Germline fail: ${N_FAIL} sequences" | tee -a ${LOG}
    N_FAIL=0
fi

echo "" | tee -a ${LOG}
echo "=============================================" | tee -a ${LOG}
echo " STEP 03 COMPLETE" | tee -a ${LOG}
echo " Finished: $(date)" | tee -a ${LOG}
echo " Pass: ${N_PASS} sequences ready for tree building" | tee -a ${LOG}
echo " NEXT: bash 04_distance_threshold.sh" | tee -a ${LOG}
echo "=============================================" | tee -a ${LOG}
