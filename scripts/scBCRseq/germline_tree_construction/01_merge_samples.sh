#!/bin/bash
# =============================================================================
# SCRIPT 01 — Merge 16 AIRR Files with Sample Barcode Prefixing
# =============================================================================
# This script:
#   1. Reads each airr_rearrangement.tsv from 00_raw_airr/
#   2. Adds sample_id prefix to every cell_id and sequence_id
#   3. Adds sample metadata columns (patient_id, tissue, group)
#   4. Merges into one single AIRR table
#
# Usage: bash 01_merge_samples.sh
# =============================================================================

set -e

# Load project config
source "$(dirname "$0")/../project_config.sh"

MANIFEST="${PROJECT_DIR}/sample_manifest.tsv"
RAW_DIR="${PROJECT_DIR}/00_raw_airr"
OUT_DIR="${PROJECT_DIR}/01_merged"
LOG="${PROJECT_DIR}/logs/01_merge.log"

echo "============================================="
echo " Step 01 — Merge 16 AIRR Samples"
echo "============================================="
echo " Started: $(date)" | tee ${LOG}
echo ""

# Check manifest exists
if [ ! -f "${MANIFEST}" ]; then
    echo "ERROR: sample_manifest.tsv not found at ${MANIFEST}"
    echo "Run 00_setup_environment.sh first"
    exit 1
fi

# Count samples
N_SAMPLES=$(tail -n +2 ${MANIFEST} | wc -l)
echo " Found ${N_SAMPLES} samples in manifest" | tee -a ${LOG}
echo ""

# Output file
MERGED_OUT="${OUT_DIR}/merged_airr_all.tsv"

# Flag to track if header has been written
HEADER_WRITTEN=false

# Process each sample
while IFS=$'\t' read -r sample_id patient_id tissue group airr_file; do
    
    # Skip header
    [ "$sample_id" = "sample_id" ] && continue

    AIRR_PATH="${RAW_DIR}/${airr_file}"

    # Check file exists
    if [ ! -f "${AIRR_PATH}" ]; then
        echo "  ⚠ WARNING: File not found — ${AIRR_PATH}" | tee -a ${LOG}
        echo "             Skipping sample ${sample_id}" | tee -a ${LOG}
        continue
    fi

    # Count input rows
    N_ROWS=$(tail -n +2 "${AIRR_PATH}" | wc -l)
    echo "  Processing ${sample_id} | ${tissue} | ${patient_id} | ${N_ROWS} contigs" | tee -a ${LOG}

    # Write header once (add new metadata columns at the end)
    if [ "$HEADER_WRITTEN" = false ]; then
        head -1 "${AIRR_PATH}" | \
            awk 'BEGIN{OFS="\t"} {print $0"\tsample_id\tpatient_id\ttissue\tgroup"}' \
            > "${MERGED_OUT}"
        HEADER_WRITTEN=true
    fi

    # Process data rows:
    #   - Prefix cell_id (col 1) with sample_id
    #   - Prefix sequence_id (col 3) with sample_id
    #   - Append metadata columns
    tail -n +2 "${AIRR_PATH}" | \
        awk -v sid="${sample_id}" \
            -v pid="${patient_id}" \
            -v tis="${tissue}" \
            -v grp="${group}" \
        'BEGIN{FS=OFS="\t"} {
            $1 = sid"_"$1    # prefix cell_id
            $3 = sid"_"$3    # prefix sequence_id
            print $0"\t"sid"\t"pid"\t"tis"\t"grp
        }' >> "${MERGED_OUT}"

done < "${MANIFEST}"

echo "" | tee -a ${LOG}

# Summary statistics
TOTAL_ROWS=$(tail -n +2 "${MERGED_OUT}" | wc -l)
TOTAL_CELLS=$(tail -n +2 "${MERGED_OUT}" | awk -F'\t' '$22=="T"' | cut -f1 | sort -u | wc -l)

echo "=============================================" | tee -a ${LOG}
echo " MERGE COMPLETE" | tee -a ${LOG}
echo "=============================================" | tee -a ${LOG}
echo " Total contigs (rows):  ${TOTAL_ROWS}" | tee -a ${LOG}
echo " Confirmed cells:       ${TOTAL_CELLS}" | tee -a ${LOG}
echo " Output file:           ${MERGED_OUT}" | tee -a ${LOG}
echo " Finished: $(date)" | tee -a ${LOG}
echo ""
echo " NEXT: Run bash 02_filter_heavy_chain.sh"
