#!/usr/bin/env Rscript
# =============================================================================
# SCRIPT 01b — Filter merged AIRR file to Seurat B cell barcodes only
# =============================================================================
# Run AFTER 01_merge_samples.sh and BEFORE 02_filter_heavy_chain.sh
#
# This removes:
#   - Non-B cells (T cells, NK, monocytes captured in VDJ but not GEX)
#   - Low quality cells filtered during Seurat QC
#   - Doublets removed in Seurat
#   - Cells not in your B cell Seurat object
#
# Result: only BCR sequences from cells with known Seurat annotations
# =============================================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
})

PROJECT_DIR <- Sys.getenv("PROJECT_DIR")

IN_FILE      <- file.path(PROJECT_DIR, "01_merged",  "merged_airr_all.tsv")
SEURAT_META  <- file.path(PROJECT_DIR, "seurat_barcodes_annotated.csv")
OUT_FILE     <- file.path(PROJECT_DIR, "01_merged",  "merged_airr_seurat_filtered.tsv")

cat("=============================================\n")
cat(" Step 01b — Filter to Seurat B Cell Barcodes\n")
cat("=============================================\n")
cat(" Started:", format(Sys.time()), "\n\n")

# =============================================================================
# Check inputs
# =============================================================================
if (!file.exists(SEURAT_META)) {
    stop(paste(
        "\nERROR: Seurat metadata not found at:", SEURAT_META,
        "\n\nExport from R with:",
        "\n  seurat_meta <- data.frame(",
        "\n      barcode   = colnames(your_b_cell_seurat),",
        "\n      cell_type = your_b_cell_seurat$immFine_cluster,",
        "\n      cluster_id= your_b_cell_seurat$seurat_clusters,",
        "\n      tissue    = your_b_cell_seurat$orig.ident",
        "\n  )",
        "\n  write.csv(seurat_meta,",
        "\n    '", file.path(PROJECT_DIR, "seurat_barcodes_annotated.csv"), "',",
        "\n    row.names=FALSE)"
    ))
}

# =============================================================================
# Load Seurat metadata
# =============================================================================
cat("[1/4] Loading Seurat metadata...\n")
seurat_meta <- read.csv(SEURAT_META, stringsAsFactors=FALSE)
cat("  Seurat cells:", nrow(seurat_meta), "\n")
cat("  Columns:", paste(colnames(seurat_meta), collapse=", "), "\n")
cat("  Cell types:\n")
print(table(seurat_meta$cell_type))
cat("\n")

# The barcode column — detect automatically
barcode_col <- intersect(c("barcode", "cell_id", "Barcode", "CellBarcode"),
                          colnames(seurat_meta))[1]
if (is.na(barcode_col)) {
    cat("  Available columns:", paste(colnames(seurat_meta), collapse=", "), "\n")
    stop("Cannot find barcode column. Expected: 'barcode' or 'cell_id'")
}
cat("  Barcode column:", barcode_col, "\n\n")

seurat_barcodes <- seurat_meta[[barcode_col]]

# =============================================================================
# Load merged AIRR
# =============================================================================
cat("[2/4] Loading merged AIRR file...\n")
airr <- read.delim(IN_FILE, stringsAsFactors=FALSE)
cat("  Total AIRR rows:", nrow(airr), "\n")
cat("  Unique cell_ids:", n_distinct(airr$cell_id), "\n\n")

# =============================================================================
# Show barcode format comparison BEFORE filtering
# =============================================================================
cat("[3/4] Checking barcode format compatibility...\n")
cat("  Example AIRR cell_id:  ", head(airr$cell_id, 3), "\n")
cat("  Example Seurat barcode:", head(seurat_barcodes, 3), "\n\n")

# Check overlap before any modification
direct_overlap <- sum(airr$cell_id %in% seurat_barcodes)
cat("  Direct overlap:", direct_overlap, "/", n_distinct(airr$cell_id), "AIRR cells\n")

if (direct_overlap == 0) {
    cat("\n  !! No direct overlap — checking barcode format...\n")
    
    # Try: AIRR has sample prefix but Seurat doesn't
    airr_stripped <- sub("^[^_]+_", "", airr$cell_id)
    overlap_stripped <- sum(airr_stripped %in% seurat_barcodes)
    cat("  After stripping AIRR prefix:", overlap_stripped, "matches\n")
    
    # Try: Seurat has prefix but AIRR doesn't  
    seurat_stripped <- sub("^[^_]+_", "", seurat_barcodes)
    overlap_seurat_strip <- sum(airr$cell_id %in% seurat_stripped)
    cat("  After stripping Seurat prefix:", overlap_seurat_strip, "matches\n")
    
    if (overlap_stripped > overlap_seurat_strip && overlap_stripped > 0) {
        cat("\n  Auto-fix: stripping prefix from AIRR cell_ids\n")
        airr$cell_id_match <- airr_stripped
        filter_col <- "cell_id_match"
    } else if (overlap_seurat_strip > 0) {
        cat("\n  Auto-fix: using prefix-stripped Seurat barcodes\n")
        seurat_barcodes <- seurat_stripped
        filter_col <- "cell_id"
    } else {
        cat("\n  !! WARNING: Cannot match barcodes automatically\n")
        cat("  Check that sample prefixes match between Seurat and AIRR files\n")
        cat("  Proceeding without filtering — output = input\n")
        file.copy(IN_FILE, OUT_FILE, overwrite=TRUE)
        cat("  Output:", OUT_FILE, "\n")
        quit(status=0)
    }
} else {
    filter_col <- "cell_id"
}

# =============================================================================
# Filter
# =============================================================================
cat("\n[4/4] Filtering AIRR to Seurat barcodes...\n")

airr_filtered <- airr %>%
    filter(.data[[filter_col]] %in% seurat_barcodes)

# Join Seurat cell type directly into AIRR table
seurat_join <- seurat_meta %>%
    select(all_of(c(barcode_col, "cell_type"))) %>%
    rename(cell_id_join = all_of(barcode_col))

# Match on the right column
if (filter_col == "cell_id_match") {
    airr_filtered <- airr_filtered %>%
        left_join(seurat_join, by=c("cell_id_match"="cell_id_join"))
} else {
    airr_filtered <- airr_filtered %>%
        left_join(seurat_join, by=c("cell_id"="cell_id_join"))
}

n_before <- nrow(airr)
n_after  <- nrow(airr_filtered)
n_removed <- n_before - n_after
pct_kept  <- round(n_after / n_before * 100, 1)

cat("  Before filter:", n_before, "rows\n")
cat("  After filter: ", n_after,  "rows (", pct_kept, "% kept)\n")
cat("  Removed:      ", n_removed,"rows\n\n")

cat("  Cell type distribution in filtered AIRR:\n")
if ("cell_type" %in% colnames(airr_filtered)) {
    print(table(airr_filtered$cell_type, useNA="ifany"))
}
cat("\n")

# Save
write.table(airr_filtered, OUT_FILE, sep="\t", row.names=FALSE, quote=FALSE)
cat("  ✓ Saved:", OUT_FILE, "\n\n")

cat("=============================================\n")
cat(" STEP 01b COMPLETE\n")
cat("=============================================\n")
cat(" Input rows: ", n_before, "\n")
cat(" Output rows:", n_after, "\n")
cat(" Finished:", format(Sys.time()), "\n\n")
cat(" IMPORTANT: Update downstream scripts to use:\n")
cat("   merged_airr_seurat_filtered.tsv\n")
cat(" instead of:\n")
cat("   merged_airr_all.tsv\n\n")
cat(" NEXT: bash 02_filter_heavy_chain.sh\n")
cat("   (after updating IN_FILE path in 02_filter_heavy_chain.R)\n")
cat("=============================================\n")
