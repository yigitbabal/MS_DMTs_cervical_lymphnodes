#!/usr/bin/env Rscript
# =============================================================================
# SCRIPT 07 — Annotate Trees with Seurat Labels & Visualize
# =============================================================================

suppressPackageStartupMessages({
    library(dowser)
    library(alakazam)
    library(dplyr)
    library(ggplot2)
    library(gridExtra)
    library(stringr)
})

# =============================================================================
# CONFIG
# =============================================================================
PROJECT_DIR      <- "/mnt/h3048/grp_laakso/Project_repository/Yigit/MS_FNA_followup/bcr_germline_project"
IGPHYML_MIN_SIZE <- 20

TREES_RDS   <- file.path(PROJECT_DIR, "05_trees",   "all_trees.rds")
SUMMARY_TSV <- file.path(PROJECT_DIR, "05_trees",   "tree_summary.tsv")
SEURAT_META <- file.path(PROJECT_DIR, "seurat_barcodes_annotated.csv")
PLOT_DIR    <- file.path(PROJECT_DIR, "06_results", "tree_plots")
TABLE_DIR   <- file.path(PROJECT_DIR, "06_results", "tables")

dir.create(PLOT_DIR, showWarnings=FALSE, recursive=TRUE)
dir.create(TABLE_DIR, showWarnings=FALSE, recursive=TRUE)

# =============================================================================
# Colour palette — matches your exact Seurat cell type labels
# =============================================================================
#CELLTYPE_COLORS <- c(
#    "Naive B"                 = "#4E79A7",
#    "SM"                      = "#F28E2B",   # Switched Memory
#    "USM"                     = "#76B7B2",   # Unswitched Memory
#    "DN"                      = "#BAB0AC",   # Double Negative
#    "IFN-stimulated B"        = "#E15759",
#    "Tbet+CD11c+"             = "#B07AA1",
#    "Pre-B"                   = "#59A14F",
#    "GC B cells (dark zone)"  = "#1B7837",   # Germinal Centre dark zone
#    "GC B cells (light zone)" = "#7FBF7B",   # Germinal Centre light zone
#    "Plasmablasts"            = "#D6604D",
#    "Unknown"                 = "#D3D3D3",
#    "Germline"                = "black"      # required by dowser
#)

CELLTYPE_COLORS <- c(
  "Naive B"                 = "#377EB8",   # Vibrant Blue
  "SM"                      = "#FF7F00",   # Bright Orange
  "USM"                     = "#4DAF4A",   # Bright Green
  "DN"                      = "#984EA3",   # Rich Purple
  "IFN-stimulated B"        = "#E41A1C",   # Strong Red
  "Tbet+CD11c+"             = "#A65628",   # Warm Brown
  "Pre-B"                   = "#F781BF",   # Hot Pink
  "GC B cells (dark zone)"  = "#1B9E77",   # Dark Teal
  "GC B cells (light zone)" = "#A6CEE3",   # Light Blue
  "Plasmablasts"            = "#E6AB02",   # Mustard / Dark Gold
  "Unknown"                 = "#CCCCCC",   # Neutral Light Grey
  "Germline"                = "black"      # Required by dowser
)

cat("=============================================\n")
cat(" Step 07 — Annotate Trees & Visualize\n")
cat("=============================================\n")
cat(" Started:", format(Sys.time()), "\n\n")

# =============================================================================
# 1. Load data
# =============================================================================
cat("[1/5] Loading trees and Seurat metadata...\n")

if (!file.exists(TREES_RDS))   stop("Trees RDS not found: ",        TREES_RDS)
if (!file.exists(SEURAT_META)) stop("Seurat metadata not found: ",  SEURAT_META)

result       <- readRDS(TREES_RDS)
seurat_meta  <- read.csv(SEURAT_META, stringsAsFactors=FALSE)
tree_summary <- read.delim(SUMMARY_TSV, stringsAsFactors=FALSE)

cat("  Trees loaded:  ", nrow(result), "\n")
cat("  Seurat cells:  ", nrow(seurat_meta), "\n")
cat("  Cell types:\n")
print(table(seurat_meta$cell_type))
cat("\n")

# =============================================================================
# 2. Barcode extraction helper
# =============================================================================
cat("[2/5] Matching barcodes to Seurat annotations...\n")

extract_cell_id <- function(seq_id) sub("_contig_[0-9]+$", "", seq_id)

seurat_lookup <- setNames(seurat_meta$cell_type, seurat_meta$barcode)

# Verify format on first clone
test_ids      <- result$data[[1]]@data$sequence_id
test_cell_ids <- extract_cell_id(test_ids)
direct_match  <- sum(test_cell_ids %in% names(seurat_lookup))
cat("  Format check (clone 1):", direct_match, "/", length(test_cell_ids), "matched\n\n")

# =============================================================================
# 3. Annotate every tree node
# =============================================================================
cat("[3/5] Annotating tree nodes...\n")

result_annotated <- result
total_cells   <- 0
matched_cells <- 0

for (i in seq_len(nrow(result_annotated))) {
    seq_ids    <- result_annotated$data[[i]]@data$sequence_id
    cell_ids   <- extract_cell_id(seq_ids)
    cell_types <- seurat_lookup[cell_ids]
    cell_types[is.na(cell_types)] <- "Unknown"
    sample_ids <- sub("_[ACGT]+-[0-9]+$", "", cell_ids)

    result_annotated$data[[i]]@data$cell_type <- cell_types
    result_annotated$data[[i]]@data$cell_id   <- cell_ids
    result_annotated$data[[i]]@data$sample_id <- sample_ids

    total_cells   <- total_cells   + length(cell_ids)
    matched_cells <- matched_cells + sum(cell_types != "Unknown")
}

cat("  Total cells:   ", total_cells, "\n")
cat("  Matched:       ", matched_cells,
    "(", round(matched_cells/total_cells*100, 1), "%)\n\n")
cat("  Cell type distribution:\n")
all_types <- unlist(lapply(seq_len(nrow(result_annotated)), function(i)
    result_annotated$data[[i]]@data$cell_type))
print(sort(table(all_types), decreasing=TRUE))
cat("\n")

# =============================================================================
# 4. Plot trees
# =============================================================================
cat("[4/5] Plotting trees...\n")

# plotTrees() returns a LIST of ggplots — always use [[1]] to extract
top_idx <- order(-tree_summary$n_cells)[1:min(20, nrow(tree_summary))]

# ── Overview: top 9 clones on one page ──
cat("  Creating overview PDF (top 9 clones)...\n")
pdf(file.path(PLOT_DIR, "00_top_clones_overview.pdf"), width=18, height=14)

plot_list <- lapply(top_idx[1:min(9, length(top_idx))], function(i) {
    tryCatch({
        plotTrees(result_annotated[i, ],
                  tips    = "cell_type",
                  tipsize = 2,
                  palette = CELLTYPE_COLORS)[[1]] +    # [[1]] extracts the ggplot
            labs(title    = paste0("Clone ", result_annotated$clone_id[i],
                                   " (n=", nrow(result_annotated$data[[i]]@data), ")"),
                 subtitle = result_annotated$data[[i]]@v_gene) +
            theme(plot.title    = element_text(size=9, face="bold"),
                  plot.subtitle = element_text(size=7),
                  legend.text   = element_text(size=7),
                  legend.title  = element_text(size=8))
    }, error=function(e) {
        cat("  Plot error clone", i, ":", conditionMessage(e), "\n")
        NULL
    })
})

plot_list <- Filter(Negate(is.null), plot_list)
if (length(plot_list) > 0) gridExtra::grid.arrange(grobs=plot_list, ncol=3)
dev.off()
cat("  ✓ Overview PDF:", file.path(PLOT_DIR, "00_top_clones_overview.pdf"), "\n")

# ── Individual 2-panel PDFs: cell_type + sample ──
cat("  Creating individual PDFs (top", length(top_idx), "clones)...\n")
n_saved <- 0

for (i in top_idx) {
    tryCatch({
        clone_id <- result_annotated$clone_id[i]
        n_cells  <- nrow(result_annotated$data[[i]]@data)
        v_gene   <- result_annotated$data[[i]]@v_gene

        p1 <- plotTrees(result_annotated[i, ],
                        tips    = "cell_type",
                        tipsize = 2,
                        palette = CELLTYPE_COLORS)[[1]] +
            labs(title    = paste0("Clone ", clone_id, " | n=", n_cells),
                 subtitle = paste0(v_gene, " | B cell type"))

        p2 <- plotTrees(result_annotated[i, ],
                        tips    = "sample_id",
                        tipsize = 2)[[1]] +
            labs(title    = paste0("Clone ", clone_id, " | n=", n_cells),
                 subtitle = paste0(v_gene, " | Sample"))

        # gridExtra instead of patchwork to avoid S3/list conflicts
        combined <- gridExtra::arrangeGrob(p1, p2, ncol=2)
        ggsave(
            file.path(PLOT_DIR,
                      paste0("clone_", clone_id, "_n", n_cells, ".pdf")),
            combined, width=14, height=8
        )
        n_saved <- n_saved + 1
    }, error=function(e) cat("  Error clone", i, ":", conditionMessage(e), "\n"))
}
cat("  ✓", n_saved, "individual PDFs saved to:", PLOT_DIR, "\n\n")

# =============================================================================
# 5. Clone composition table
# =============================================================================
cat("[5/5] Saving clone composition table...\n")

clone_composition <- do.call(rbind, lapply(seq_len(nrow(result_annotated)), function(i) {
    tryCatch({
        d <- result_annotated$data[[i]]@data
        data.frame(
            clone_id      = as.character(result_annotated$clone_id[i]),
            n_cells       = nrow(d),
            v_gene        = result_annotated$data[[i]]@v_gene,
            j_gene        = result_annotated$data[[i]]@j_gene,
            junc_len      = result_annotated$data[[i]]@junc_len,
            method        = ifelse(nrow(d) >= IGPHYML_MIN_SIZE, "igphyml", "pml"),
            dominant_type = names(sort(table(d$cell_type), decreasing=TRUE))[1],
            cell_types    = paste(sort(unique(d$cell_type)), collapse=";"),
            n_samples     = n_distinct(d$sample_id),
            samples       = paste(sort(unique(d$sample_id)), collapse=";"),
            stringsAsFactors = FALSE
        )
    }, error=function(e) NULL)
}))

clone_composition <- clone_composition[order(-clone_composition$n_cells), ]

write.table(clone_composition,
            file.path(TABLE_DIR, "clone_composition.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)
cat("  ✓ Saved:", file.path(TABLE_DIR, "clone_composition.tsv"), "\n\n")

# =============================================================================
# MS-relevant summary
# =============================================================================
cat("=============================================\n")
cat(" MS-relevant findings\n")
cat("=============================================\n\n")

cat("Top 15 largest clones:\n")
print(head(clone_composition[, c("clone_id","n_cells","v_gene",
                                  "dominant_type","n_samples")], 15))

cat("\nCross-sample clones (same clone across multiple samples):\n")
cross <- clone_composition[clone_composition$n_samples > 1, ]
cat("  Found:", nrow(cross), "\n")
if (nrow(cross) > 0) {
    print(head(cross[, c("clone_id","n_cells","n_samples",
                          "dominant_type","v_gene","samples")], 10))
}

cat("\nCell type composition in large (igphyml) clones:\n")
large <- clone_composition[clone_composition$method == "igphyml", ]
cat("  Large clones:", nrow(large), "\n")
print(sort(table(large$dominant_type), decreasing=TRUE))

cat("\n=============================================\n")
cat(" STEP 07 COMPLETE\n")
cat("=============================================\n")
cat(" Finished:", format(Sys.time()), "\n")
cat(" Plots:  ", PLOT_DIR, "\n")
cat(" Tables: ", TABLE_DIR, "\n\n")
cat(" Key output files:\n")
cat("   00_top_clones_overview.pdf  — 9-clone overview\n")
cat("   clone_XXX_nYYY.pdf          — individual clone trees\n")
cat("   clone_composition.tsv       — full clone annotation table\n")
