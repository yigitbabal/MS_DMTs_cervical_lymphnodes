#!/usr/bin/env Rscript
# =============================================================================
# SCRIPT 04 — Distance-to-Nearest Analysis (Clone Threshold Selection)
# =============================================================================
# Before defining clones, you MUST determine the correct CDR3 distance
# threshold. This script:
#   1. Calculates nearest-neighbor CDR3 distances
#   2. Plots the bimodal distribution
#   3. Suggests a threshold (you must visually confirm)
#
# The valley between the two peaks = your threshold
#   Left peak  = sequences in the SAME clone
#   Right peak = sequences in DIFFERENT clones
#
# Usage: Rscript 04_distance_threshold.R
# =============================================================================

suppressPackageStartupMessages({
    library(shazam)
    library(alakazam)
    library(ggplot2)
    library(gridExtra)
    library(dplyr)
})

PROJECT_DIR <- Sys.getenv("PROJECT_DIR",
                           unset = file.path(Sys.getenv("HOME"),
                                             "bcr_germline_project"))

IN_FILE  <- file.path(PROJECT_DIR, "03_germlines", "germlined_germ-pass.tsv")
OUT_DIR  <- file.path(PROJECT_DIR, "04_clones")
PLOT_DIR <- file.path(PROJECT_DIR, "06_results", "qc_plots")

dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=============================================\n")
cat(" Step 04 — Distance Threshold Analysis\n")
cat("=============================================\n")
cat(" Started:", format(Sys.time()), "\n\n")

# =============================================================================
# 1. LOAD GERMLINED TABLE
# =============================================================================
cat("[1/3] Loading germlined table...\n")
db <- readChangeoDb(IN_FILE)
cat("  Sequences loaded:", nrow(db), "\n\n")

# =============================================================================
# 2. CALCULATE DISTANCE TO NEAREST NEIGHBOR
# =============================================================================
cat("[2/3] Calculating nearest-neighbor distances...\n")
cat("  (This may take 5-15 minutes for large datasets)\n\n")

# distToNearest groups sequences by:
#   - Same V gene family
#   - Same J gene
#   - Same CDR3 length
# Then calculates hamming distance between CDR3 sequences
dist_db <- distToNearest(
    db,
    sequenceColumn = "junction",
    vCallColumn    = "v_call",
    jCallColumn    = "j_call",
    model          = "ham",        # Hamming distance
    normalize      = "len",        # Normalize by CDR3 length
    nproc          = 1,            # increase if you have multiple cores
    fields         = NULL          # calculate across all samples combined
)

cat("  ✓ Distance calculation complete\n\n")

# =============================================================================
# 3. PLOT AND SUGGEST THRESHOLD
# =============================================================================
cat("[3/3] Plotting distance distribution...\n")

# Main distance histogram
p_dist <- ggplot(dist_db %>% filter(!is.na(dist_nearest)),
                 aes(x = dist_nearest)) +
    geom_histogram(bins = 60, fill = "steelblue", color = "white", alpha = 0.8) +
    geom_vline(xintercept = 0.3,  linetype = "dashed", color = "red",   linewidth=1) +
    geom_vline(xintercept = 0.32, linetype = "dashed", color = "orange", linewidth=1) +
    geom_vline(xintercept = 0.35,  linetype = "dashed", color = "green",  linewidth=1) +
    annotate("text", x=0.3, y=Inf, label="0.10", vjust=2, color="red",    size=3) +
    annotate("text", x=0.32, y=Inf, label="0.15", vjust=2, color="orange", size=3) +
    annotate("text", x=0.35, y=Inf, label="0.20", vjust=2, color="green",  size=3) +
    labs(
        title    = "Distance to Nearest Neighbor",
        subtitle = "Valley between peaks = your clone threshold\nRed=0.3, Orange=0.32, Green=0.35",
        x        = "Normalized Hamming Distance",
        y        = "Count"
    ) +
    theme_bw(base_size = 12)

# Per-tissue breakdown
p_tissue <- ggplot(dist_db %>% filter(!is.na(dist_nearest)),
                   aes(x = dist_nearest)) +
    geom_histogram(bins = 60, alpha = 0.6, position = "identity") +
    labs(
        title = "Distance Distribution by Tissue",
        x     = "Normalized Hamming Distance",
        y     = "Count",
        fill  = "Tissue"
    ) +
    theme_bw(base_size = 12)

# Save plots
ggsave(file.path(PLOT_DIR, "04_distance_to_nearest.pdf"),
       gridExtra::arrangeGrob(p_dist, p_tissue, ncol=1), width = 10, height = 12)

# Auto-detect valley (simple method — visual confirmation required)
hist_data   <- hist(dist_db$dist_nearest[!is.na(dist_db$dist_nearest)],
                    breaks = 60, plot = FALSE)
# Find the valley in range 0.05-0.35
idx_range   <- which(hist_data$mids >= 0.05 & hist_data$mids <= 0.35)
valley_idx  <- idx_range[which.min(hist_data$counts[idx_range])]
auto_thresh <- round(hist_data$mids[valley_idx], 3)

# Save distance table for manual inspection
write.table(dist_db, 
            file.path(OUT_DIR, "dist_to_nearest.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

cat("\n=============================================\n")
cat(" DISTANCE ANALYSIS COMPLETE\n")
cat("=============================================\n")
cat(" Auto-detected threshold: ", auto_thresh, "\n")
cat("\n")
cat(" !! IMPORTANT: Open this plot and inspect:\n")
cat("   ", file.path(PLOT_DIR, "04_distance_to_nearest.pdf"), "\n")
cat("\n")
cat(" Look for the VALLEY between two peaks:\n")
cat("   Left peak  = sequences in the SAME clone\n")
cat("   Right peak = sequences in DIFFERENT clones\n")
cat("   Valley     = your threshold\n")
cat("\n")
cat(" Typical human BCR threshold: 0.10 - 0.20\n")
cat(" For MS B cells:              often ~0.15\n")
cat("\n")
cat(" Edit CLONE_DIST in project_config.sh then run:\n")
cat(" bash 05_define_clones.sh\n")
cat("=============================================\n")
