#!/usr/bin/env Rscript
# =============================================================================
# SCRIPT 02 — Filter Heavy Chains & QC (FIXED — robust column handling)
# =============================================================================

suppressPackageStartupMessages({
    library(alakazam)
    library(dplyr)
    library(ggplot2)
    library(stringr)
    library(patchwork)
})

PROJECT_DIR <- Sys.getenv("PROJECT_DIR",
                           unset = file.path(Sys.getenv("HOME"),
                                             "bcr_germline_project"))

IN_FILE <- file.path(PROJECT_DIR, "01_merged", "merged_airr_seurat_filtered.tsv")
OUT_DIR  <- file.path(PROJECT_DIR, "02_filtered")
PLOT_DIR <- file.path(PROJECT_DIR, "06_results", "qc_plots")

dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=============================================\n")
cat(" Step 02 — Filter Heavy Chains & QC\n")
cat("=============================================\n")
cat(" Started:", format(Sys.time()), "\n\n")

# =============================================================================
# 1. LOAD DATA
# =============================================================================
cat("[1/5] Loading merged AIRR table...\n")
db <- read.delim(IN_FILE, stringsAsFactors = FALSE)
cat("  Total contigs loaded:", nrow(db), "\n")
cat("  Columns:", ncol(db), "\n")
cat("  Column names:\n ")
cat(paste(colnames(db), collapse="\n  "), "\n\n")

# =============================================================================
# 2. NORMALISE is_cell — 10x uses TRUE/FALSE, "T"/"F", or 1/0 depending on 
#    Cell Ranger version
# =============================================================================
cat("[2/5] Normalising column formats...\n")

if ("is_cell" %in% colnames(db)) {
    cat("  is_cell raw values:", paste(unique(db$is_cell), collapse=", "), "\n")
    cat("  is_cell class:     ", class(db$is_cell), "\n")
    
    db$is_cell_flag <- db$is_cell %in% c(TRUE, "T", "TRUE", "true", 1, "1")
    cat("  is_cell TRUE count:", sum(db$is_cell_flag), "/", nrow(db), "\n")
} else {
    cat("  is_cell not found — treating all rows as cells\n")
    db$is_cell_flag <- TRUE
}

# Normalise locus — derive from v_call if missing
if (!"locus" %in% colnames(db)) {
    cat("  locus column missing — deriving from v_call\n")
    db$locus <- case_when(
        grepl("^IGHV", db$v_call) ~ "IGH",
        grepl("^IGKV", db$v_call) ~ "IGK",
        grepl("^IGLV", db$v_call) ~ "IGL",
        TRUE ~ NA_character_
    )
} else {
    cat("  locus values:", paste(unique(db$locus), collapse=", "), "\n")
}

# Normalise productive
if ("productive" %in% colnames(db)) {
    cat("  productive raw values:", paste(unique(db$productive), collapse=", "), "\n")
    db$productive_flag <- db$productive %in% c(TRUE, "T", "TRUE", "true", 1, "1")
} else {
    db$productive_flag <- TRUE
}

cat("\n")

# =============================================================================
# 3. FILTER STEP BY STEP
# =============================================================================
cat("[3/5] Filtering...\n")
n_start <- nrow(db)

# is_cell
db <- db %>% filter(is_cell_flag == TRUE)
cat("  After is_cell filter:          ", nrow(db),
    " (removed:", n_start - nrow(db), ")\n")

if (nrow(db) == 0) {
    stop(paste(
        "\nERROR: 0 rows remain after is_cell filter.",
        "\nYour is_cell values were:", paste(unique(db$is_cell), collapse=", "),
        "\nRun: bash diagnose_airr.sh  to inspect the merged file."
    ))
}

# IGH only
n_before  <- nrow(db)
db_heavy  <- db %>% filter(locus == "IGH")
cat("  After IGH filter:              ", nrow(db_heavy),
    " (removed:", n_before - nrow(db_heavy), " light chains)\n")

if (nrow(db_heavy) == 0) {
    stop(paste(
        "\nERROR: 0 IGH rows found.",
        "\nLocus values present:", paste(unique(db$locus), collapse=", "),
        "\nRun: bash diagnose_airr.sh"
    ))
}

# Productive
n_before  <- nrow(db_heavy)
db_heavy  <- db_heavy %>% filter(productive_flag == TRUE)
cat("  After productive filter:       ", nrow(db_heavy),
    " (removed:", n_before - nrow(db_heavy), ")\n")

# V and J assigned
n_before  <- nrow(db_heavy)
db_heavy  <- db_heavy %>%
    filter(!is.na(v_call), v_call != "",
           !is.na(j_call), j_call != "")
cat("  After V+J assigned filter:     ", nrow(db_heavy),
    " (removed:", n_before - nrow(db_heavy), ")\n")

# Junction length — calculate if column missing
if (!"junction_length" %in% colnames(db_heavy) & "junction" %in% colnames(db_heavy)) {
    db_heavy <- db_heavy %>% mutate(junction_length = nchar(junction))
}

if ("junction_length" %in% colnames(db_heavy)) {
    n_before <- nrow(db_heavy)
    db_heavy <- db_heavy %>%
        filter(!is.na(junction_length),
               junction_length >= 9,
               junction_length <= 150)
    cat("  After junction length filter:  ", nrow(db_heavy),
        " (removed:", n_before - nrow(db_heavy), ")\n")
}

# One IGH per cell — keep highest UMI/consensus support
umi_col <- intersect(c("consensus_count", "umi_count", "duplicate_count"),
                     colnames(db_heavy))[1]
if (!is.na(umi_col)) {
    n_before <- nrow(db_heavy)
    db_heavy <- db_heavy %>%
        group_by(cell_id) %>%
        slice_max(order_by = .data[[umi_col]], n = 1, with_ties = FALSE) %>%
        ungroup()
    cat("  After 1 IGH per cell:          ", nrow(db_heavy),
        " (removed:", n_before - nrow(db_heavy), " multi-heavy)\n")
}

cat("\n  FINAL heavy chain table:", nrow(db_heavy), "cells\n\n")

# =============================================================================
# 4. DERIVED COLUMNS
# =============================================================================
cat("[4/5] Adding derived columns...\n")

db_heavy <- db_heavy %>%
    mutate(
        v_family      = str_extract(v_call, "IGHV[0-9]+"),
        isotype       = str_extract(c_call, "IGH[A-Z][0-9]?"),
        class_switched = !isotype %in% c("IGHM", "IGHD", NA)
    )

cat("  Isotype distribution:\n")
print(table(db_heavy$isotype, useNA = "ifany"))
cat("\n")

# =============================================================================
# 5. QC PLOTS
# =============================================================================
cat("[5/5] Generating QC plots...\n")

plots <- list()

if ("sample_id" %in% colnames(db_heavy)) {
    plots$cells <- db_heavy %>%
        count(sample_id) %>%
        ggplot(aes(reorder(sample_id, n), n)) +
        geom_col(fill="steelblue") + coord_flip() +
        labs(title="IGH Cells Per Sample", x="Sample", y="Count") +
        theme_bw(base_size=9)

    plots$iso <- db_heavy %>%
        filter(!is.na(isotype)) %>%
        ggplot(aes(sample_id, fill=isotype)) +
        geom_bar(position="fill") +
        scale_y_continuous(labels=scales::percent) +
        labs(title="Isotype Distribution", x=NULL, y="Proportion") +
        theme_bw(base_size=9) +
        theme(axis.text.x=element_text(angle=45,hjust=1))
}

if ("junction_length" %in% colnames(db_heavy)) {
    plots$cdr3 <- ggplot(db_heavy, aes(junction_length)) +
        geom_histogram(bins=40, fill="steelblue", color="white") +
        labs(title="CDR3 Length", x="Length (nt)", y="Count") +
        theme_bw(base_size=9)
}

plots$vgene <- db_heavy %>%
    filter(!is.na(v_family)) %>%
    count(v_family) %>%
    mutate(prop = n/sum(n)) %>%
    filter(prop > 0.02) %>%
    ggplot(aes(reorder(v_family, prop), prop)) +
    geom_col(fill="steelblue") + coord_flip() +
    scale_y_continuous(labels=scales::percent) +
    labs(title="V Family Usage (>2%)", x=NULL, y="Proportion") +
    theme_bw(base_size=9)

if (length(plots) >= 2) {
    ggsave(file.path(PLOT_DIR, "02_qc_overview.pdf"),
           wrap_plots(plots, ncol=2), width=14, height=10)
    cat("  ✓ QC plots saved\n\n")
}

# =============================================================================
# SAVE
# =============================================================================
OUT_FILE <- file.path(OUT_DIR, "heavy_chain_filtered.tsv")
write.table(db_heavy, OUT_FILE, sep="\t", row.names=FALSE, quote=FALSE)

cat("=============================================\n")
cat(" FILTER COMPLETE\n")
cat("=============================================\n")
cat(" Output:", OUT_FILE, "\n")
cat(" Cells: ", nrow(db_heavy), "\n\n")

if ("sample_id" %in% colnames(db_heavy)) {
    cat(" Per-sample summary:\n")
    db_heavy %>%
        group_by(sample_id) %>%
        summarise(n_cells=n(),
                  pct_switched=round(mean(class_switched, na.rm=TRUE)*100,1),
                  .groups="drop") %>%
        print(n=20)
}

cat("\n Finished:", format(Sys.time()), "\n")
cat(" NEXT: bash 03_create_germlines.sh\n\n")
