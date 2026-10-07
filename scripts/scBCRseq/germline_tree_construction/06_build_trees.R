#!/usr/bin/env Rscript
# =============================================================================
# SCRIPT 06 — Build Germline Trees (v10 — tiered strategy)
# =============================================================================
# Strategy:
#   Large clones (>= IGPHYML_MIN_SIZE cells): igphyml (slow, SHM-aware)
#   Small clones (2 to IGPHYML_MIN_SIZE-1):   pml    (fast, good enough)
#
# Biologically: only large clones have meaningful tree topology anyway.
# 2-cell clones produce trivial 2-tip trees regardless of method.
# =============================================================================

suppressPackageStartupMessages({
    library(dowser)
    library(alakazam)
    library(dplyr)
    library(ape)
    library(stringr)
    library(Biostrings)
})

# ── EDIT THESE TWO LINES FOR YOUR SYSTEM ──────────────────────────────
PROJECT_DIR      <- "/mnt/h3048/grp_laakso/Project_repository/Yigit/MS_FNA_followup/bcr_germline_project"
NPROC            <- 20        # set to your number of CPU cores
IGPHYML_EXEC     <- "/home/babayigi/igphyml/src/igphyml"  # update if different
# ───────────────────────────────────────────────────────────────────────
MIN_CLONE_SIZE   <- 2
IGPHYML_MIN_SIZE <- 20    # clones >= this size use igphyml; smaller use pml

IN_FILE <- file.path(PROJECT_DIR, "04_clones", "cloned_germlined_germ-pass.tsv")
OUT_DIR <- file.path(PROJECT_DIR, "05_trees")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=============================================\n")
cat(" Building Germline Trees (v10 — tiered)\n")
cat("=============================================\n")
cat(" igphyml for clones >=", IGPHYML_MIN_SIZE, "cells\n")
cat(" pml     for clones <",  IGPHYML_MIN_SIZE, "cells\n")
cat(" Started:", format(Sys.time()), "\n\n")

# =============================================================================
# Helper: detect stop codons
# =============================================================================
has_stop_codon <- function(seq_nt) {
    if (is.na(seq_nt) || nchar(seq_nt) < 3) return(FALSE)
    seq_clean <- gsub("[\\.-]", "", seq_nt)
    trim_len  <- floor(nchar(seq_clean) / 3) * 3
    if (trim_len < 3) return(FALSE)
    tryCatch({
        aa <- as.character(
            Biostrings::translate(Biostrings::DNAString(substr(seq_clean, 1, trim_len)),
                                  if.fuzzy.codon="solve"))
        grepl("\\*", aa)
    }, error=function(e) FALSE)
}

# =============================================================================
# 1. Load
# =============================================================================
cat("[1/7] Loading data...\n")
db <- readChangeoDb(IN_FILE)
cat("  Sequences:", nrow(db), "\n")
cat("  Clones:   ", n_distinct(db$clone_id), "\n\n")

# =============================================================================
# 2. Germline column
# =============================================================================
germ_col <- grep("germline_alignment_d_mask", colnames(db), value=TRUE)[1]
if (is.na(germ_col)) germ_col <- "germline_alignment"
cat("[2/7] Germline column:", germ_col, "\n\n")

# =============================================================================
# 3. Filter by clone size
# =============================================================================
cat("[3/7] Filtering clones (min =", MIN_CLONE_SIZE, ")...\n")
clone_sizes  <- db %>% count(clone_id, name="n_cells")
large_clones <- clone_sizes %>% filter(n_cells >= MIN_CLONE_SIZE)
db_filtered  <- db %>% filter(clone_id %in% large_clones$clone_id)

size_table <- clone_sizes %>%
    mutate(bin = case_when(
        n_cells == 1  ~ "singleton",
        n_cells == 2  ~ "2",
        n_cells <= 5  ~ "3-5",
        n_cells <= 10 ~ "6-10",
        n_cells <= 20 ~ "11-20",
        n_cells <= 50 ~ "21-50",
        TRUE          ~ ">50")) %>%
    count(bin)
print(size_table)
cat("  Clones:", nrow(large_clones), "| Sequences:", nrow(db_filtered), "\n\n")

write.table(clone_sizes %>% arrange(desc(n_cells)),
            file.path(OUT_DIR, "clone_sizes.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

# =============================================================================
# 4. Remove stop codon sequences
# =============================================================================
cat("[4/7] Checking for stop codons...\n")
stop_flags <- sapply(db_filtered$sequence_alignment, has_stop_codon)
n_stop <- sum(stop_flags)
cat("  Stop codon sequences found:", n_stop, "\n")

if (n_stop > 0) {
    db_filtered <- db_filtered[!stop_flags, ]
    clone_sizes2  <- db_filtered %>% count(clone_id, name="n_cells")
    large_clones2 <- clone_sizes2 %>% filter(n_cells >= MIN_CLONE_SIZE)
    db_filtered   <- db_filtered %>% filter(clone_id %in% large_clones2$clone_id)
    cat("  Sequences after removal:", nrow(db_filtered), "\n")
    cat("  Clones after removal:   ", nrow(large_clones2), "\n")
} else {
    cat("  None found ✓\n")
}
cat("\n")

# =============================================================================
# 5. Split into large (igphyml) and small (pml) clone sets
# =============================================================================
cat("[5/7] Splitting by clone size...\n")
final_sizes <- db_filtered %>% count(clone_id, name="n_cells")

large_ids <- final_sizes %>% filter(n_cells >= IGPHYML_MIN_SIZE) %>% pull(clone_id)
small_ids <- final_sizes %>% filter(n_cells >= MIN_CLONE_SIZE,
                                     n_cells <  IGPHYML_MIN_SIZE) %>% pull(clone_id)

db_large <- db_filtered %>% filter(clone_id %in% large_ids)
db_small <- db_filtered %>% filter(clone_id %in% small_ids)

cat("  igphyml clones (>=", IGPHYML_MIN_SIZE, "cells):", length(large_ids),
    "| sequences:", nrow(db_large), "\n")
cat("  pml     clones (2 to", IGPHYML_MIN_SIZE-1, "cells):", length(small_ids),
    "| sequences:", nrow(db_small), "\n\n")

# =============================================================================
# 6. formatClones + getTrees for LARGE clones with igphyml
# =============================================================================
result_large <- NULL
if (length(large_ids) > 0) {
    cat("[6/7] Building trees for LARGE clones with igphyml...\n")
    cat("  Clones:", length(large_ids), "\n")
    cat("  Started:", format(Sys.time()), "\n\n")

    clones_large <- suppressWarnings(
        formatClones(db_large,
                     seq="sequence_alignment", germ=germ_col,
                     clone="clone_id", v_call="v_call", j_call="j_call",
                     junc_len="junction_length", minseq=2, nproc=1)
    )
    cat("  Formatted:", nrow(clones_large), "clones\n")

    if (nrow(clones_large) > 0) {
        result_large <- tryCatch(
            getTrees(clones_large, build="igphyml",
                     exec=IGPHYML_EXEC, nproc=NPROC),
            error=function(e) {
                cat("  igphyml error:", conditionMessage(e), "\n")
                cat("  Falling back to pml for large clones...\n")
                getTrees(clones_large, build="pml", nproc=1)
            }
        )
        n_ok <- sum(sapply(result_large$trees, function(t) !is.null(t)))
        cat("  igphyml trees succeeded:", n_ok, "/", nrow(result_large), "\n")
        cat("  Finished:", format(Sys.time()), "\n\n")
    }
} else {
    cat("[6/7] No large clones — skipping igphyml\n\n")
}

# =============================================================================
# 7. formatClones + getTrees for SMALL clones with pml
# =============================================================================
result_small <- NULL
if (length(small_ids) > 0) {
    cat("[7/7] Building trees for SMALL clones with pml...\n")
    cat("  Clones:", length(small_ids), "\n")
    cat("  Started:", format(Sys.time()), "\n\n")

    clones_small <- suppressWarnings(
        formatClones(db_small,
                     seq="sequence_alignment", germ=germ_col,
                     clone="clone_id", v_call="v_call", j_call="j_call",
                     junc_len="junction_length", minseq=2, nproc=1)
    )
    cat("  Formatted:", nrow(clones_small), "clones\n")

    if (nrow(clones_small) > 0) {
        result_small <- getTrees(clones_small, build="pml", nproc=NPROC)
        n_ok <- sum(sapply(result_small$trees, function(t) !is.null(t)))
        cat("  pml trees succeeded:", n_ok, "/", nrow(result_small), "\n")
        cat("  Finished:", format(Sys.time()), "\n\n")
    }
} else {
    cat("[7/7] No small clones — skipping pml\n\n")
}

# =============================================================================
# Combine and save
# =============================================================================
cat("Combining and saving results...\n")

# Combine tibbles (both have same structure from getTrees)
result <- bind_rows(result_large, result_small)
n_success <- sum(sapply(result$trees, function(t) !is.null(t)))
cat("  Total trees:", n_success, "/", nrow(result), "\n\n")

# Save full result tibble
saveRDS(result, file.path(OUT_DIR, "all_trees.rds"))
cat("  ✓ RDS:", file.path(OUT_DIR, "all_trees.rds"), "\n")

# Newick files
newick_dir <- file.path(OUT_DIR, "newick_files")
dir.create(newick_dir, showWarnings=FALSE)
n_nwk <- 0
for (i in seq_len(nrow(result))) {
    tryCatch({
        t <- result$trees[[i]]
        if (!is.null(t)) {
            cid <- result$clone_id[i]
            nc  <- nrow(result$data[[i]])
            write.tree(t, file=file.path(newick_dir,
                                          paste0("clone_", cid, "_n", nc, ".nwk")))
            n_nwk <- n_nwk + 1
        }
    }, error=function(e) NULL)
}
cat("  ✓ Newick files:", n_nwk, "\n")

# Summary table — uses S4 airrClone @data slot directly
summary_rows <- lapply(seq_len(nrow(result)), function(i) {
    tryCatch({
        if (is.null(result$trees[[i]])) return(NULL)
        obj <- result$data[[i]]   # airrClone S4 object
        nc  <- nrow(obj@data)
        data.frame(
            clone_id  = as.character(result$clone_id[i]),
            n_cells   = nc,
            method    = ifelse(nc >= IGPHYML_MIN_SIZE, "igphyml", "pml"),
            n_samples = n_distinct(obj@data$sample_id),
            v_gene    = obj@v_gene,
            j_gene    = obj@j_gene,
            junc_len  = obj@junc_len,
            stringsAsFactors=FALSE
        )
    }, error=function(e) NULL)
})
tree_summary <- do.call(rbind, Filter(Negate(is.null), summary_rows))
tree_summary <- tree_summary[order(-tree_summary$n_cells), ]

write.table(tree_summary,
            file.path(OUT_DIR, "tree_summary.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)
cat("  ✓ Summary:", file.path(OUT_DIR, "tree_summary.tsv"), "\n\n")

cat("  Top 15 clones:\n")
print(head(tree_summary, 15))

cat("\n=============================================\n")
cat(" TREE BUILDING COMPLETE\n")
cat("=============================================\n")
cat(" igphyml clones (>=", IGPHYML_MIN_SIZE, "cells):", 
    ifelse(is.null(result_large), 0,
           sum(sapply(result_large$trees, function(t) !is.null(t)))), "\n")
cat(" pml     clones (<",  IGPHYML_MIN_SIZE, "cells):",
    ifelse(is.null(result_small), 0,
           sum(sapply(result_small$trees, function(t) !is.null(t)))), "\n")
cat(" Total trees:", n_success, "\n")
cat(" Finished:", format(Sys.time()), "\n")
cat(" NEXT: bash 07_annotate_visualize.sh\n")
cat("=============================================\n")
