setwd("P:/grp_laakso/Project_repository/Yigit/MS_FNA_followup")

library(tidyverse)
library(Seurat)
library(data.table)
library(scRepertoire)


#seurat objects
seurat.fna <- read_rds("seurat_objects/FNA/seurat_fna_curated_with_subset.rds")
seurat.pbmc <- read_rds("seurat_objects/PBMC/seurat_pbmc_curated_with_subset.rds")


# tcr files
tcr.files <- list.files(path = "rawdata/TCR/FNA/", pattern = "filtered_contig_annotations.csv", full.names = T)
tcr.samples <- list.files(path = "rawdata/TCR/FNA/", pattern = "filtered_contig_annotations.csv", full.names = F)
tcr.samples <- str_split_fixed(tcr.samples, pattern = "_filtered_contig_annotations.csv", n = 2)[,1]

tcr.list <- lapply(tcr.files, fread)
names(tcr.list) <- tcr.samples
tcr.samples <- gsub(pattern = "-", replacement = "_", x = tcr.samples)

## add cellid column in tcr.list as tcr.samples_barcode
tcr.list <- lapply(1:25, function(x){
  tcr.list[[x]]$cellid <- paste0(tcr.samples[x], "_", tcr.list[[x]]$barcode)
  return(tcr.list[[x]])
})
  

sample.meta <- fread("sample_meta_data.csv")
sample.meta <- sample.meta[match(tcr.samples, sample.meta$ID),]


tcr.list <- lapply(tcr.list, function(x){
  x = x[x$cellid %in% colnames(seurat.fna),]
  ## add cell type info from seurat object
  x$celltype <- seurat.fna$Celltypes[match(x$cellid, colnames(seurat.fna))]
  ## add dmt info
  x$dmt <- seurat.fna$DMT[match(x$cellid, colnames(seurat.fna))]
  return(x)
})
names(tcr.list) <- tcr.samples



combined_TCR <- combineTCR(tcr.list, 
                           samples = sample.meta$ID)


p.unique.clone.fna <- clonalQuant(combined_TCR, 
            cloneCall="CTstrict", 
            chain = "both", 
            scale = F)

p.unique.clone.fna + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


fna.all.tcr.percent <- clonalQuant(combined_TCR, 
            cloneCall="CTstrict", 
            chain = "both",  
            scale = TRUE, exportTable = T)

fna.all.tcr.percent$dmt <- sample.meta$DMT[match(fna.all.tcr.percent$values, sample.meta$ID)]
fna.all.tcr.percent$dmt = factor(fna.all.tcr.percent$dmt, levels = c("Control", "BeforeTreatment","antiCD20", "DMF", "NTZ"))

fna.all.tcr.percent %>%
  ggplot(aes(x = dmt, y = scaled, fill = dmt)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 3) +
  labs(title = "All TCR Clonotypes", x = "DMT", y = "Frequency of All Clonotypes") +
  ggsignif::geom_signif(comparisons = list(c("BeforeTreatment", "Control"), c("BeforeTreatment", "antiCD20"), c("BeforeTreatment", "DMF")), 
                        map_signif_level = TRUE, test = "wilcox.test", color = "black", step_increase = 0.1) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) + 
  theme(legend.position = "none")

fna.all.tcr.diversity <- clonalDiversity(combined_TCR, 
            cloneCall = "CTstrict", 
            metric = "norm.entropy", 
            exportTable = T)

fna.all.tcr.diversity$dmt <- sample.meta$DMT[match(fna.all.tcr.diversity$Group, sample.meta$ID)]
fna.all.tcr.diversity$dmt = factor(fna.all.tcr.diversity$dmt, levels = c("Control", "BeforeTreatment","antiCD20", "DMF", "NTZ"))

fna.all.tcr.diversity %>%
  ggplot(aes(x = dmt, y = value, fill = dmt)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 3) +
  labs(title = "All TCR Clonotypes", x = "DMT", y = "Normalized Shannon Diversity") +
  ggsignif::geom_signif(comparisons = list(c("BeforeTreatment", "Control"), c("BeforeTreatment", "antiCD20"), c("BeforeTreatment", "DMF")), 
                        map_signif_level = TRUE, test = "wilcox.test", color = "black", step_increase = 0.1) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) + 
  theme(legend.position = "none")




combined_TCR <- lapply(combined_TCR, function(x){
  x$celltype <- seurat.fna$Celltypes_curated_detailed_withsub[match(x$barcode, colnames(seurat.fna))]
  x$dmt <- seurat.fna$DMT[match(x$barcode, colnames(seurat.fna))]
  return(x)
})

celltypes <- c("Naive CD4", "Effector CD4", "Naive CD8", "Early activated CD8",
               "Tfr", "Tfh-like","Treg", "IFN-stimulated CD4", "GC Tfh", "Th17",
               "Intermediate/Primed naive CD4 T","Early activated CD4",
               "Memory CD8", "NKlike.CD8")

# create combine_TCR.celltype list for each celltype
combined_TCR.celltype <- lapply(celltypes, function(x){
  combined_TCR.celltype <- lapply(combined_TCR, function(y){
    y <- y[y$celltype == x,]
    return(y)
  })
})


t.cells <- c("Naive CD4", "Effector CD4", "Naive CD8", "Early activated CD8",
             "Tfr", "Tfh-like","Treg", "IFN-stimulated CD4", "GC Tfh", "Th17",
             "Intermediate/Primed naive CD4 T","Early activated CD4",
             "Memory CD8", "NKlike.CD8")


# 1. Initialize an empty list to hold your new cell-type specific lists
tcr_by_celltype <- list()

# 2. Loop through each cell type
for (ct in t.cells) {
  
  # 3. Use lapply to subset each sample's data frame within the combined_TCR list
  ct_specific_list <- lapply(combined_TCR, function(sample_df) {
    # Keep only the rows where the celltype column matches the current cell type
    subset_df <- sample_df[sample_df$celltype == ct, , drop = FALSE]
    return(subset_df)
  })
  
  # 4. Optional but highly recommended: 
  # Remove any samples (list elements) that ended up with 0 cells for this cell type
  # scRepertoire functions can sometimes throw errors on empty data frames
  ct_specific_list <- ct_specific_list[sapply(ct_specific_list, nrow) > 0]
  
  # 5. Save the resulting list into the master list, named after the cell type
  tcr_by_celltype[[ct]] <- ct_specific_list
}


fna.unique.clonetype <- lapply(tcr_by_celltype, function(x){
  p = clonalQuant(x, 
                  cloneCall="CTstrict", 
                  chain = "both",  
                  scale = TRUE, exportTable = T)
  
  p$dmt = sample.meta$DMT[match(p$values, sample.meta$ID)]
  return(p)
})

write_rds(fna.unique.clonetype, file = "results/TCR/fna.unique.clonetype.rds")

plot.unique.clonetype <- lapply(1:length(fna.unique.clonetype), function(x){
  p = fna.unique.clonetype[[x]]
  p <- p[p$dmt %in% c("Control", "BeforeTreatment","antiCD20"),]
  p$dmt = factor(p$dmt, levels = c("Control", "BeforeTreatment","antiCD20"))
  
  
  plot = ggplot(p, aes(x = dmt, y = scaled, fill = dmt)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 3, height = 0) +
    labs(title = names(fna.unique.clonetype)[x], x = "DMT", y = "% of unique TCR clones in sample") +
    ggsignif::geom_signif(comparisons = list(c("BeforeTreatment", "Control"), c("BeforeTreatment", "antiCD20")), 
                          map_signif_level = TRUE, test = "t.test", color = "black", step_increase = 0.1) +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) + 
    theme(legend.position = "none")
  
  return(plot)
})

names(plot.unique.clonetype) <- names(fna.unique.clonetype)


patchwork::wrap_plots(plot.unique.clonetype)


## for clone diversity
fna.diversity <- lapply(tcr_by_celltype, function(x){
  p = clonalDiversity(x, 
                      cloneCall = "CTstrict", 
                      metric = "norm.entropy", 
                      exportTable = T)
  
  p$dmt = sample.meta$DMT[match(p$Group, sample.meta$ID)]
  return(p)
})


plot.fna.diversity <- lapply(1:length(fna.diversity), function(x){
  p = fna.diversity[[x]]
  p <- p[p$dmt %in% c("Control", "BeforeTreatment","antiCD20"),]
  p$dmt = factor(p$dmt, levels = c("Control", "BeforeTreatment","antiCD20"))
  
  plot = ggplot(p, aes(x = dmt, y = value, fill = dmt)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 3, height = 0) +
    labs(title = names(fna.unique.clonetype)[x], x = "DMT", y = "Normalized Shannon Diversity") +
    ggsignif::geom_signif(comparisons = list(c("BeforeTreatment", "Control"), c("BeforeTreatment", "antiCD20")), 
                          map_signif_level = TRUE, test = "t.test", color = "black", step_increase = 0.1) +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) + 
    theme(legend.position = "none")
  
  return(plot)
})

names(plot.fna.diversity) <- names(fna.diversity)


patchwork::wrap_plots(plot.fna.diversity)


# save objects
write_rds(combined_TCR, file = "results/TCR/FNA_combined_TCR.rds")
write_rds(tcr_by_celltype, file = "results/TCR/FNA_combined_TCR_by_celltype.rds")
write_rds(fna.unique.clonetype, file = "results/TCR/FNA_unique_clonetype_by_celltype.rds")
write_rds(fna.diversity, file = "results/TCR/FNA_diversity_by_celltype.rds")


tcr_by_celltype <- read_rds("results/TCR/FNA_combined_TCR_by_celltype.rds")




# ============================================================
# TCR Clonal Expansion Analysis from scRepertoire combineTCR
# Structure: tcr_by_celltype[[celltype]][[sample]] = data.table
# ============================================================

library(dplyr)
library(tidyr)
library(purrr)

# -----------------------------------------------------------
# STEP 1: Define clonotype column to use
# Options: "CTstrict" (most specific, gene+nt), 
#          "CTgene"   (gene-level),
#          "CTaa"     (amino acid CDR3)
# CTstrict is recommended for expansion analysis
# -----------------------------------------------------------
CLONE_COL <- "CTstrict"

# -----------------------------------------------------------
# STEP 2: Core function — expansion stats per sample
# -----------------------------------------------------------
get_expansion_stats <- function(sample_dt, clone_col = CLONE_COL) {
  
  # Remove cells with no TCR (NA clonotype)
  df <- sample_dt[!is.na(get(clone_col)) & get(clone_col) != ""]
  
  total_cells    <- nrow(df)
  
  if (total_cells == 0) {
    return(data.frame(
      total_cells         = 0,
      n_unique_clones     = 0,
      n_expanded_clones   = 0,   # unique clone IDs that appear ≥2x
      n_cells_expanded    = 0,   # cells belonging to expanded clones
      pct_cells_expanded  = NA_real_
    ))
  }
  
  clone_counts <- df %>%
    dplyr::count(.data[[clone_col]], name = "clone_size")
  
  expanded_clones <- clone_counts %>% filter(clone_size >= 2)
  
  n_cells_expanded <- sum(expanded_clones$clone_size)
  
  data.frame(
    total_cells        = total_cells,
    n_unique_clones    = nrow(clone_counts),
    n_expanded_clones  = nrow(expanded_clones),       # ≥2 unique clones
    n_cells_expanded   = n_cells_expanded,             # cells in expanded clones
    pct_cells_expanded = (n_cells_expanded / total_cells) * 100
  )
}

# -----------------------------------------------------------
# STEP 3: Iterate over all cell types and samples
# -----------------------------------------------------------
expansion_df <- imap_dfr(tcr_by_celltype, function(sample_list, celltype) {
  
  imap_dfr(sample_list, function(sample_dt, sample_name) {
    
    stats <- get_expansion_stats(sample_dt)
    
    data.frame(
      cell_type  = celltype,
      sample     = sample_name,
      stats,
      stringsAsFactors = FALSE
    )
  })
})








####
CLONE_COL <- "CTstrict"

# -----------------------------------------------------------
# STEP 2: Core function — expansion stats per sample (per cell type)
# -----------------------------------------------------------
get_expansion_stats <- function(sample_dt, clone_col = CLONE_COL) {
  
  df <- sample_dt[!is.na(get(clone_col)) & get(clone_col) != ""]
  
  total_cells <- nrow(df)
  
  if (total_cells == 0) {
    return(data.frame(
      total_cells        = 0,
      n_unique_clones    = 0,
      n_expanded_clones  = 0,
      n_cells_expanded   = 0
    ))
  }
  
  clone_counts <- df %>%
    dplyr::count(.data[[clone_col]], name = "clone_size")
  
  expanded_clones <- clone_counts %>% filter(clone_size >= 2)
  
  n_cells_expanded <- sum(expanded_clones$clone_size)
  
  data.frame(
    total_cells       = total_cells,
    n_unique_clones   = nrow(clone_counts),
    n_expanded_clones = nrow(expanded_clones),
    n_cells_expanded  = n_cells_expanded
  )
}

# -----------------------------------------------------------
# STEP 3: Iterate over all cell types and samples
# -----------------------------------------------------------
expansion_df <- imap_dfr(tcr_by_celltype, function(sample_list, celltype) {
  imap_dfr(sample_list, function(sample_dt, sample_name) {
    stats <- get_expansion_stats(sample_dt)
    data.frame(
      cell_type = celltype,
      sample    = sample_name,
      stats,
      stringsAsFactors = FALSE
    )
  })
})

# -----------------------------------------------------------
# STEP 4: pct_cells_expanded = this cell type's expanded cells
#         ÷ total expanded cells across ALL cell types, per sample
# -----------------------------------------------------------
sample_totals <- expansion_df %>%
  group_by(sample) %>%
  summarise(total_expanded_in_sample = sum(n_cells_expanded), .groups = "drop")

expansion_df <- expansion_df %>%
  left_join(sample_totals, by = "sample") %>%
  mutate(pct_cells_expanded = ifelse(total_expanded_in_sample > 0,
                                     (n_cells_expanded / total_expanded_in_sample) * 100,
                                     NA_real_))

####













# -----------------------------------------------------------
# STEP 4: Add disease group annotation
# (adjust the pattern to match your sample naming)
# -----------------------------------------------------------
expansion_df <- expansion_df %>%
  mutate(
    group = case_when(
      grepl("^CTRL", sample) ~ "Control",
      grepl("^MS",   sample) ~ "MS",
      TRUE                   ~ "Unknown"
    )
  ) %>%
  select(cell_type, sample, group, everything())

expansion_df$group = sample.meta$DMT[match(expansion_df$sample, sample.meta$ID)]

# -----------------------------------------------------------
# STEP 5: Wide format — % expanded per cell type × sample
# (useful for heatmaps or group comparisons)
# -----------------------------------------------------------
expansion_wide <- expansion_df %>%
  select(cell_type, sample, group, pct_cells_expanded) %>%
  pivot_wider(
    names_from  = cell_type,
    values_from = pct_cells_expanded
  )

# -----------------------------------------------------------
# STEP 6: Summary table — mean ± SD per group per cell type
# -----------------------------------------------------------
expansion_summary <- expansion_df %>%
  group_by(cell_type, group) %>%
  summarise(
    n_samples          = n(),
    mean_pct_expanded  = mean(pct_cells_expanded, na.rm = TRUE),
    sd_pct_expanded    = sd(pct_cells_expanded,   na.rm = TRUE),
    median_pct_expanded = median(pct_cells_expanded, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------------------------------------
# STEP 7: Export results
# -----------------------------------------------------------
write.csv(expansion_df,      "results/TCR/TCR_expansion_per_sample_celltype.csv",  row.names = FALSE)
#write.csv(expansion_wide,    "results/TCR/TCR_expansion_wide_format.csv",           row.names = FALSE)
#write.csv(expansion_summary, "results/TCR/TCR_expansion_summary_group.csv",         row.names = FALSE)

# -----------------------------------------------------------
# STEP 8: Quick visualisation (ggplot2)
# -----------------------------------------------------------
library(ggplot2)

expansion_df <- read.csv("results/TCR/TCR_expansion_per_sample_celltype.csv")
expansion_df$group <- sample.meta$DMT[match(expansion_df$sample, sample.meta$ID)]
expansion_df$group <- factor(expansion_df$group , levels = c("Control", "BeforeTreatment", "antiCD20", "DMF", "NTZ"))
expansion_df$patient <- gsub(pattern = "_FU", "", expansion_df$sample)
expansion_df$dmt_detail <- sample.meta$DMT_detail[match(expansion_df$sample, sample.meta$ID)]
expansion_df$relapse <- "Non-relapsing"
expansion_df[expansion_df$patient == "MS022", "relapse"] <- "Relapse-during-DMT"

expansion_df %>%
  filter(cell_type == "NKlike.CD8") %>%
  filter(group %in% c("Control", "BeforeTreatment", "antiCD20")) %>%
  ggplot(aes(x = group, y = pct_cells_expanded)) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.6) +
  geom_line(aes(group = patient),
            color = "grey50", alpha = 0.6, linewidth = 0.4) +
  geom_point(aes(color = dmt_detail, shape = relapse),
              size = 2) +
  scale_color_manual(values = c(
    "Control" = "#8C8C8C",
    "BeforeTreatment" = "#4D4D4D",
    "antiCD20 (rituximab)" = "#0073C2",
    "antiCD20 (ofatumumab)" = "#EFC000"
  )) +
  labs(
    title    = "Expanded NKlike CD8 T cells",
    x        = "",
    y        = "% Expanded Cells",
    fill     = "Group", color = "DMT detail", shape = "Relapse status"
  ) +
  guides(fill = "none") +
  theme_bw(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))


# -----------------------------------------------------------
# STEP 1: Define parameters
# -----------------------------------------------------------
CLONE_COL    <- "CTstrict"
#TARGET_CT    <- "NKlike.CD8"   # must match name in tcr_by_celltype exactly
TARGET_CT <- c("Early activated CD8", "Memory CD8", "Naive CD8")
EXPAND_CUTOFF <- 2             # clone size threshold

# -----------------------------------------------------------
# STEP 2: Extract expanded clones + barcodes per sample
# -----------------------------------------------------------
#memory_cd8_list <- tcr_by_celltype[[TARGET_CT]]
memory_cd8_list <- tcr_by_celltype[TARGET_CT]

expanded_barcodes_list <- imap(memory_cd8_list, function(ct_list, ct_name) {
  
  imap(ct_list, function(sample_dt, sample_name) {
    
    df <- sample_dt[!is.na(get(CLONE_COL)) & get(CLONE_COL) != ""]
    
    if (nrow(df) == 0) {
      message("  [SKIP] ", ct_name, " / ", sample_name, " — no TCR-matched cells")
      return(NULL)
    }
    
    df <- df %>%
      group_by(.data[[CLONE_COL]]) %>%
      mutate(clone_size = n()) %>%
      ungroup()
    
    expanded <- df %>% filter(clone_size >= EXPAND_CUTOFF)
    
    n_total    <- nrow(df)
    n_expanded <- nrow(expanded)
    
    message(sprintf("  [%s / %s] total=%d | expanded cells=%d (%.1f%%)",
                    ct_name, sample_name, n_total, n_expanded,
                    100 * n_expanded / n_total))
    
    if (n_expanded == 0) return(NULL)
    
    expanded
  })
})

# Remove NULL (samples with no expanded cells)
expanded_barcodes_list <- Filter(Negate(is.null), expanded_barcodes_list)

cat("\nSamples with expanded Memory CD8 clones:", 
    length(expanded_barcodes_list), "\n")

#write_rds(expanded_barcodes_list, file = "results/TCR/tcrdist3_dawit/expanded_barcodes_list.rds")
write_rds(expanded_barcodes_list, file = "results/TCR/tcrdist3_dawit/expanded_barcodes_list_nkcd8.rds")


# Get all unique sample names across the three cell-type lists
all_samples <- unique(unlist(map(expanded_barcodes_list, names)))

# For each sample, pull its data from each cell type and stack them
combined_list <- map(all_samples, function(s) {
  dfs <- map(expanded_barcodes_list, s) %>% compact()   # grab sample `s` from each cell-type list, drop NULLs
  if (length(dfs) == 0) return(NULL)
  bind_rows(dfs, .id = "cell_type")
}) %>%
  set_names(all_samples) %>%
  compact()

write_rds(combined_list, file = "results/TCR/tcrdist3_dawit/expanded_barcodes_list_allCD8.rds")



# -----------------------------------------------------------
# STEP 3: Combined flat dataframe (all samples together)
# —— useful for overview & barcode extraction for Seurat
# -----------------------------------------------------------
expanded_all <- bind_rows(combined_list, .id = "sample_source")

write.csv(expanded_all, file = "results/TCR/allCD8_expanded_all.csv", row.names = F)
#write.csv(expanded_all, file = "results/TCR/NKlikecd8_expanded_all.csv", row.names = F)
#write.csv(expanded_all, file = "results/TCR/cd8_expanded_all.csv", row.names = F)

# All unique barcodes across samples (for Seurat subsetting)
expanded_barcodes_all <- expanded_all$barcode
#cat("Total expanded Memory CD8 cells:", length(expanded_barcodes_all), "\n")
#cat("Total expanded NK like CD8 cells:", length(expanded_barcodes_all), "\n")
cat("Total expanded CD8 cells:", length(expanded_barcodes_all), "\n")

# -----------------------------------------------------------
# STEP 4: Format each sample for VDJdb matching
# VDJdb expects: CDR3aa, V gene, J gene (alpha and/or beta)
# We work with beta chain (TRB) as primary — most VDJdb entries
# -----------------------------------------------------------

format_for_vdjdb <- function(sample_dt, sample_name) {
  
  df <- as.data.frame(sample_dt)
  
  # --- Parse CTgene to extract TRA / TRB gene info ---
  # CTgene format: "TRAV12-1.TRAJ6.TRAC_TRBV6-2..TRBJ1-1.TRBC1"
  # Split on "_" to separate alpha/beta
  df <- df %>%
    separate(CTgene, into = c("TRA_genes", "TRB_genes"), 
             sep = "_", extra = "merge", fill = "right", remove = FALSE) %>%
    # Parse alpha chain genes
    separate(TRA_genes, into = c("TRAV", "TRAJ", "TRAC"), 
             sep = "\\.", extra = "drop", fill = "right", remove = FALSE) %>%
    # Parse beta chain genes  
    separate(TRB_genes, into = c("TRBV", "TRBJ_raw", "TRBC"),
             sep = "\\.", extra = "drop", fill = "right", remove = FALSE) %>%
    # CTaa format: "CDR3aa_TRA_CDR3aa_TRB"
    separate(CTaa, into = c("CDR3aa_TRA", "CDR3aa_TRB"),
             sep = "_", extra = "merge", fill = "right", remove = FALSE) %>%
    # CTnt format: same structure as CTaa
    separate(CTnt, into = c("CDR3nt_TRA", "CDR3nt_TRB"),
             sep = "_", extra = "merge", fill = "right", remove = FALSE) %>%
    mutate(
      sample     = sample_name,
      # Standardize gene names (remove allele info if present, e.g. TRBV6-2*01 → TRBV6-2)
      TRBV_clean = sub("\\*.*", "", TRBV),
      TRAV_clean = sub("\\*.*", "", TRAV),
      TRAJ_clean = sub("\\*.*", "", TRAJ),
      TRBJ_clean = sub("\\*.*", "", TRBJ_raw)
    )
  
  # --- Build VDJdb-ready table (beta chain primary) ---
  vdjdb_ready <- df %>%
    select(
      barcode,
      sample,
      group     = dmt,          # Control / MS label already in TCR data
      clone_id  = all_of(CLONE_COL),
      clone_size,
      # Beta chain (most VDJdb entries)
      CDR3_beta  = CDR3aa_TRB,
      V_beta     = TRBV_clean,
      J_beta     = TRBJ_clean,
      CDR3nt_beta = CDR3nt_TRB,
      # Alpha chain (for paired analysis)
      CDR3_alpha  = CDR3aa_TRA,
      V_alpha     = TRAV_clean,
      J_alpha     = TRAJ_clean,
      CDR3nt_alpha = CDR3nt_TRA,
      # Full clonotype strings
      CTstrict, CTgene, CTaa, CTnt
    ) %>%
    filter(!is.na(CDR3_beta) & CDR3_beta != "" & CDR3_beta != "NA") %>%
    # Add species/chain columns required by VDJdb schema
    mutate(
      Species    = "HomoSapiens",
      Gene       = "TRB"
    )
  
  vdjdb_ready
}

#vdjdb_list <- imap(expanded_barcodes_list, format_for_vdjdb)
vdjdb_list <- imap(combined_list, format_for_vdjdb)

vdjdb_combined <- bind_rows(vdjdb_list)


clone_unique <- vdjdb_combined %>%
  group_by(CDR3_beta, V_beta, J_beta, Species, Gene) %>%
  summarise(
    n_cells      = n(),
    n_samples    = n_distinct(sample),
    samples      = paste(unique(sample), collapse = ";"),
    groups       = paste(unique(group),  collapse = ";"),
    clone_ids    = paste(unique(clone_id), collapse = ";"),
    CDR3_alpha   = paste(unique(na.omit(CDR3_alpha)), collapse = ";"),
    V_alpha      = paste(unique(na.omit(V_alpha)),    collapse = ";"),
    .groups = "drop"
  ) %>%
  arrange(desc(n_cells))

#write.csv(clone_unique, "results/TCR/NKlike_CD8_unique_clones_vdjdb_input.csv", row.names = FALSE)
write.csv(clone_unique, "results/TCR/all_CD8_unique_clones_vdjdb_input.csv", row.names = FALSE)






















































# persistent clones

paired.patients <- c("MS004", "MS005", "MS010", "MS018", "MS020", "MS021", "MS022")

persistent.clones <- lapply(paired.patients, function(x){
  res = find_persistent_clones(tcr_list = combined_TCR, 
                               pre_sample_name = x, 
                               post_sample_name = paste0(x, "_FU"))
  return(res)
})

names(persistent.clones) <- paired.patients


persistent.cells <- c(persistent.clones$MS004$persistent_barcodes, persistent.clones$MS005$persistent_barcodes, persistent.clones$MS010$persistent_barcodes, 
                      persistent.clones$MS018$persistent_barcodes, persistent.clones$MS020$persistent_barcodes, 
                      persistent.clones$MS021$persistent_barcodes, persistent.clones$MS022$persistent_barcodes)

t.cellid <- rownames(seurat.fna@meta.data[seurat.fna$Celltypes %in% t.cells,])
persistent.t.cells <- persistent.cells[persistent.cells %in% t.cellid]

seurat.fna$persistent_T_clone <- ifelse(colnames(seurat.fna) %in% persistent.t.cells, "Persistent Clone", NA)
unique(seurat.fna$persistent_T_clone)

DimPlot(seurat.fna, cells.highlight = persistent.t.cells)

persistent.df <- as.data.frame(table(seurat.fna$persistent_T_clone, seurat.fna$sample))
persistent.df <- persistent.df[persistent.df$Freq > 0,]
persistent.df$dmt <- sample.meta$DMT[match(persistent.df$Var2, sample.meta$ID)]
persistent.df$patient <- gsub(pattern = "_FU", replacement = "", x = persistent.df$Var2)


seurat.fna.tcells <- subset(seurat.fna, subset = Celltypes %in% t.cells)

t.cell.total.number <- as.data.frame(table(seurat.fna.tcells$sample))

persistent.df$total_t_cells <- t.cell.total.number$Freq[match(persistent.df$Var2, t.cell.total.number$Var1)]
persistent.df$percent_persistent <- persistent.df$Freq / persistent.df$total_t_cells * 100

## perform wilcox test with ggsignif for paired samples
persistent.df %>%
  filter(patient %in% c("MS010", "MS018", "MS021", "MS022")) %>%
  ggplot(aes(x = dmt, y = percent_persistent, fill = dmt)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 3) +
  geom_line(aes(group = patient), color = "black") +
  geom_text(aes(label = patient), vjust = -0.5, size = 3) +
  ggsignif::geom_signif(comparisons = list(c("BeforeTreatment", "antiCD20")), 
                        map_signif_level = TRUE, test = "wilcox.test", test.args = list(paired = T), color = "black", step_increase = 0.1) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) + 
  theme(legend.position = "none")

  
persistent.df %>%
  filter(patient %in% c("MS004", "MS005")) %>%
  ggplot(aes(x = dmt, y = percent_persistent, fill = dmt)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 3) +
  geom_line(aes(group = patient), color = "black") +
  geom_text(aes(label = patient), vjust = -0.5, size = 3) +
  ggsignif::geom_signif(comparisons = list(c("BeforeTreatment", "DMF")), 
                        map_signif_level = TRUE, test = "wilcox.test", test.args = list(paired = T), color = "black", step_increase = 0.1) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) + 
  theme(legend.position = "none")

