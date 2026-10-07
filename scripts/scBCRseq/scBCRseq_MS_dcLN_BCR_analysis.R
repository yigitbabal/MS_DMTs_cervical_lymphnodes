setwd("P:/grp_laakso/Project_repository/Yigit/MS_FNA_followup")

library(tidyverse)
library(Seurat)
library(data.table)
library(scRepertoire)

source("scripts/00_functions.R")

#seurat objects
seurat.fna <- read_rds("seurat_objects/FNA/seurat_fna_curated_with_subset.rds")
seurat.fna$patient <- str_split(seurat.fna$sample, pattern = "_", simplify = T)[,1]

bcr.files <- list.files(path = "rawdata/BCR/FNA/", pattern = "filtered_contig_annotations.csv", full.names = T)
bcr.samples <- list.files(path = "rawdata/BCR/FNA/", pattern = "filtered_contig_annotations.csv", full.names = F)
bcr.samples <- str_split_fixed(bcr.samples, pattern = "_filtered_contig_annotations.csv", n = 2)[,1]

bcr.files <- bcr.files[!grepl(pattern = "MS001", x = bcr.files)]
bcr.samples <- bcr.samples[!grepl(pattern = "MS001", x = bcr.samples)]


bcr.list <- lapply(bcr.files, fread)
names(bcr.list) <- bcr.samples
bcr.samples <- gsub(pattern = "-", replacement = "_", x = bcr.samples)

## add cellid column in bcr.list as tcr.samples_barcode
bcr.list <- lapply(1:25, function(x){
  bcr.list[[x]]$cellid <- paste0(bcr.samples[x], "_", bcr.list[[x]]$barcode)
  return(bcr.list[[x]])
})

sample.meta <- fread("sample_meta_data.csv")
sample.meta <- sample.meta[match(bcr.samples, sample.meta$ID),]


bcr.list <- lapply(bcr.list, function(x){
  x = x[x$cellid %in% colnames(seurat.fna),]
  ## add cell type info from seurat object
  x$celltype <- seurat.fna$Celltypes_curated_detailed_withsub[match(x$cellid, colnames(seurat.fna))]
  ## add dmt info
  x$dmt <- seurat.fna$DMT[match(x$cellid, colnames(seurat.fna))]
  return(x)
})
names(bcr.list) <- bcr.samples


combined_BCR <- combineBCR(bcr.list, 
                           samples = sample.meta$ID, 
                           removeNA = T, removeMulti = T, 
                           filterMulti = F, filterNonproductive = T)


combined_BCR <- lapply(combined_BCR, function(x){
  x$celltype <- seurat.fna$Celltypes_curated_detailed_withsub[match(x$barcode, colnames(seurat.fna))]
  x$dmt <- seurat.fna$DMT[match(x$barcode, colnames(seurat.fna))]
  return(x)
})


celltypes <- c("DN","GC B cells (dark zone)", "Plasmablasts",
             "Tbet+CD11c+","GC B cells (light zone)",
             "USM", "Pre-B", "Naive B", "SM", "IFN-stimulated B")

# create combine_TCR.celltype list for each celltype
combined_BCR.celltype <- lapply(celltypes, function(x){
  combined_BCR.celltype <- lapply(combined_BCR, function(y){
    y <- y[y$celltype == x,]
    return(y)
  })
})

names(combined_BCR.celltype) <- celltypes

# 1. Initialize an empty list to hold your new cell-type specific lists
bcr_by_celltype <- list()

# 2. Loop through each cell type
for (ct in celltypes) {
  
  # 3. Use lapply to subset each sample's data frame within the combined_TCR list
  ct_specific_list <- lapply(combined_BCR, function(sample_df) {
    # Keep only the rows where the celltype column matches the current cell type
    subset_df <- sample_df[sample_df$celltype == ct, , drop = FALSE]
    return(subset_df)
  })
  
  # 4. Optional but highly recommended: 
  # Remove any samples (list elements) that ended up with 0 cells for this cell type
  # scRepertoire functions can sometimes throw errors on empty data frames
  ct_specific_list <- ct_specific_list[sapply(ct_specific_list, nrow) > 0]
  
  # 5. Save the resulting list into the master list, named after the cell type
  bcr_by_celltype[[ct]] <- ct_specific_list
}


fna.unique.clonetype <- lapply(bcr_by_celltype, function(x){
  p = clonalQuant(x, 
                  cloneCall="CTstrict", 
                  chain = "both",  
                  scale = T, exportTable = T)
  
  p$dmt = sample.meta$DMT[match(p$values, sample.meta$ID)]
  return(p)
})

write_rds(fna.unique.clonetype, file = "results/BCR/fna.unique.clonetype.rds")


plot.unique.clonetype <- lapply(1:length(fna.unique.clonetype), function(x){
  p = fna.unique.clonetype[[x]]
  p$dmt = factor(p$dmt, levels = c("Control", "BeforeTreatment", "antiCD20", "DMF", "NTZ"))
  
  my_comparisons <- list(
    c("BeforeTreatment", "Control"),
    c("BeforeTreatment", "antiCD20"),
    c("BeforeTreatment", "DMF")
  )
  
  # Filter to only comparisons where both groups have >= 3 observations
  valid_comparisons <- Filter(function(comp) {
    n1 = sum(p$dmt == comp[1], na.rm = TRUE)
    n2 = sum(p$dmt == comp[2], na.rm = TRUE)
    n1 >= 2 & n2 >= 2
  }, my_comparisons)
  
  plot = ggplot(p, aes(x = dmt, y = scaled, fill = dmt)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 3, height = 0) +
    labs(
      title = names(fna.unique.clonetype)[x],
      x = "DMT",
      y = "% of Unique BCR Clones in sample"
    ) +
    theme_bw() +
    theme(
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      legend.position = "none"
    )
  
  # Only add significance brackets if valid comparisons exist
  if (length(valid_comparisons) > 0) {
    plot = plot + stat_compare_means(
      comparisons = valid_comparisons,
      method = "t.test",
      paired = FALSE,
      label = "p.signif",
      step.increase = 0.1,
      color = "black"
    )
  }
  
  return(plot)
})

names(plot.unique.clonetype) <- names(fna.unique.clonetype)
patchwork::wrap_plots(plot.unique.clonetype)


plot.unique.clonetype <- lapply(1:length(fna.unique.clonetype), function(x){
  p = fna.unique.clonetype[[x]]
  p <- filter(p, dmt %in% c("Control", "BeforeTreatment", "antiCD20"))
  p$dmt = factor(p$dmt, levels = c("Control", "BeforeTreatment", "antiCD20"))
  
  my_comparisons <- list(
    c("BeforeTreatment", "Control"),
    c("BeforeTreatment", "antiCD20")
  )
  
  # Filter to only comparisons where both groups have >= 3 observations
  valid_comparisons <- Filter(function(comp) {
    n1 = sum(p$dmt == comp[1], na.rm = TRUE)
    n2 = sum(p$dmt == comp[2], na.rm = TRUE)
    n1 >= 2 & n2 >= 2
  }, my_comparisons)
  
  plot = ggplot(p, aes(x = dmt, y = scaled, fill = dmt)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 3, height = 0) +
    labs(
      title = names(fna.unique.clonetype)[x],
      x = "DMT",
      y = "% of Unique BCR Clones in sample"
    ) +
    theme_bw() +
    theme(
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      legend.position = "none"
    )
  
  # Only add significance brackets if valid comparisons exist
  if (length(valid_comparisons) > 0) {
    plot = plot + stat_compare_means(
      comparisons = valid_comparisons,
      method = "t.test",
      paired = FALSE,
      label = "p.signif",
      step.increase = 0.1,
      color = "black"
    )
  }
  
  return(plot)
})

names(plot.unique.clonetype) <- names(fna.unique.clonetype)
patchwork::wrap_plots(plot.unique.clonetype)


## for clone diversity
fna.diversity <- lapply(bcr_by_celltype, function(x){
  p = clonalDiversity(x, 
                      cloneCall = "CTstrict", 
                      metric = "norm.entropy",
                      skip.boots = T, 
                      exportTable = T)
  
  p$dmt = sample.meta$DMT[match(p$Group, sample.meta$ID)]
  return(p)
})

write_rds(fna.diversity, file = "results/BCR/fna.diversity.rds")
fna.diversity <- read_rds("results/BCR/fna.diversity.rds")

plot.diversity <- lapply(1:length(fna.diversity), function(x){
  p = fna.diversity[[x]]
  p <- filter(p, dmt %in% c("Control", "BeforeTreatment", "antiCD20"))
  p$dmt = factor(p$dmt, levels = c("Control", "BeforeTreatment", "antiCD20"))
  
  p$dmt_detail = sample.meta$DMT_detail[match(p$Group, sample.meta$ID)]
  p$relapse = "Non-relapsing"
  p$patient = gsub(pattern = "_FU", replacement = "", p$Group)
  p[p$patient == "MS022", "relapse"] = "Relapse-during-DMT"
  
  my_comparisons <- list(
    c("BeforeTreatment", "Control"),
    c("BeforeTreatment", "antiCD20")
  )
  
  # Filter to only comparisons where both groups have >= 3 observations
  valid_comparisons <- Filter(function(comp) {
    n1 = sum(p$dmt == comp[1], na.rm = TRUE)
    n2 = sum(p$dmt == comp[2], na.rm = TRUE)
    n1 >= 2 & n2 >= 2
  }, my_comparisons)
  
  # same seed/width used for both layers so line endpoints match the dots exactly
  jitter_pos <- position_jitter(width = 0.15, height = 0, seed = 42)
  
  plot = ggplot(p, aes(x = dmt, y = value)) +
    geom_boxplot(aes(fill = dmt), outlier.shape = NA) +
    geom_line(aes(group = patient), position = jitter_pos,
              color = "grey50", alpha = 0.6, linewidth = 0.4) +
    geom_point(aes(color = dmt_detail, shape = relapse),
               position = jitter_pos, size = 3) +
    scale_color_manual(values = c(
      "Control" = "#8C8C8C",
      "BeforeTreatment" = "#4D4D4D",
      "antiCD20 (rituximab)" = "#0073C2",
      "antiCD20 (ofatumumab)" = "#EFC000"
    )) +
    guides(fill = "none") +
    labs(
      title = names(fna.diversity)[x],
      x = "DMT",
      y = "Normalized Shannon diversity"
    ) +
    theme_bw() +
    theme(
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    )
  
  # Only add significance brackets if valid comparisons exist
  if (length(valid_comparisons) > 0) {
    plot = plot + stat_compare_means(
      comparisons = valid_comparisons,
      method = "t.test",
      paired = FALSE,
      label = "p.signif",
      step.increase = 0.1,
      color = "black"
    )
  }
  
  return(plot)
})
names(plot.diversity) <- names(fna.diversity)
patchwork::wrap_plots(plot.diversity)



















### create meta.data for germline tree
seurat.fna <- read_rds("seurat_objects/FNA/seurat_fna_curated_with_subset.rds")


seurat.fna.germline <- subset(seurat.fna, subset = Celltypes_curated_detailed_withsub %in% c("DN","GC B cells (dark zone)", "Plasmablasts",
                                                                                             "Tbet+CD11c+","GC B cells (light zone)",
                                                                                             "USM", "Pre-B", "Naive B", "SM", "IFN-stimulated B"))

unique(seurat.fna.germline$Celltypes_curated_detailed_withsub)

seurat_meta <- data.frame(
  barcode   = colnames(seurat.fna.germline),
  cell_type = seurat.fna.germline$Celltypes_curated_detailed_withsub,
  cluster_id= seurat.fna.germline$seurat_clusters,
  tissue    = "FNA"
)

seurat_meta$barcode <- gsub(pattern = "_FU", replacement = "-FU", seurat_meta$barcode)

write.csv(seurat_meta, 
          "bcr_germline_project/seurat_barcodes_annotated.csv",
          row.names = FALSE)
