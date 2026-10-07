setwd("P:/grp_laakso/Project_repository/Yigit/MS_FNA_followup")

library(tidyverse)
library(Seurat)
library(data.table)
library(scRepertoire)

source("scripts/00_functions.R")


clone_germline <- fread("bcr_germline_project/04_clones/cloned_germlined_germ-pass.tsv")

seurat.fna <- read_rds("seurat_objects/FNA/seurat_fna_curated_with_subset.rds")
seurat.fna$patient <- str_split(seurat.fna$sample, pattern = "_", simplify = T)[,1]

b.cells <- c("DN","GC B cells (dark zone)", "Plasmablasts",
             "Tbet+CD11c+","GC B cells (light zone)",
             "USM", "Pre-B", "Naive B", "SM", "IFN-stimulated B")

seurat.fna.bcell <- subset(seurat.fna, subset = Celltypes_curated_detailed_withsub %in% b.cells)

clone_germline$cell_id <- gsub(pattern = "-FU", replacement = "_FU", clone_germline$cell_id)

clone.table <- as.data.frame(table(clone_germline$clone_id))
clone.table <- clone.table[clone.table$Freq > 50,]

clone_germline.filtered <- clone_germline[clone_germline$clone_id %in% clone.table$Var1,]

seurat.fna.bcell$clone_id <- clone_germline.filtered$clone_id[match(colnames(seurat.fna.bcell), clone_germline.filtered$cell_id)]

DefaultAssay(seurat.fna.bcell) <- "RNA"
seurat.fna.bcell[["RNA"]] <- split(seurat.fna.bcell[["RNA"]], f = seurat.fna.bcell$orig.ident)
seurat.fna.bcell
#seurat.fna.bcell <- NormalizeData(seurat.fna.bcell, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.fna.bcell <- FindVariableFeatures(seurat.fna.bcell, nfeatures = 2000)
seurat.fna.bcell <- ScaleData(seurat.fna.bcell, features = VariableFeatures(seurat.fna.bcell))
seurat.fna.bcell <- RunPCA(seurat.fna.bcell, features = VariableFeatures(seurat.fna.bcell))
seurat.fna.bcell <- IntegrateLayers(
  object = seurat.fna.bcell, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = TRUE, k.weight = 62
)
seurat.fna.bcell <- FindNeighbors(seurat.fna.bcell, reduction = "integrated.rpca",  dims = 1:50)
#seurat.fna.bcell <- FindClusters(seurat.fna.bcell, resolution = 0.7, cluster.name = "bcell_cluster_0.7")
seurat.fna.bcell <- RunUMAP(seurat.fna.bcell, reduction = "integrated.rpca", dims = 1:50)
#seurat.fna.bcell <- RunTSNE(seurat.fna.bcell, reduction = "integrated.rpca", dims = 1:50)

DimPlot(seurat.fna.bcell, group.by = "Celltypes_curated_detailed_withsub")
seurat.fna.bcell$clone_id <- as.character(seurat.fna.bcell$clone_id)
cluster_states(seurat.fna.bcell, "Celltypes_curated_detailed_withsub","clone_id")

library(SCP)

seurat.fna.bcell$clone_id_patient <- paste(seurat.fna.bcell$patient, seurat.fna.bcell$clone_id, sep = "_")


p.list <- lapply(b.cells, function(x){
  seurat.fna.bcell.stat <- subset(seurat.fna.bcell, subset = Celltypes_curated_detailed_withsub == x & DMT %in% c("BeforeTreatment", "antiCD20") & !is.na(clone_id))
  seurat.fna.bcell.stat$DMT <- factor(seurat.fna.bcell.stat$DMT, levels = c(c("BeforeTreatment", "antiCD20")))
  p <- CellStatPlot(seurat.fna.bcell.stat, stat.by = "patient", group.by = "DMT", plot_type = "trend", title = paste0(x, " - BCR clonal share"))
  return(p)
})


seurat.fna.bcell.stat <- subset(seurat.fna.bcell, subset = Celltypes_curated_detailed_withsub == "SM" & DMT %in% c("BeforeTreatment", "antiCD20") & !is.na(clone_id))
seurat.fna.bcell.stat$DMT <- factor(seurat.fna.bcell.stat$DMT, levels = c(c("BeforeTreatment", "antiCD20")))
p.sm <- CellStatPlot(seurat.fna.bcell.stat, stat.by = "clone_id", group.by = "DMT", split.by = "patient", plot_type = "trend", title = "SM - BCR clonal share")
p.sm



# --- build ONE global palette, before calling the function at all ---

library(RColorBrewer)
pal_qual <- c(brewer.pal(12, "Paired"), brewer.pal(8, "Set2"))  # 20 distinct colors, you need 14
clone_palette <- setNames(pal_qual[seq_along(all_clone_ids)], all_clone_ids)

bcr_alluvial <- function(seurat, celltype, conditions, clone_palette){
  s1 <- subset(seurat, subset = Celltypes_curated_detailed_withsub == celltype & DMT %in% conditions & !is.na(clone_id))
  meta.bcell <- s1@meta.data
  
  df <- meta.bcell %>%
    mutate(clone_id = as.character(clone_id)) %>%
    dplyr::count(patient, DMT, clone_id, name = "n") %>%
    group_by(patient, DMT) %>%
    mutate(Percentage = 100 * n / sum(n)) %>%
    ungroup()
  
  df2 <- df %>%
    mutate(
      DMT = factor(DMT, levels = conditions),
      x_group = interaction(patient, DMT, sep = "_", drop = TRUE),
      x_group = factor(x_group, levels = interaction(
        rep(sort(unique(patient)), each = length(conditions)),
        rep(conditions, times = length(unique(patient))),
        sep = "_")),
      flow_id = interaction(patient, clone_id, drop = TRUE),
      clone_id = factor(clone_id, levels = names(clone_palette))  # <- locks in the global order
    )
  
  plot <- ggplot(df2, aes(x = x_group, y = Percentage,
                          alluvium = flow_id, stratum = flow_id, fill = clone_id)) +
    geom_flow(stat = "alluvium", color = "grey30", linewidth = 0.2, alpha = 0.85) +
    geom_stratum(color = "black", linewidth = 0.4) +
    scale_fill_manual(values = clone_palette, limits = names(clone_palette),
                      drop = TRUE, na.translate = FALSE) +
    labs(title = paste0(celltype, " - BCR clonal share"), x = NULL, y = "Percentage") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "right")
  
  return(plot)
}

b.cells.list <- c("DN","GC B cells (dark zone)", "Plasmablasts",
             "Tbet+CD11c+","GC B cells (light zone)",
             "USM", "Naive B", "SM", "IFN-stimulated B")

plot.alluvial.list <- lapply(b.cells.list, function(x){
  plot = bcr_alluvial(seurat.fna.bcell, x, c("BeforeTreatment", "antiCD20"), clone_palette)
  return(plot)
})

patchwork::wrap_plots(plot.alluvial.list)


library(ggplot2)
library(ggalluvial)
library(dplyr)

seurat.fna.bcell.stat <- subset(seurat.fna.bcell, subset = Celltypes_curated_detailed_withsub == "SM" & DMT %in% c("BeforeTreatment", "antiCD20") & !is.na(clone_id))
meta.bcell <- seurat.fna.bcell.stat@meta.data

df <- meta.bcell %>%
  mutate(clone_id = as.character(clone_id)) %>%
  dplyr::count(patient, DMT, clone_id, name = "n") %>%
  group_by(patient, DMT) %>%
  mutate(Percentage = 100 * n / sum(n)) %>%
  ungroup()

df2 <- df %>%
  mutate(
    DMT = factor(DMT, levels = c("BeforeTreatment", "antiCD20")),
    x_group = interaction(patient, DMT, sep = "_", drop = TRUE),
    x_group = factor(x_group, levels = interaction(
      rep(sort(unique(patient)), each = 2),
      rep(c("BeforeTreatment", "antiCD20"), times = length(unique(patient))),
      sep = "_")),
    flow_id = interaction(patient, clone_id, drop = TRUE)
  )

ggplot(df2, aes(x = x_group, y = Percentage,
                alluvium = flow_id, stratum = flow_id, fill = clone_id)) +
  geom_flow(stat = "alluvium", color = "grey30", linewidth = 0.2, alpha = 0.85) +
  geom_stratum(color = "black", linewidth = 0.4) +
  labs(title = "SM - BCR clonal share", x = NULL, y = "Percentage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "right")





unique(seurat.fna.bcell$Celltypes_curated_detailed_withsub)

seurat.fna.bcell.stat <- subset(seurat.fna.bcell, subset = Celltypes_curated_detailed_withsub == "DN" & DMT %in% c("BeforeTreatment", "antiCD20") & !is.na(clone_id))
seurat.fna.bcell.stat$DMT <- factor(seurat.fna.bcell.stat$DMT, levels = c(c("BeforeTreatment", "antiCD20")))
p.dn <- CellStatPlot(seurat.fna.bcell.stat, stat.by = "patient", group.by = "DMT", plot_type = "trend", title = "DN - BCR clonal share")

seurat.fna.bcell.stat <- subset(seurat.fna.bcell, subset = Celltypes_curated_detailed_withsub == "SM" & DMT %in% c("BeforeTreatment", "antiCD20") & !is.na(clone_id))
seurat.fna.bcell.stat$DMT <- factor(seurat.fna.bcell.stat$DMT, levels = c(c("BeforeTreatment", "antiCD20")))
p.sm <- CellStatPlot(seurat.fna.bcell.stat, stat.by = "patient", group.by = "DMT", plot_type = "trend", title = "SM - BCR clonal share")

seurat.fna.bcell.stat <- subset(seurat.fna.bcell, subset = Celltypes_curated_detailed_withsub == "Naive B" & DMT %in% c("BeforeTreatment", "antiCD20") & !is.na(clone_id))
seurat.fna.bcell.stat$DMT <- factor(seurat.fna.bcell.stat$DMT, levels = c(c("BeforeTreatment", "antiCD20")))
p.naive <- CellStatPlot(seurat.fna.bcell.stat, stat.by = "patient", group.by = "DMT", plot_type = "trend", title = "Naive B - BCR clonal share")

seurat.fna.bcell.stat <- subset(seurat.fna.bcell, subset = Celltypes_curated_detailed_withsub == "USM" & DMT %in% c("BeforeTreatment", "antiCD20") & !is.na(clone_id))
seurat.fna.bcell.stat$DMT <- factor(seurat.fna.bcell.stat$DMT, levels = c(c("BeforeTreatment", "antiCD20")))
p.usm <- CellStatPlot(seurat.fna.bcell.stat, stat.by = "patient", group.by = "DMT", plot_type = "trend", title = "USM - BCR clonal share")

seurat.fna.bcell.stat <- subset(seurat.fna.bcell, subset = Celltypes_curated_detailed_withsub == "GC B cells (light zone)" & DMT %in% c("BeforeTreatment", "antiCD20") & !is.na(clone_id))
seurat.fna.bcell.stat$DMT <- factor(seurat.fna.bcell.stat$DMT, levels = c(c("BeforeTreatment", "antiCD20")))
p.gclz <- CellStatPlot(seurat.fna.bcell.stat, stat.by = "patient", group.by = "DMT", plot_type = "trend", title = "GC B cells (light zone) - BCR clonal share")

seurat.fna.bcell.stat <- subset(seurat.fna.bcell, subset = Celltypes_curated_detailed_withsub == "IFN-stimulated B" & DMT %in% c("BeforeTreatment", "antiCD20") & !is.na(clone_id))
seurat.fna.bcell.stat$DMT <- factor(seurat.fna.bcell.stat$DMT, levels = c(c("BeforeTreatment", "antiCD20")))
p.ifnb <- CellStatPlot(seurat.fna.bcell.stat, stat.by = "patient", group.by = "DMT", plot_type = "trend", title = "IFN-stimulated B - BCR clonal share")

seurat.fna.bcell.stat <- subset(seurat.fna.bcell, subset = Celltypes_curated_detailed_withsub == "Plasmablasts" & DMT %in% c("BeforeTreatment", "antiCD20") & !is.na(clone_id))
seurat.fna.bcell.stat$DMT <- factor(seurat.fna.bcell.stat$DMT, levels = c(c("BeforeTreatment", "antiCD20")))
p.plasma <- CellStatPlot(seurat.fna.bcell.stat, stat.by = "patient", group.by = "DMT", plot_type = "trend", title = "Plasmablasts - BCR clonal share")

seurat.fna.bcell.stat <- subset(seurat.fna.bcell, subset = Celltypes_curated_detailed_withsub == "GC B cells (dark zone)" & DMT %in% c("BeforeTreatment", "antiCD20") & !is.na(clone_id))
seurat.fna.bcell.stat$DMT <- factor(seurat.fna.bcell.stat$DMT, levels = c(c("BeforeTreatment", "antiCD20")))
p.gcdz <- CellStatPlot(seurat.fna.bcell.stat, stat.by = "patient", group.by = "DMT", plot_type = "trend", title = "GC B cells (dark zone) - BCR clonal share")

seurat.fna.bcell.stat <- subset(seurat.fna.bcell, subset = Celltypes_curated_detailed_withsub == "Tbet+CD11c+" & DMT %in% c("BeforeTreatment", "antiCD20") & !is.na(clone_id))
seurat.fna.bcell.stat$DMT <- factor(seurat.fna.bcell.stat$DMT, levels = c(c("BeforeTreatment", "antiCD20")))
p.tbet <- CellStatPlot(seurat.fna.bcell.stat, stat.by = "patient", group.by = "DMT", plot_type = "trend", title = "Tbet+CD11c+ - BCR clonal share")

patchwork::wrap_plots(list(p.dn, p.sm, p.naive, p.usm, p.gclz, p.gcdz, p.ifnb, p.tbet, p.plasma), ncol = 3)




library(Seurat)
library(circlize)
library(dplyr)

# Extract metadata
meta <- seurat.fna.bcell@meta.data


meta.clones <- meta[!is.na(meta$clone_id),]

unique(meta.clones[meta.clones$Celltypes_curated_detailed_withsub == "Naive B" & meta.clones$DMT == "BeforeTreatment", "sample"])
unique(meta.clones[meta.clones$Celltypes_curated_detailed_withsub == "Naive B" & meta.clones$DMT == "antiCD20", "sample"])

unique(meta.clones[meta.clones$Celltypes_curated_detailed_withsub == "DN" & meta.clones$DMT == "BeforeTreatment", "sample"])
unique(meta.clones[meta.clones$Celltypes_curated_detailed_withsub == "DN" & meta.clones$DMT == "antiCD20", "sample"])

unique(meta.clones[meta.clones$Celltypes_curated_detailed_withsub == "SM" & meta.clones$DMT == "BeforeTreatment", "sample"])
unique(meta.clones[meta.clones$Celltypes_curated_detailed_withsub == "SM" & meta.clones$DMT == "antiCD20", "sample"])



## DN check
unique(meta[meta$Celltypes_curated_detailed_withsub == "DN" & meta$DMT == "BeforeTreatment", "clone_id"])
unique(meta[meta$Celltypes_curated_detailed_withsub == "DN" & meta$DMT == "antiCD20", "clone_id"])

intersect(unique(meta[meta$Celltypes_curated_detailed_withsub == "DN" & meta$DMT == "BeforeTreatment", "clone_id"]),
          unique(meta[meta$Celltypes_curated_detailed_withsub == "DN" & meta$DMT == "antiCD20", "clone_id"]))

meta[meta$Celltypes_curated_detailed_withsub == "DN" & meta$DMT == "BeforeTreatment" & meta$clone_id %in% c("10798", "4934") , "cellid"]
meta[meta$Celltypes_curated_detailed_withsub == "DN" & meta$DMT == "antiCD20" & meta$clone_id %in% c("10798", "4934") , "cellid"]



meta <- seurat.fna.bcell@meta.data
#meta <- meta[meta$DMT == "Control",]
meta <- meta[meta$DMT == "BeforeTreatment",]
meta <- meta[meta$DMT == "antiCD20",]

meta <- meta[meta$DMT == "BeforeTreatment" & meta$patient == "MS022",]
meta <- meta[meta$DMT == "antiCD20" & meta$patient == "MS022",]

meta <- meta[meta$DMT %in% c("BeforeTreatment","antiCD20"),]
meta <- meta[meta$patient == "MS022",]

meta$celltype_dmt <- paste(meta$Celltypes_curated_detailed_withsub, meta$DMT, sep = "-")



# Pull columns out as plain vectors using [[ ]]
clone_vec    <- as.character(unlist(meta[["clone_id"]]))
#celltype_vec <- as.character(unlist(meta[["Celltypes_curated_detailed_withsub"]]))
celltype_vec <- as.character(unlist(meta[["celltype_dmt"]]))


# Build a clean plain data.frame (not tibble)
meta_clean <- data.frame(
  clone_id   = clone_vec,
  celltype   = celltype_vec,
  stringsAsFactors = FALSE
)

# Filter
meta_filtered <- meta_clean[
  !is.na(meta_clean$clone_id) &
    meta_clean$clone_id != ""   &
    meta_clean$clone_id != "None", 
]

# Check it worked
str(meta_filtered)
head(meta_filtered)

# Base R equivalent of count()
counts <- as.data.frame(table(
  clone_id = meta_filtered$clone_id,
  celltype = meta_filtered$celltype
))

counts <- counts[counts$Freq > 0, ]
head(counts)


# Build matrix
mat <- table(meta_filtered$clone_id, meta_filtered$celltype)
mat <- as.matrix(mat)

# Cross-product for celltype-sharing view
mat_shared <- t(mat) %*% mat
diag(mat_shared) <- 0

# Plot
library(circlize)
library(RColorBrewer)

celltypes <- rownames(mat_shared)

# High-contrast "Paired" + "Dark2" combined
n <- length(celltypes)
colors_hc <- setNames(
  colorRampPalette(c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#984EA3",
    "#A65628", "#F781BF", "#00CED1", "#FFD700", "#1B9E77",
    "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02"
  ))(n),
  celltypes
)

circos.clear()
chordDiagram(mat_shared,
             grid.col = colors_hc,
             transparency = 0.3,
             symmetric = TRUE,
             annotationTrack = "grid",
             preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim <- get.cell.meta.data("xlim")
  circos.text(mean(xlim), 1, get.cell.meta.data("sector.index"),
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
}, bg.border = NA)
circos.clear()

