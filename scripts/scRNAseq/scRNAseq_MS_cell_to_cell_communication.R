setwd("P:/grp_laakso/Project_repository/Yigit/MS_FNA_followup")

library(tidyverse)
library(Seurat)
library(data.table)
library(CellChat)
library(future)


source("scripts/00_functions.R")

#seurat objects
seurat.fna <- read_rds("seurat_objects/FNA/seurat_fna_curated_with_subset.rds")
seurat.fna$patient <- str_split(seurat.fna$sample, pattern = "_", simplify = T)[,1]

seurat.pbmc <- read_rds("seurat_objects/PBMC/seurat_pbmc_curated_with_subset.rds")
seurat.pbmc$patient <- str_split(seurat.pbmc$sample, pattern = "-", simplify = T)[,1]
seurat.pbmc$patient <- gsub(pattern = "FU", replacement = "", x = seurat.pbmc$patient)


seurat.fna$samples <- as.factor(seurat.fna$sample)

unique(seurat.fna$Celltypes_curated_detailed_withsub)

plan("multisession", workers = 5)
options(future.globals.maxSize = 1000 * 10240^2)

fna_cellchat_list <- run_comparative_cellchat_lapply(
  seurat_obj = seurat.fna, 
  celltype_col = "Celltypes_curated_detailed_withsub",
  condition_col = "DMT"
)

# Merge the list into a single comparative object
cellchat.merged.fna <- mergeCellChat(fna_cellchat_list, add.names = names(fna_cellchat_list))
write_rds(cellchat.merged.fna, file = "results/cellchat/FNA_merged_cellchat.rds")
cellchat.merged.fna <- read_rds("results/cellchat/FNA_merged_cellchat.rds")

seurat.pbmc$samples <- as.factor(seurat.pbmc$sample)

pbmc_cellchat_list <- run_comparative_cellchat_lapply(
  seurat_obj = seurat.pbmc, 
  celltype_col = "Celltypes_curated_detailed_withsub",
  condition_col = "DMT"
)

cellchat.merged.pbmc <- mergeCellChat(pbmc_cellchat_list, add.names = names(pbmc_cellchat_list))
write_rds(cellchat.merged.pbmc, file = "results/cellchat/PBMC_merged_cellchat.rds")





## cellchat viz for FNA

unique(seurat.fna$Celltypes_curated_detailed_withsub)

b.cells <- c("DN","GC B cells (dark zone)", "Plasmablasts",
                "Tbet+CD11c+","GC B cells (light zone)",
                "USM", "Pre-B", "Naive B", "SM", "IFN-stimulated B")

t.cells <- c("Naive CD4", "Effector CD4", "Naive CD8", "Early activated CD8",
             "Tfr", "Tfh-like","Treg", "IFN-stimulated CD4", "GC Tfh", "Th17",
             "Intermediate/Primed naive CD4 T","Early activated CD4",
             "Memory CD8", "NKlike.CD8")

LR.fna.before.ctrl <- netVisual_bubble(cellchat.merged.fna, sources.use = t.cells, targets.use = b.cells,  comparison = c(2, 1), angle.x = 45, remove.isolate = F, return.data = T)
LR.fna.before.ctrl.data <- LR.fna.before.ctrl$communication

LR.fna.anticd20.before <- netVisual_bubble(cellchat.merged.fna, sources.use = t.cells, targets.use = b.cells,  comparison = c(4, 2), angle.x = 45, remove.isolate = F, return.data = T)
LR.fna.anticd20.before.data <- LR.fna.anticd20.before$communication

LR.fna.dmf.before <- netVisual_bubble(cellchat.merged.fna, sources.use = t.cells, targets.use = b.cells,  comparison = c(3, 2), angle.x = 45, remove.isolate = F, return.data = T)
LR.fna.dmf.before.data <- LR.fna.dmf.before$communication

LR.fna.ntz.before <- netVisual_bubble(cellchat.merged.fna, sources.use = t.cells, targets.use = b.cells,  comparison = c(5, 2), angle.x = 45, remove.isolate = F, return.data = T)
LR.fna.ntz.before.data <- LR.fna.ntz.before$communication

LR.fna.plots <- lapply(list(control = LR.fna.before.ctrl.data,
                         antiCD20 = LR.fna.anticd20.before.data,
                         DMF = LR.fna.dmf.before.data,
                         NTZ = LR.fna.ntz.before.data), function(x){
                           
                           plot = plot_split_LR_bubble(
                             LR.data = x, 
                             seurat_obj = seurat.fna, 
                             celltype_col = "Celltypes_curated_detailed_withsub", 
                             condition_col = "DMT"
                           )
                           
                           return(plot)
                         })


LR.fna.plots$antiCD20
ggsave(filename = "figures/cellchat/FNA_LR_antiCD20_update.pdf", width = 680, height = 780, units = "mm")
LR.fna.plots$control
ggsave(filename = "figures/cellchat/FNA_LR_beforetreatment_update.pdf", width = 680, height = 780, units = "mm")
LR.fna.plots$DMF
ggsave(filename = "figures/cellchat/FNA_LR_DMF_update.pdf", width = 680, height = 780, units = "mm")
LR.fna.plots$NTZ
ggsave(filename = "figures/cellchat/FNA_LR_NTZ_update.pdf", width = 680, height = 780, units = "mm")

#write_rds(LR.fna.plots, file = "results/cellchat/FNA_LR_bubble_plots.rds")


## cellchat viz for pbmc

unique(seurat.pbmc$Celltypes_curated_detailed_withsub)

b.cells <- c("SM", "USM", "Naive B", "Plasmablasts", "Tbet+CD11c+")

t.cells <- c("Memory CD8","Naive CD4","Vδ2+ Gamma Delta T Cell","Naive CD8",
             "Treg","Vδ1+ Gamma Delta T Cell","Th17","Effector CD4")

LR.pbmc.before.anticd20 <- netVisual_bubble(cellchat.merged.pbmc, sources.use = t.cells, targets.use = b.cells,  comparison = c(2, 1), angle.x = 45, remove.isolate = F, return.data = T)
LR.pbmc.before.anticd20.data <- LR.pbmc.before.anticd20$communication

LR.pbmc.before.ntz <- netVisual_bubble(cellchat.merged.pbmc, sources.use = t.cells, targets.use = b.cells,  comparison = c(3, 1), angle.x = 45, remove.isolate = F, return.data = T)
LR.pbmc.before.ntz.data <- LR.pbmc.before.ntz$communication

LR.pbmc.plots <- lapply(list(antiCD20 = LR.pbmc.before.anticd20.data,
                            NTZ = LR.pbmc.before.ntz.data), function(x){
                              
                              plot = plot_split_LR_bubble(
                                LR.data = x, 
                                seurat_obj = seurat.pbmc, 
                                celltype_col = "Celltypes_curated_detailed_withsub", 
                                condition_col = "DMT"
                              )
                              
                              return(plot)
                            })
LR.pbmc.plots$antiCD20
ggsave(filename = "figures/cellchat/PBMC_LR_antiCD20.pdf", width = 680, height = 780, units = "mm")
LR.pbmc.plots$NTZ
ggsave(filename = "figures/cellchat/PBMC_LR_ntz.pdf", width = 680, height = 780, units = "mm")

#write_rds(LR.pbmc.plots, file = "results/cellchat/PBMC_LR_bubble_plots.rds")




#### exhausted T cells and B cell interactions
seurat.fna$exhausted_included <- seurat.fna$Celltypes_curated_detailed_withsub
seurat.fna.cd8 <- read_rds("seurat_objects/seurat.fna.cd8.rds")
seurat.fna.cd8.exhaust <- read_rds("seurat_objects/seurat.fna.cd8.exhaust.rds")

seurat.fna.cd8.exhaust$subtypes <- seurat.fna.cd8$subtypes[match(colnames(seurat.fna.cd8.exhaust), names(seurat.fna.cd8$subtypes))]
unique(seurat.fna.cd8.exhaust$subtypes)



seurat.fna$exhausted_included[colnames(seurat.fna.cd8.exhaust)] <- as.character(seurat.fna.cd8.exhaust$subtypes)
unique(seurat.fna$exhausted_included)

seurat.fna.anticd20 <- subset(seurat.fna, subset = DMT %in% c("Control","BeforeTreatment", "antiCD20"))
unique(seurat.fna.anticd20$DMT)
unique(seurat.fna.anticd20$exhausted_included)

exhausted_cellchat_list <- run_comparative_cellchat_lapply(
  seurat_obj = seurat.fna.anticd20, 
  celltype_col = "exhausted_included",
  condition_col = "DMT"
)


cellchat.merged.exhausted <- mergeCellChat(exhausted_cellchat_list, add.names = names(exhausted_cellchat_list))
write_rds(cellchat.merged.exhausted, file = "results/cellchat/exhausted_merged_cellchat.rds")
cellchat.merged.exhausted <- read_rds("results/cellchat/exhausted_merged_cellchat.rds")

b.cells <- c("DN","GC B cells (dark zone)", "Plasmablasts",
             "Tbet+CD11c+","GC B cells (light zone)",
             "USM", "Pre-B", "Naive B", "SM", "IFN-stimulated B")

t.cells <- c("TRM-EX (c5)", "TRM-TEX (c6)", "TRM-TEX (c13)", "TEX (c15)")


LR.fna.before.ctrl <- netVisual_bubble(cellchat.merged.exhausted, sources.use = t.cells, targets.use = b.cells,  comparison = c(2, 1), angle.x = 45, remove.isolate = F, return.data = T)
LR.fna.before.ctrl.data <- LR.fna.before.ctrl$communication

LR.fna.anticd20.before <- netVisual_bubble(cellchat.merged.exhausted, sources.use = t.cells, targets.use = b.cells,  comparison = c(3, 2), angle.x = 45, remove.isolate = F, return.data = T)
LR.fna.anticd20.before.data <- LR.fna.anticd20.before$communication

LR.fna.plots <- lapply(list(control = LR.fna.before.ctrl.data,
                            antiCD20 = LR.fna.anticd20.before.data), 
                       function(x){
                              
                              plot = plot_split_LR_bubble(
                                LR.data = x, 
                                seurat_obj = seurat.fna.anticd20, 
                                celltype_col = "exhausted_included", 
                                condition_col = "DMT"
                              )
                              
                              return(plot)
                            })


LR.fna.plots$antiCD20
ggsave(filename = "figures/cellchat/exhausted_FNA_LR_antiCD20.pdf", width = 780, height = 880, units = "mm")
LR.fna.plots$control
ggsave(filename = "figures/cellchat/exhausted_FNA_LR_beforetreatment.pdf", width = 780, height = 880, units = "mm")
