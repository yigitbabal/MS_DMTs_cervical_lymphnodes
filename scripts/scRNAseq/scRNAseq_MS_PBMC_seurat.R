setwd("P:/grp_laakso/Project_repository/Yigit/MS_FNA_followup")

library(tidyverse)
library(Seurat)
library(org.Hs.eg.db)
library(data.table)

source("scripts/00_functions.R")

samples.names <- list.dirs(path = "rawdata/GEX/PBMC/", full.names = F)
samples.names <- samples.names[-1]

seurat.list.PBMC <- counts_to_seurat(dir.path = "rawdata/GEX/PBMC/",
                                    mt.pattern = "^MT-",
                                    remove.ncRNA = T,
                                    remove.snRNA = T,
                                    mincells = 3,
                                    minfeatures = 200)

write_rds(seurat.list.PBMC, "seurat_objects/seurat_list_PMBC.rds")
#seurat.list.PBMC <- read_rds("seurat_objects/seurat_list_PMBC.rds")

preQC <- lapply(seurat.list.PBMC, function(x){
  p <- preQCSeurat(object = x, 
                   sample_col = "sample", 
                   mito_pattern = "^MT-", 
                   k_mad = 2.5, 
                   show_guides = T, 
                   show_ribo = F, 
                   saveplot = T,
                   plot_path = "figures/preQC_PBMC")
  return(p)
})

seurat.list.PBMC <- lapply(1:length(samples.names), function(x){
  print(paste0("Filtering ", samples.names[x], " ..."))
  s = seurat.list.PBMC[[x]]
  s = subset(s, subset = nFeature_RNA > preQC[[x]]$thresholds$nFeature_low & nFeature_RNA < preQC[[x]]$thresholds$nFeature_high & 
               nCount_RNA > preQC[[x]]$thresholds$nCount_low & nCount_RNA < preQC[[x]]$thresholds$nCount_high & 
               percent_mito < preQC[[x]]$thresholds$percent_mito_high)
  return(s)
})

names(seurat.list.PBMC) <- samples.names

postQC <- lapply(seurat.list.PBMC, function(x){
  p <- preQCSeurat(object = x, 
                   sample_col = "sample", 
                   show_guides = F, 
                   show_ribo = F, 
                   saveplot = T,
                   plot_path = "figures/postQC_PBMC")
  return(p)
})

### important: MS022 has high mitochondrial perchange and bimodel distribution on nFeature_RNA
### Therefore additional manual filtering is applied
seurat.list.PBMC$`MS022-PBMC` <- subset(seurat.list.PBMC$`MS022-PBMC`, subset = percent_mito < 7 & nFeature_RNA > 1000)

postQC <- lapply(seurat.list.PBMC, function(x){
  p <- preQCSeurat(object = x, 
                   sample_col = "sample", 
                   show_guides = F, 
                   show_ribo = F, 
                   saveplot = T,
                   plot_path = "figures/postQC_PBMC_after_manual")
  return(p)
})




seurat.pbmc <- process_seurat(obj_list = seurat.list.PBMC,
                             project_name = "FNA",
                             integrate = TRUE, 
                             use_jackstraw = TRUE,
                             res_range = c(0.2, 0.5, 0.8, 1.2),
                             js_threshold = 0.05)


seurat.pbmc

write_rds(seurat.pbmc, "seurat_objects/seurat_pbmc_integrated.rds")

# read sample meta data
sample.meta <- fread("sample_meta_data.csv")

seurat.pbmc$group <- sample.meta$Group[match(seurat.pbmc$sample, sample.meta$PBMC_ID)]
seurat.pbmc$DMT <- sample.meta$DMT[match(seurat.pbmc$sample, sample.meta$PBMC_ID)]
seurat.pbmc$DMT_detail <- sample.meta$DMT_detail[match(seurat.pbmc$sample, sample.meta$PBMC_ID)]

DimPlot(seurat.pbmc, reduction = "umap", group.by = "DMT")
ggsave("figures/PBMC/seurat_pbmc_integrated_umap_by_DMT.png", width = 150, height = 150, units = "mm")

DimPlot(seurat.pbmc, reduction = "umap", group.by = "RNA_snn_res.0.2", label = T)
ggsave("figures/PBMC/seurat_pbmc_integrated_umap_by_cluster_res_0.2.png", width = 150, height = 150, units = "mm")


DimPlot(seurat.pbmc, reduction = "umap", group.by = "RNA_snn_res.0.5", label = T)
ggsave("figures/PBMC/seurat_pbmc_integrated_umap_by_cluster_res_0.5.png", width = 150, height = 150, units = "mm")

# looks good
DimPlot(seurat.pbmc, reduction = "umap", group.by = "RNA_snn_res.0.8", label = T)
ggsave("figures/PBMC/seurat_pbmc_integrated_umap_by_cluster_res_0.8.png", width = 150, height = 150, units = "mm")

DimPlot(seurat.pbmc, reduction = "umap", group.by = "RNA_snn_res.1.2", label = T)
ggsave("figures/PBMC/seurat_pbmc_integrated_umap_by_cluster_res_1.2.png", width = 150, height = 150, units = "mm")


annotationResuls = annotateCells(seurat.pbmc, cluster_level = F)

seurat.pbmc[["hpcaMain"]] <- annotationResuls[[1]]
seurat.pbmc[["hpcaFine"]] <- annotationResuls[[2]]
seurat.pbmc[["immMain"]] <- annotationResuls[[3]]
seurat.pbmc[["immFine"]] <- annotationResuls[[4]]

write_rds(seurat.pbmc, "seurat_objects/seurat_pbmc_integrated_singleR_annotated.rds")

DimPlot(seurat.pbmc, reduction = "umap", group.by = "hpcaMain", label = T)
ggsave("figures/PBMC/seurat_pbmc_integrated_umap_by_hpcaMain.png", width = 150, height = 150, units = "mm")

DimPlot(seurat.pbmc, reduction = "umap", group.by = "hpcaFine", label = T) + NoLegend()
ggsave("figures/PBMC/seurat_pbmc_integrated_umap_by_hpcaFine.png", width = 150, height = 150, units = "mm")

DimPlot(seurat.pbmc, reduction = "umap", group.by = "immMain", label = T)
ggsave("figures/PBMC/seurat_pbmc_integrated_umap_by_immMain.png", width = 150, height = 150, units = "mm")

DimPlot(seurat.pbmc, reduction = "umap", group.by = "immFine", label = T) + NoLegend()
ggsave("figures/PBMC/seurat_pbmc_integrated_umap_by_immFine.png", width = 150, height = 150, units = "mm")


# cluster level SingleR annotation
annotationResuls.cluster = annotateCells(seurat.pbmc, cluster_level = T, cluster_col = "RNA_snn_res.0.8")

seurat.pbmc[["hpcaMain_cluster"]] <- annotationResuls.cluster[[1]]
seurat.pbmc[["hpcaFine_cluster"]] <- annotationResuls.cluster[[2]]
seurat.pbmc[["immMain_cluster"]] <- annotationResuls.cluster[[3]]
seurat.pbmc[["immFine_cluster"]] <- annotationResuls.cluster[[4]]

DimPlot(seurat.pbmc, reduction = "umap", group.by = "hpcaMain_cluster", label = T)
ggsave("figures/PBMC/seurat_pbmc_integrated_umap_by_hpcaMain_cluster_level.png", width = 150, height = 150, units = "mm")

DimPlot(seurat.pbmc, reduction = "umap", group.by = "hpcaFine_cluster", label = T)
ggsave("figures/PBMC/seurat_pbmc_integrated_umap_by_hpcaFine_cluster_level.png", width = 300, height = 300, units = "mm")

DimPlot(seurat.pbmc, reduction = "umap", group.by = "immMain_cluster", label = T)
ggsave("figures/PBMC/seurat_pbmc_integrated_umap_by_immMain_cluster_level.png", width = 150, height = 150, units = "mm")

DimPlot(seurat.pbmc, reduction = "umap", group.by = "immFine_cluster", label = T)
ggsave("figures/PBMC/seurat_pbmc_integrated_umap_by_immFine_cluster_level.png", width = 150, height = 150, units = "mm")

write_rds(seurat.pbmc, "seurat_objects/seurat_pbmc_integrated_singleR_annotated.rds")
seurat.pbmc <- read_rds("seurat_objects/PBMC/seurat_pbmc_integrated_singleR_annotated.rds")


## marker genes for cluster RNA_snn_res.0.8
seurat.pbmc <- JoinLayers(seurat.pbmc)

markers.res.0.8 <- FindAllMarkers(seurat.pbmc,
                                  assay = "RNA",
                                  slot = "data",
                                  test.use = "wilcox",
                                  only.pos = T,
                                  min.pct = 0.25,
                                  logfc.threshold = 0.25,
                                  group.by = "RNA_snn_res.0.8")

fwrite(markers.res.0.8, "results/seurat_pbmc_markers_res_0.8.csv")

markers.res.0.8.sig <- markers.res.0.8 %>% filter(p_val_adj < 0.05)
fwrite(markers.res.0.8.sig, "results/seurat_pbmc_markers_res_0.8_sig.csv")


seurat.pbmc@meta.data %>%
  dplyr::select(hpcaFine_cluster,hpcaMain_cluster, immMain_cluster, immFine_cluster, RNA_snn_res.0.8) %>%
  group_by(RNA_snn_res.0.8) %>%
  summarise(hpcaFine_cluster = unique(hpcaFine_cluster),
            hpcaMain_cluster = unique(hpcaMain_cluster),
            immMain_cluster = unique(immMain_cluster),
            immFine_cluster = unique(immFine_cluster)) -> cluster_annotation_summary

fwrite(cluster_annotation_summary, "results/seurat_pbmc_cluster_annotation_summary.csv")

curated.annotation <- xlsx::read.xlsx("results/FNA/annocation_curated/Cluster_annotation_curated.xlsx", sheetIndex = 2)



seurat.pbmc$Celltypes_curated_main <- plyr::mapvalues(seurat.pbmc$RNA_snn_res.0.8, from = curated.annotation$Cluster, to = curated.annotation$Main)
unique(seurat.pbmc$Celltypes_curated_main)
seurat.pbmc$Celltypes_curated_detailed <- plyr::mapvalues(seurat.pbmc$RNA_snn_res.0.8, from = curated.annotation$Cluster, to = curated.annotation$Detailed)
unique(seurat.pbmc$Celltypes_curated_detailed)

DimPlot(seurat.pbmc, reduction = "umap", group.by = "Celltypes_curated_main", label = T, repel = T, label.box = T)
DimPlot(seurat.pbmc, reduction = "umap", group.by = "Celltypes_curated_detailed", label = T, repel = T, label.box = T)

write_rds(seurat.pbmc, "seurat_objects/PBMC/seurat_pbmc_integrated_singleR_annotated_curated.rds")

##############################










marker.summary <- fread("results/PBMC_Cluster Annotation Summary Table.csv")

# make vector from marker genes
marker.summary %>%
  summarise(markers = paste(Markers, collapse = ", ")) %>%
  as.character() %>%
  strsplit(split = ", ") %>%
  unlist() %>%
  unique() -> marker.genes


library(scCustomize)

Idents(seurat.pbmc) <- seurat.pbmc$RNA_snn_res.0.8

p.marker <- Clustered_DotPlot(seurat_object = seurat.pbmc, flip = T, features = marker.genes, k = 13, plot_km_elbow = F, colors_use_exp = Seurat::BlueAndRed())

pdf(file = "figures/pbmc_marker_genes_dotplot.pdf", width = 20, height = 10)
p.marker <- Clustered_DotPlot(seurat_object = seurat.pbmc, flip = T, features = marker.genes, k = 13, plot_km_elbow = F, colors_use_exp = Seurat::BlueAndRed())
dev.off()

# assign cluster names based on the annotation summary table
seurat.pbmc$lineage <- marker.summary$Compartment[match(seurat.pbmc$RNA_snn_res.0.8, marker.summary$`Cluster ID`)]
seurat.pbmc$Celltypes <- marker.summary$Annotation[match(seurat.pbmc$RNA_snn_res.0.8, marker.summary$`Cluster ID`)]


DimPlot(seurat.pbmc, reduction = "umap", group.by = "lineage", label = T)
DimPlot(seurat.pbmc, reduction = "umap", group.by = "Celltypes", label = T, raster = F, repel = T, label.box = T)
ggsave("figures/PBMC/seurat_pbmc_integrated_umap_by_celltype.png", width = 400, height = 350, units = "mm")

VlnPlot(seurat.pbmc, features = c("nFeature_RNA","nCount_RNA", "percent_mito"), ncol = 1, group.by = "sample", pt.size = 0, raster = F)
ggsave("figures/PBMC/PBMC_overall_QC_after_merge.pdf", width = 300, height = 300, units = "mm")


write_rds(seurat.pbmc, "seurat_objects/seurat_pbmc_integrated_singleR_annotated.rds")

