setwd("P:/grp_laakso/Project_repository/Yigit/MS_FNA_followup")

library(tidyverse)
library(Seurat)
library(org.Hs.eg.db)
library(data.table)

source("scripts/00_functions.R")


samples.names <- list.dirs(path = "rawdata/GEX/FNA/", full.names = F)
samples.names <- samples.names[-1]

seurat.list.FNA <- counts_to_seurat(dir.path = "rawdata/GEX/FNA/",
                                    mt.pattern = "^MT-",
                                    remove.ncRNA = T,
                                    remove.snRNA = T,
                                    mincells = 3,
                                    minfeatures = 200)

write_rds(seurat.list.FNA, "seurat_objects/FNA/seurat_list_FNA.rds")


preQC <- lapply(seurat.list.FNA, function(x){
  p <- preQCSeurat(object = x, 
                   sample_col = "sample", 
                   mito_pattern = "^MT-", 
                   k_mad = 2.5, 
                   show_guides = T, 
                   show_ribo = F, 
                   saveplot = T,
                   plot_path = "figures/preQC_FNA")
  return(p)
})

seurat.list.FNA.filtered <- lapply(1:length(samples.names), function(x){
  print(paste0("Filtering ", samples.names[x], " ..."))
  s = seurat.list.FNA[[x]]
  s = subset(s, subset = nFeature_RNA > preQC[[x]]$thresholds$nFeature_low & nFeature_RNA < preQC[[x]]$thresholds$nFeature_high & 
               nCount_RNA > preQC[[x]]$thresholds$nCount_low & nCount_RNA < preQC[[x]]$thresholds$nCount_high & 
               percent_mito < preQC[[x]]$thresholds$percent_mito_high)
  return(s)
})

names(seurat.list.FNA.filtered) <- samples.names

postQC <- lapply(seurat.list.FNA.filtered, function(x){
  p <- preQCSeurat(object = x, 
                   sample_col = "sample", 
                   show_guides = F, 
                   show_ribo = F, 
                   saveplot = T,
                   plot_path = "figures/postQC_FNA")
  return(p)
})



### Manual adjustment for filtering
seurat.fna <- merge(x = seurat.list.FNA.filtered[[1]], y = seurat.list.FNA.filtered[-1])
seurat.fna # 232960 
VlnPlot(seurat.fna, features = c("nFeature_RNA","nCount_RNA", "percent_mito"), ncol = 1, group.by = "sample", pt.size = 0, raster = F)
seurat.fna <- subset(seurat.fna, subset = nCount_RNA < 20000 & percent_mito < 10)
seurat.fna # 232272 
VlnPlot(seurat.fna, features = c("nFeature_RNA","nCount_RNA", "percent_mito"), ncol = 1, group.by = "sample", pt.size = 0, raster = F)
ggsave("figures/FNA/FNA_overall_QC_after_merge.pdf", width = 300, height = 300, units = "mm")



seurat.fna <- process_seurat(obj = seurat.fna,
                             project_name = "FNA",
                             integrate = TRUE, 
                             use_jackstraw = TRUE,
                             res_range = c(0.2, 0.4, 0.5, 0.7, 0.8, 1.2),
                             js_threshold = 0.05)


seurat.fna

write_rds(seurat.fna, "seurat_objects/FNA/seurat_fna_integrated.rds")

# read sample meta data
sample.meta <- fread("sample_meta_data.csv")

seurat.fna$group <- sample.meta$Group[match(seurat.fna$sample, sample.meta$ID)]
seurat.fna$DMT <- sample.meta$DMT[match(seurat.fna$sample, sample.meta$ID)]
seurat.fna$DMT_detail <- sample.meta$DMT_detail[match(seurat.fna$sample, sample.meta$ID)]
seurat.fna$Protocol <- sample.meta$Protocol[match(seurat.fna$sample, sample.meta$ID)]

DimPlot(seurat.fna, reduction = "umap", group.by = "Protocol")
ggsave("figures/FNA/seurat_fna_integrated_umap_by_protocol.png", width = 150, height = 150, units = "mm")

DimPlot(seurat.fna, reduction = "umap", group.by = "DMT")
ggsave("figures/FNA/seurat_fna_integrated_umap_by_DMT.png", width = 150, height = 150, units = "mm")

DimPlot(seurat.fna, reduction = "umap", group.by = "RNA_snn_res.0.2", label = T)
ggsave("figures/FNA/seurat_fna_integrated_umap_by_cluster_res_0.2.png", width = 150, height = 150, units = "mm")

DimPlot(seurat.fna, reduction = "umap", group.by = "RNA_snn_res.0.4", label = T)
ggsave("figures/FNA/seurat_fna_integrated_umap_by_cluster_res_0.4.png", width = 150, height = 150, units = "mm")

DimPlot(seurat.fna, reduction = "umap", group.by = "RNA_snn_res.0.5", label = T)
ggsave("figures/FNA/seurat_fna_integrated_umap_by_cluster_res_0.5.png", width = 150, height = 150, units = "mm")

DimPlot(seurat.fna, reduction = "umap", group.by = "RNA_snn_res.0.7", label = T)
ggsave("figures/FNA/seurat_fna_integrated_umap_by_cluster_res_0.7.png", width = 150, height = 150, units = "mm")

DimPlot(seurat.fna, reduction = "umap", group.by = "RNA_snn_res.0.8", label = T)
ggsave("figures/FNA/seurat_fna_integrated_umap_by_cluster_res_0.8.png", width = 150, height = 150, units = "mm")

# look good
DimPlot(seurat.fna, reduction = "umap", group.by = "RNA_snn_res.1.2", label = T)
ggsave("figures/FNA/seurat_fna_integrated_umap_by_cluster_res_1.2.png", width = 150, height = 150, units = "mm")



annotationResuls = annotateCells(seurat.fna, cluster_level = F)
seurat.fna[["hpcaMain"]] <- annotationResuls[[1]]
seurat.fna[["hpcaFine"]] <- annotationResuls[[2]]
seurat.fna[["immMain"]] <- annotationResuls[[3]]
seurat.fna[["immFine"]] <- annotationResuls[[4]]

DimPlot(seurat.fna, reduction = "umap", group.by = "hpcaMain", label = T)
ggsave("figures/FNA/seurat_fna_integrated_umap_by_hpcaMain.png", width = 150, height = 150, units = "mm")

DimPlot(seurat.fna, reduction = "umap", group.by = "hpcaFine", label = T) + NoLegend()
ggsave("figures/FNA/seurat_fna_integrated_umap_by_hpcaFine.png", width = 150, height = 150, units = "mm")

DimPlot(seurat.fna, reduction = "umap", group.by = "immMain", label = T)
ggsave("figures/FNA/seurat_fna_integrated_umap_by_hpcaFine.png", width = 150, height = 150, units = "mm")

DimPlot(seurat.fna, reduction = "umap", group.by = "immFine", label = T) + NoLegend()
ggsave("figures/FNA/seurat_fna_integrated_umap_by_immFine.png", width = 150, height = 150, units = "mm")

# cluster level SingleR annotation
annotationResuls.cluster = annotateCells(seurat.fna, cluster_level = T, cluster_col = "RNA_snn_res.1.2")

seurat.fna[["hpcaMain_cluster"]] <- annotationResuls.cluster[[1]]
seurat.fna[["hpcaFine_cluster"]] <- annotationResuls.cluster[[2]]
seurat.fna[["immMain_cluster"]] <- annotationResuls.cluster[[3]]
seurat.fna[["immFine_cluster"]] <- annotationResuls.cluster[[4]]

DimPlot(seurat.fna, reduction = "umap", group.by = "hpcaMain_cluster", label = T)
ggsave("figures/FNA/seurat_fna_integrated_umap_by_hpcaMain_cluster_level.png", width = 150, height = 150, units = "mm")

DimPlot(seurat.fna, reduction = "umap", group.by = "hpcaFine_cluster", label = T)
ggsave("figures/FNA/seurat_fna_integrated_umap_by_hpcaFine_cluster_level.png", width = 300, height = 300, units = "mm")

DimPlot(seurat.fna, reduction = "umap", group.by = "immMain_cluster", label = T)
ggsave("figures/FNA/seurat_fna_integrated_umap_by_immMain_cluster_level.png", width = 150, height = 150, units = "mm")

DimPlot(seurat.fna, reduction = "umap", group.by = "immFine_cluster", label = T)
ggsave("figures/FNA/seurat_fna_integrated_umap_by_immMain_cluster_level.png", width = 150, height = 150, units = "mm")

write_rds(seurat.fna, "seurat_objects/FNA/seurat_fna_integrated_singleR_annotated.rds")
seurat.fna <- read_rds("seurat_objects/FNA/seurat_fna_integrated_singleR_annotated.rds")

########################################
### annotation from old samples
load("annotation_validation/msFNAmerged_oldAndFU.RData")
sce_ref <- as.SingleCellExperiment(msFNAmerged.oldAndFU)
sce_query <- as.SingleCellExperiment(seurat.fna)
predictions <- SingleR(test = sce_query, 
                       ref = sce_ref, 
                       labels = sce_ref$refAnnMan)

seurat.fna$SingleR_labels_old <- predictions$pruned.labels
DimPlot(seurat.fna, reduction = "umap", group.by = "SingleR_labels_old", label = T)
DimPlot(seurat.fna, reduction = "umap", group.by = "Celltypes", label = T)


write_rds(seurat.fna, "seurat_objects/FNA/seurat_fna_integrated_singleR_annotated_with_old_labels.rds")
seurat.fna <- read_rds("seurat_objects/FNA/seurat_fna_integrated_singleR_annotated_with_old_labels.rds")
###################################



patchwork::wrap_plots(
  DimPlot(seurat.fna, reduction = "umap", group.by = "RNA_snn_res.0.2", label = T),
  DimPlot(seurat.fna, reduction = "umap", group.by = "RNA_snn_res.0.4", label = T),
  DimPlot(seurat.fna, reduction = "umap", group.by = "RNA_snn_res.0.5", label = T),
  DimPlot(seurat.fna, reduction = "umap", group.by = "RNA_snn_res.0.7", label = T),
  DimPlot(seurat.fna, reduction = "umap", group.by = "RNA_snn_res.0.8", label = T),
  DimPlot(seurat.fna, reduction = "umap", group.by = "RNA_snn_res.1.2", label = T)
)


#########################
#### Cluster stability analysis and tree
library(clustree)
library(ape)
clustree.res <- clustree(seurat.fna, prefix = "RNA_snn_res.")
fwrite(clustree.res$data, "results/FNA/seurat_fna_clustree_results.csv")
stability_summary <- clustree.res$data %>%
  # Group the data by your resolution column
  group_by(RNA_snn_res.) %>%
  
  # Calculate the metrics
  summarise(
    Total_Clusters = n_distinct(cluster),
    Mean_Stability = mean(sc3_stability, na.rm = TRUE),
    Median_Stability = median(sc3_stability, na.rm = TRUE),
    Unstable_Clusters = sum(sc3_stability < 0.1, na.rm = TRUE),
    Highly_Stable_Clusters = sum(sc3_stability > 0.5, na.rm = TRUE)
  ) %>%
  # Arrange by resolution for clean viewing
  arrange(RNA_snn_res.)

Idents(seurat.fna) <- seurat.fna$RNA_snn_res.0.5
seurat.fna <- BuildClusterTree(object = seurat.fna, reduction = "integrated.rpca", dims = 1:50, reorder = T)
PlotClusterTree(seurat.fna)


cluster_tree <- seurat.fna@tools[["BuildClusterTree"]]
hclust_tree <- as.hclust(cluster_tree)
k_groups <- 10
super_clusters <- cutree(hclust_tree, k = k_groups)
original_clusters <- as.character(seurat.fna$RNA_snn_res.0.5)
cell_super_clusters <- super_clusters[as.character(seurat.fna$RNA_snn_res.0.5)]
names(cell_super_clusters) <- colnames(seurat.fna)
seurat.fna$SuperCluster <- cell_super_clusters
cluster_mapping_df <- seurat.fna@meta.data %>%
  dplyr::select(Celltypes, RNA_snn_res.0.5, SuperCluster) %>%
  distinct() %>% # Removes duplicate rows to leave only the unique mapping
  arrange(SuperCluster, as.numeric(as.character(RNA_snn_res.0.5))) # Sorts neatly

DimPlot(seurat.fna, reduction = "umap", group.by = "SuperCluster", label = T)
#########################

## marker genes for cluster RNA_snn_res.1.2
seurat.fna <- JoinLayers(seurat.fna)

markers.res.1.2 <- FindAllMarkers(seurat.fna,
                                 assay = "RNA",
                                 slot = "data",
                                 test.use = "wilcox",
                                 only.pos = T,
                                 min.pct = 0.25,
                                 logfc.threshold = 0.25,
                                 group.by = "RNA_snn_res.1.2")

fwrite(markers.res.1.2, "results/FNA/seurat_fna_markers_res_1.2.csv")
markers.res.1.2 <- fread("results/FNA/seurat_fna_markers_res_1.2.csv")

markers.res.1.2.sig <- markers.res.1.2 %>% filter(p_val_adj < 0.05)
fwrite(markers.res.1.2.sig, "results/FNA/seurat_fna_markers_res_1.2_sig.csv")


seurat.fna@meta.data %>%
  dplyr::select(hpcaFine_cluster,hpcaMain_cluster, immMain_cluster, immFine_cluster, SingleR_labels_old, RNA_snn_res.1.2) %>%
  group_by(RNA_snn_res.1.2) %>%
  summarise(hpcaFine_cluster = unique(hpcaFine_cluster),
            hpcaMain_cluster = unique(hpcaMain_cluster),
            immMain_cluster = unique(immMain_cluster),
            immFine_cluster = unique(immFine_cluster),
            SingleR_labels_old = unique(SingleR_labels_old)) -> cluster_annotation_summary


fwrite(cluster_annotation_summary, "results/FNA/seurat_fna_cluster_annotation_summary.csv")


##### subclusters
Idents(seurat.fna) <- seurat.fna$RNA_snn_res.1.2
Idents(seurat.fna)
seurat.fna <- FindSubCluster(seurat.fna, cluster = 27, graph.name = "RNA_snn", subcluster.name = "27_sub", resolution = 0.07)
unique(seurat.fna$`27_sub`)

Idents(seurat.fna) <- seurat.fna$`27_sub`
DimPlot(seurat.fna, reduction = "umap", label = T, label.box = T)

# marker genes for 27_0 and 27_1
i2 <- unique(seurat.fna$`27_sub`)[!unique(seurat.fna$`27_sub`) %in% "27_0"]
marker_27_0 <- FindMarkers(seurat.fna, assay = "RNA", slot = "data", test.use = "wilcox", only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, ident.1 = "27_0") 
marker_27_0.sig <- marker_27_0 %>% filter(p_val_adj < 0.05)
marker_27_0.sig$gene <- rownames(marker_27_0.sig)
fwrite(marker_27_0.sig, "results/FNA/seurat_fna_markers_cluster_27_0.csv")

marker_27_1 <- FindMarkers(seurat.fna, assay = "RNA", slot = "data", test.use = "wilcox", only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, ident.1 = "27_1") 
marker_27_1.sig <- marker_27_1 %>% filter(p_val_adj < 0.05)
marker_27_1.sig$gene <- rownames(marker_27_1.sig)
fwrite(marker_27_1.sig, "results/FNA/seurat_fna_markers_cluster_27_1.csv")

#####

# annotation curated
curated.annotation <- xlsx::read.xlsx("results/FNA/annocation_curated/Cluster_annotation_curated.xlsx", sheetIndex = 1)

seurat.fna$Celltypes_curated_main <- plyr::mapvalues(seurat.fna$seurat_clusters, from = curated.annotation$Cluster, to = curated.annotation$Main)
unique(seurat.fna$Celltypes_curated_main)
seurat.fna$Celltypes_curated_detailed <- plyr::mapvalues(seurat.fna$seurat_clusters, from = curated.annotation$Cluster, to = curated.annotation$Detailed)
unique(seurat.fna$Celltypes_curated_detailed)

DimPlot(seurat.fna, reduction = "umap", group.by = "Celltypes_curated_main", label = T, repel = T, label.box = T)
DimPlot(seurat.fna, reduction = "umap", group.by = "Celltypes_curated_detailed", label = T, repel = T, label.box = T)

write_rds(seurat.fna, "seurat_objects/FNA/seurat_fna_integrated_singleR_annotated_curated.rds")


curated.annotation$Marker <- gsub(", ", ",", curated.annotation$Marker)
curated.annotation %>%
  summarise(markers = paste(Marker, collapse = ",")) %>%
  as.character() %>%
  strsplit(split = ",") %>%
  unlist() %>%
  unique() -> marker.genes

marker.genes

library(scCustomize)

Idents(seurat.fna) <- seurat.fna$RNA_snn_res.1.2

p.marker <- Clustered_DotPlot(seurat_object = seurat.fna, flip = T, features = marker.genes, k = 13, plot_km_elbow = F, colors_use_exp = Seurat::BlueAndRed())

pdf(file = "results/FNA/annocation_curated/marker_genes_curated_dotplot.pdf", width = 30, height = 10)
p.marker <- Clustered_DotPlot(seurat_object = seurat.fna, flip = T, features = marker.genes, plot_km_elbow = F, colors_use_exp = Seurat::BlueAndRed())
dev.off()


Idents(seurat.fna) <- seurat.fna$Celltypes_curated_detailed

pdf(file = "results/FNA/annocation_curated/marker_genes_curated_dotplot_celltype.pdf", width = 30, height = 10)
p.marker <- Clustered_DotPlot(seurat_object = seurat.fna, flip = T, features = marker.genes, plot_km_elbow = F, colors_use_exp = Seurat::BlueAndRed())
dev.off()



######


marker.summary <- fread("results/FNA/Cluster Annotation Summary Table.csv")

# make vector from marker genes
marker.summary %>%
  summarise(markers = paste(Markers, collapse = ", ")) %>%
  as.character() %>%
  strsplit(split = ", ") %>%
  unlist() %>%
  unique() -> marker.genes

marker.genes

library(scCustomize)

Idents(seurat.fna) <- seurat.fna$RNA_snn_res.1.2

p.marker <- Clustered_DotPlot(seurat_object = seurat.fna, flip = T, features = marker.genes, k = 13, plot_km_elbow = F, colors_use_exp = Seurat::BlueAndRed())

pdf(file = "figures/FNA/marker_genes_dotplot.pdf", width = 20, height = 10)
p.marker <- Clustered_DotPlot(seurat_object = seurat.fna, flip = T, features = marker.genes, k = 13, plot_km_elbow = F, colors_use_exp = Seurat::BlueAndRed())
dev.off()

# assign cluster names based on the annotation summary table
seurat.fna$lineage <- marker.summary$Lineage[match(seurat.fna$RNA_snn_res.1.2, marker.summary$`Cluster ID`)]
seurat.fna$Celltypes <- marker.summary$Annotation[match(seurat.fna$RNA_snn_res.1.2, marker.summary$`Cluster ID`)]


DimPlot(seurat.fna, reduction = "umap", group.by = "lineage", label = T)
DimPlot(seurat.fna, reduction = "umap", group.by = "Celltypes", label = T, raster = F, repel = T, label.box = T)
ggsave("figures/FNA/seurat_fna_integrated_umap_by_celltype.png", width = 400, height = 350, units = "mm")

write_rds(seurat.fna, "seurat_objects/FNA/seurat_fna_integrated_singleR_annotated.rds")


Idents(seurat.fna) <- seurat.fna$Celltypes

pdf(file = "figures/FNA/marker_genes_dotplot_celltype.pdf", width = 30, height = 10)
p.marker.celltype <- Clustered_DotPlot(seurat_object = seurat.fna, flip = T, features = marker.genes, k = 13, plot_km_elbow = F, colors_use_exp = Seurat::BlueAndRed())
dev.off()


cluster_states(seurat.fna, "sample", "Celltypes")




