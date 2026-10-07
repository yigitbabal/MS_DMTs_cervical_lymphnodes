setwd("P:/grp_laakso/Project_repository/Yigit/MS_FNA_followup")

library(tidyverse)
library(Seurat)
library(data.table)
library(SingleR)
library(SingleCellExperiment)
library(plyr)

seurat.pbmc <- read_rds("seurat_objects/PBMC/seurat_pbmc_integrated_singleR_annotated_curated.rds")

seurat.pbmc$patient <- str_split(seurat.pbmc$sample, pattern = "_", simplify = T)[,1]

seurat.pbmc$cellid <- colnames(seurat.pbmc)
seurat.pbmc$cellid <- gsub("-", "_", seurat.pbmc$cellid)


unique(seurat.pbmc$Celltypes_curated_detailed)

b.cells <- c("Naive B cells", "Memory B cells")

seurat.pbmc.bcell <- subset(seurat.pbmc, subset = Celltypes_curated_detailed %in% c(b.cells))

seurat.pbmc.bcell
DefaultAssay(seurat.pbmc.bcell) <- "RNA"
seurat.pbmc.bcell[["RNA"]] <- split(seurat.pbmc.bcell[["RNA"]], f = seurat.pbmc.bcell$sample)
seurat.pbmc.bcell
seurat.pbmc.bcell <- NormalizeData(seurat.pbmc.bcell, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.pbmc.bcell <- FindVariableFeatures(seurat.pbmc.bcell, nfeatures = 2000)
seurat.pbmc.bcell <- ScaleData(seurat.pbmc.bcell, features = VariableFeatures(seurat.pbmc.bcell))
seurat.pbmc.bcell <- RunPCA(seurat.pbmc.bcell, features = VariableFeatures(seurat.pbmc.bcell))
seurat.pbmc.bcell <- IntegrateLayers(
  object = seurat.pbmc.bcell, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = TRUE,  dims = 1:7 
)
seurat.pbmc.bcell <- FindNeighbors(seurat.pbmc.bcell, reduction = "integrated.rpca",  dims = 1:20)
seurat.pbmc.bcell <- FindClusters(seurat.pbmc.bcell, resolution = 0.4, cluster.name = "bcell_cluster_0.4")
seurat.pbmc.bcell <- RunUMAP(seurat.pbmc.bcell, reduction = "integrated.rpca", dims = 1:20)

DimPlot(seurat.pbmc.bcell, group.by = "bcell_cluster_0.4", label = T)
FeaturePlot(seurat.pbmc.bcell, features = c("CD3D", "CD3E", "CD4", "CD8A", "CD8B"))
VlnPlot(seurat.pbmc.bcell, features = c("CD3D", "CD3E", "CD4", "CD8A", "CD8B"), pt.size = 0)

seurat.pbmc.bcell <- subset(seurat.pbmc.bcell, subset = bcell_cluster_0.4 != 4 & 
                              bcell_cluster_0.4 != 5)

seurat.pbmc.bcell <- FindNeighbors(seurat.pbmc.bcell, reduction = "integrated.rpca",  dims = 1:20)
seurat.pbmc.bcell <- FindClusters(seurat.pbmc.bcell, resolution = 0.7, cluster.name = "bcell_cluster_0.7")
seurat.pbmc.bcell <- FindClusters(seurat.pbmc.bcell, resolution = 0.8, cluster.name = "bcell_cluster_0.8")
seurat.pbmc.bcell <- FindClusters(seurat.pbmc.bcell, resolution = 1, cluster.name = "bcell_cluster_1")
seurat.pbmc.bcell <- FindClusters(seurat.pbmc.bcell, resolution = 1.2, cluster.name = "bcell_cluster_1.2")

seurat.pbmc.bcell <- RunUMAP(seurat.pbmc.bcell, reduction = "integrated.rpca", dims = 1:20)

DimPlot(seurat.pbmc.bcell, group.by = "bcell_cluster_0.7", label = T)
DimPlot(seurat.pbmc.bcell, group.by = "bcell_cluster_0.8", label = T)
DimPlot(seurat.pbmc.bcell, group.by = "bcell_cluster_1", label = T)

seurat.pbmc.bcell$seurat_clusters <- seurat.pbmc.bcell$bcell_cluster_0.7
seurat.pbmc.bcell <- FindSubCluster(seurat.pbmc.bcell, cluster = 4, graph.name = "RNA_snn", subcluster.name = "sub_4", resolution = 0.3)

DimPlot(seurat.pbmc.bcell, group.by = "sub_4", label = T)
FeaturePlot(seurat.pbmc.bcell, features = c("TBX21", "ITGAX"))


DimPlot(seurat.pbmc.bcell, group.by = "Celltypes_curated_detailed", label = T)
FeaturePlot(seurat.pbmc.bcell, features = c("IGHD", "CD27", "IGHM"))
VlnPlot(seurat.pbmc.bcell, features = c("IGHD", "CD27", "IGHM"), pt.size = 0)

bcell.markers.res.0.7 <- FindAllMarkers(seurat.pbmc.bcell,
                                        assay = "RNA",
                                        slot = "data",
                                        test.use = "wilcox",
                                        only.pos = T,
                                        min.pct = 0.1,
                                        logfc.threshold = 0.25,
                                        group.by = "sub_4")


bcell.markers.res.0.7.sig <- bcell.markers.res.0.7 %>% filter(p_val_adj < 0.05)
fwrite(bcell.markers.res.0.7.sig, file = "results/PBMC/bcell_markers_res_0.7_sig.csv")

subclusters <- c("Naive B" = "0",
                 "Naive B" = "1",
                 "Naive B" = "2",
                 "Naive B" = "6",
                 "Naive B" = "8",
                 "Naive B" = "7",
                 "USM" = "3",
                 "USM" = "4_0",
                 "Tbet+CD11c+" = "4_1",
                 "SM" = "5")

seurat.pbmc.bcell$sub_4 <- as.character(seurat.pbmc.bcell$sub_4)
seurat.pbmc.bcell$Bsubclusters <- plyr::mapvalues(seurat.pbmc.bcell$sub_4, from = subclusters, to = names(subclusters))

DimPlot(seurat.pbmc.bcell, group.by = "Bsubclusters", label = T)

VlnPlot(seurat.pbmc.bcell, features = c("MS4A1", "CD19" , "CD79A","IGHD", 
                                       "CD27", "IGHM", "IGHA1", "IGHA2", 
                                       "CR2", "CXCR5", "IGHG1", 
                                       "IGHG2", "IGHG3", "IGHG4", 
                                       "IGHE", "TBX21", "ITGAX"), 
        group.by = "Bsubclusters", pt.size = 0)


write_rds(seurat.pbmc.bcell, file = "seurat_objects/PBMC/seurat_pbmc_bcell_subset.rds")
