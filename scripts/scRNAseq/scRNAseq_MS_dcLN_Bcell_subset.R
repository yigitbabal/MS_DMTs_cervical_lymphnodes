setwd("P:/grp_laakso/Project_repository/Yigit/MS_FNA_followup")

library(tidyverse)
library(Seurat)
library(data.table)
library(SingleR)
library(SingleCellExperiment)
library(plyr)

seurat.fna <- read_rds("seurat_objects/FNA/seurat_fna_integrated_singleR_annotated_curated.rds")
seurat.fna$patient <- str_split(seurat.fna$sample, pattern = "_", simplify = T)[,1]

seurat.fna$cellid <- colnames(seurat.fna)
seurat.fna$cellid <- gsub("-", "_", seurat.fna$cellid)

unique(seurat.fna$SingleR_labels_old)
unique(seurat.fna$Celltypes_curated_detailed)

b.cells <- c("Naive B cells", "Memory B cells", "Pre B cells")

seurat.fna.bcell <- subset(seurat.fna, subset = Celltypes_curated_detailed %in% c(b.cells))

seurat.fna.bcell
DefaultAssay(seurat.fna.bcell) <- "RNA"
seurat.fna.bcell[["RNA"]] <- split(seurat.fna.bcell[["RNA"]], f = seurat.fna.bcell$orig.ident)
seurat.fna.bcell
seurat.fna.bcell <- NormalizeData(seurat.fna.bcell, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.fna.bcell <- FindVariableFeatures(seurat.fna.bcell, nfeatures = 2000)
seurat.fna.bcell <- ScaleData(seurat.fna.bcell, features = VariableFeatures(seurat.fna.bcell))
seurat.fna.bcell <- RunPCA(seurat.fna.bcell, features = VariableFeatures(seurat.fna.bcell))
seurat.fna.bcell <- IntegrateLayers(
  object = seurat.fna.bcell, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = TRUE, k.weight = 62
)
seurat.fna.bcell <- FindNeighbors(seurat.fna.bcell, reduction = "integrated.rpca",  dims = 1:20)
seurat.fna.bcell <- FindClusters(seurat.fna.bcell, resolution = 0.7, cluster.name = "bcell_cluster_0.7")
seurat.fna.bcell <- RunUMAP(seurat.fna.bcell, reduction = "integrated.rpca", dims = 1:20)

DimPlot(seurat.fna.bcell, group.by = "bcell_cluster_0.7", label = T)

## T cell contamination
FeaturePlot(seurat.fna.bcell, features = c("CD3D", "CD3E", "CD4", "CD8A", "CD8B"))
seurat.fna.bcell <- subset(seurat.fna.bcell, subset = bcell_cluster_0.7 != 5 & 
                             bcell_cluster_0.7 != 7 & bcell_cluster_0.7 != 10 & 
                             bcell_cluster_0.7 != 13)

seurat.fna.bcell <- RunUMAP(seurat.fna.bcell, reduction = "integrated.rpca", dims = 1:20)
seurat.fna.bcell <- FindNeighbors(seurat.fna.bcell, reduction = "integrated.rpca",  dims = 1:20)
seurat.fna.bcell <- FindClusters(seurat.fna.bcell, resolution = 0.7, cluster.name = "bcell_cluster_0.7")
DimPlot(seurat.fna.bcell, group.by = "bcell_cluster_0.7", label = T)
seurat.fna.bcell <- FindSubCluster(seurat.fna.bcell, cluster = 6, graph.name = "RNA_snn", subcluster.name = "sub_6", resolution = 0.2)
DimPlot(seurat.fna.bcell, group.by = "sub_6", label = T)

FeaturePlot(seurat.fna.bcell, features = c("IGHD", "CD27", "IGHM"))

bcell.markers.res.0.7 <- FindAllMarkers(seurat.fna.bcell,
                                        assay = "RNA",
                                        slot = "data",
                                        test.use = "wilcox",
                                        only.pos = T,
                                        min.pct = 0.1,
                                        logfc.threshold = 0.25,
                                        group.by = "sub_6")

bcell.markers.res.0.7.sig <- bcell.markers.res.0.7[bcell.markers.res.0.7$p_val_adj < 0.05,]
bcell.markers.res.0.7.sig$cluster <- as.character(bcell.markers.res.0.7.sig$cluster)
fwrite(bcell.markers.res.0.7.sig, file = "results/FNA/annocation_curated/bcell_markers_res_0.7.csv")

FeaturePlot(seurat.fna.bcell, features = c("IGHD", "CD27", "IGHM"))
FeaturePlot(seurat.fna.bcell, features = c("TBX21", "ITGAX"))
FeaturePlot(seurat.fna.bcell, features = c("IFIT1", "IFI44L"))

FeaturePlot(seurat.fna.bcell, features = c("IGHD", "CD27", "IGHM","CR2",
                                           "TBX21", "ITGAX",
                                           "IFIT1", "IFI44L", "MME"))


subclusters <- c("Naive B" = "1",
                 "Naive B" = "3",
                 "Naive B" = "4",
                 "Naive B" = "5",
                 "Naive B" = "7",
                 "IFN-stimulated B" = "8",
                 "USM" = "6_0",
                 "USM" = "0",
                 "USM" = "9",
                 "SM" = "2",
                 "SM" = "6_1",
                 "Tbet+CD11c+" = "6_2",
                 "DN" = "10",
                 "Pre-B" = "12",
                 "Naive B" = "11")

Idents(seurat.fna.bcell) <- as.character(seurat.fna.bcell$sub_6)
seurat.fna.bcell$Bsubclusters <- plyr::mapvalues(seurat.fna.bcell$sub_6, from = subclusters, to = names(subclusters))

DimPlot(seurat.fna.bcell, group.by = "Bsubclusters", label = T)

VlnPlot(seurat.fna.bcell, features = c("MS4A1", "CD19" , "CD79A","IGHD", 
                                       "CD27", "IGHM", "IGHA1", "IGHA2", 
                                       "CR2", "CXCR5", "MME", "IGHG1", 
                                       "IGHG2", "IGHG3", "IGHG4", "IGHE",
                                       "IFIT1", "IFI44L", "TBX21", "ITGAX"), 
        group.by = "Bsubclusters", pt.size = 0)


write_rds(seurat.fna.bcell, file = "seurat_objects/FNA/seurat_fna_bcell_subset.rds")
