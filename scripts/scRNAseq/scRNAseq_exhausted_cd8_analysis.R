setwd("P:/grp_laakso/Project_repository/Yigit/MS_FNA_followup")

library(tidyverse)
library(Seurat)
library(data.table)
library(RColorBrewer)

source("scripts/00_functions.R")

#seurat objects
seurat.fna <- read_rds("seurat_objects/FNA/seurat_fna_curated_with_subset.rds")
seurat.pbmc <- read_rds("seurat_objects/PBMC/seurat_pbmc_curated_with_subset.rds")

exhaustion.markers <- c("TIGIT", "HAVCR2", "PDCD1", "CD244", "LAG3", "TCF7", "CTLA4") 


seurat.fna.cd8 <- subset(seurat.fna, subset = Celltypes_curated_detailed_withsub %in% c("Naive CD8", "Early activated CD8","Memory CD8"))
seurat.fna.cd8 <- AddModuleScore(seurat.fna.cd8, assay = "RNA",
                                 features = list(exhaustion.markers),ctrl = 100,
                                 name = 'Exhaustion_score')
seurat.fna.cd8 <- RunUMAP(seurat.fna.cd8, reduction = "integrated.rpca", dims = 1:50)
seurat.fna.cd8 <- FindClusters(seurat.fna.cd8, resolution = 0.7, cluster.name = "CD8_res_0.7")
DimPlot(seurat.fna.cd8, group.by = "CD8_res_0.7", label = T)

DimPlot(seurat.fna.cd8, group.by = "Celltypes_curated_detailed_withsub")
FeaturePlot(seurat.fna.cd8, features = "Exhaustion_score1", cols=c("lightblue","blue","darkblue")) 

library(scCustomize)
library(viridis)
pal <- viridis(n = 10, option = "D")
FeaturePlot_scCustom(seurat_object = seurat.fna.cd8, features = "Exhaustion_score1", colors_use = c("lightblue","blue","darkblue"), order = T)

DotPlot(seurat.fna.cd8, features = exhaustion.markers, group.by = "CD8_res_0.7")


seurat.fna.cd8.exhaust <- subset(seurat.fna.cd8, subset = CD8_res_0.7 %in% c(5,6,13,15))


nExhaustionPerCD8Cluster <- seurat.fna.cd8.exhaust@meta.data %>% group_by(sample) %>% 
  summarise(n=ncol(seurat.fna.cd8),nEX = sum(Exhaustion_score1 > 0)) %>% mutate(pEX = (nEX/n) * 100) %>% arrange(desc(pEX))

sample.meta <- fread("sample_meta_data.csv")

nExhaustionPerCD8Cluster$DMT = sample.meta$DMT[match(nExhaustionPerCD8Cluster$sample, sample.meta$ID)]
nExhaustionPerCD8Cluster$DMT = factor(nExhaustionPerCD8Cluster$DMT, levels = c("Control", "BeforeTreatment","antiCD20", "DMF", "NTZ"))


nExhaustionPerCD8Cluster %>%
  ggplot(aes(x = DMT, y = pEX, fill = DMT)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 3) +
  labs(title = "", x = "", y = "Exhausted CD8 of all CD8 T cells (%)") +
  ggsignif::geom_signif(comparisons = list(c("BeforeTreatment", "Control"), c("BeforeTreatment", "antiCD20"), c("BeforeTreatment", "DMF")), 
                        map_signif_level = TRUE, test = "t.test", color = "black", step_increase = 0.1) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) + 
  theme(legend.position = "none")


nExhaustionPerCD8Cluster %>%
  filter(DMT %in% c("Control", "BeforeTreatment", "antiCD20")) %>%
  ggplot(aes(x = DMT, y = pEX, fill = DMT)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 3) +
  labs(title = "", x = "", y = "Exhausted CD8 of all CD8 T cells (%)") +
  ggsignif::geom_signif(comparisons = list(c("BeforeTreatment", "Control"), c("BeforeTreatment", "antiCD20")), 
                        map_signif_level = TRUE, test = "t.test", color = "black", step_increase = 0.1) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) + 
  theme(legend.position = "none")

write_rds(seurat.fna.cd8, file = "seurat_objects/seurat.fna.cd8.rds")
write_rds(seurat.fna.cd8.exhaust, file = "seurat_objects/seurat.fna.cd8.exhaust.rds")


seurat.fna.cd8 <- read_rds("seurat_objects/seurat.fna.cd8.rds")


cd8.sub.markers <- FindAllMarkers(seurat.fna.cd8, assay = "RNA",
                                  only.pos = TRUE, 
                                  min.pct = 0.25, 
                                  logfc.threshold = 0.25)

cd8.sub.markers.sig <- cd8.sub.markers[cd8.sub.markers$p_val_adj < 0.05,]
fwrite(cd8.sub.markers.sig, file = "results/cd8_exhaust/cd8_subcluster_markers.csv", row.names = T)

# Memory CD8 T cell subcluster marker gene vectors
# Clusters: 8, 9, 11, 14 (non-exhausted) | 5, 6, 13, 15 (exhausted)

c8.marker.genes <- c("SELL", "CCR7", "TCF7", "LEF1",
                     "CCL5", "GZMK", "GZMA", "NKG7",
                     "KLRG1", "EOMES", "TRDV1", "TRDV2", "TRDC", 
                     "CTLA4", "HAVCR2", "LAG3", "PDCD1", "TIGIT")

c8.subtypes <- c("Naive (c0)",
                 "Naive (c1)",
                 "Naive (c2)",
                 "Naive (c3)",
                 "Naive (c4)",
                 "TRM-EX (c5)",
                 "TRM-TEX (c6)",
                 "Naive (c7)",
                 "TEM (c8)",
                 "TEM (c9)",
                 "Naive (c10)",
                 "γδ T (c11)",
                 "Early activated (c12)",
                 "TRM-TEX (c13)",
                 "γδ T (c14)",
                 "TEX (c15)"
)

names(c8.subtypes) <- levels(seurat.fna.cd8)
seurat.fna.cd8 <- RenameIdents(seurat.fna.cd8, c8.subtypes)
seurat.fna.cd8$subtypes <- Idents(seurat.fna.cd8)

write_rds(seurat.fna.cd8, file = "seurat_objects/seurat.fna.cd8.rds")
seurat.fna.cd8 <- read_rds("seurat_objects/seurat.fna.cd8.rds")
seurat.fna.cd8.exhaust <- read_rds("seurat_objects/seurat.fna.cd8.exhaust.rds")

seurat.fna.cd8.exhaust$subtypes <- seurat.fna.cd8$subtypes[match(colnames(seurat.fna.cd8.exhaust), names(seurat.fna.cd8$subtypes))]
unique(seurat.fna.cd8.exhaust$subtypes)

nExhaustionPerCD8Cluster <- seurat.fna.cd8.exhaust@meta.data %>% group_by(sample) %>% 
  summarise(n=nrow(seurat.fna.cd8@meta.data[seurat.fna.cd8@meta.data$Celltypes_curated_detailed_withsub == "Memory CD8",]),nEX = sum(Exhaustion_score1 > 0)) %>% mutate(pEX = (nEX/n) * 100) %>% arrange(desc(pEX))

sample.meta <- fread("sample_meta_data.csv")

nExhaustionPerCD8Cluster$DMT = sample.meta$DMT[match(nExhaustionPerCD8Cluster$sample, sample.meta$ID)]
nExhaustionPerCD8Cluster$DMT = factor(nExhaustionPerCD8Cluster$DMT, levels = c("Control", "BeforeTreatment","antiCD20", "DMF", "NTZ"))


# total Memory CD8 cells per sample — single shared denominator across all subtypes
total_memcd8_per_sample <- seurat.fna.cd8@meta.data %>%
  #filter(Celltypes_curated_detailed_withsub == "Memory CD8") %>%
  group_by(sample) %>%
  summarise(n = n(), .groups = "drop")

# exhausted cells per sample per subtype (numerator)
nExhaustionPerCD8Subtype <- seurat.fna.cd8.exhaust@meta.data %>%
  group_by(sample, subtypes) %>%
  summarise(nEX = sum(Exhaustion_score1 > 0), .groups = "drop") %>%
  left_join(total_memcd8_per_sample, by = "sample") %>%
  mutate(pEX = (nEX / n) * 100) %>%
  arrange(desc(pEX))

sample.meta <- fread("sample_meta_data.csv")
nExhaustionPerCD8Subtype$DMT = sample.meta$DMT[match(nExhaustionPerCD8Subtype$sample, sample.meta$ID)]
nExhaustionPerCD8Subtype$DMT = factor(nExhaustionPerCD8Subtype$DMT, levels = c("Control", "BeforeTreatment","antiCD20", "DMF", "NTZ"))
nExhaustionPerCD8Subtype$DMT_detail = sample.meta$DMT_detail[match(nExhaustionPerCD8Subtype$sample, sample.meta$ID)]
nExhaustionPerCD8Subtype$patient <- gsub(pattern = "_FU", "", x = nExhaustionPerCD8Subtype$sample)
nExhaustionPerCD8Subtype$relapse <- "Non-relapsing"
nExhaustionPerCD8Subtype[nExhaustionPerCD8Subtype$patient == "MS022", "relapse"] <- "Relapse-during-DMT"


nExhaustionPerCD8Subtype %>%
  filter(DMT %in% c("Control", "BeforeTreatment", "antiCD20")) %>%
  ggplot(aes(x = DMT, y = pEX)) +
  geom_boxplot(aes(fill = DMT), outlier.shape = NA) +
  geom_line(aes(group = patient),
            color = "grey50", alpha = 0.6, linewidth = 0.4) +
  geom_point(aes(color = DMT_detail, shape = relapse),
             size = 3) +
  scale_color_manual(values = c(
    "Control" = "#8C8C8C",
    "BeforeTreatment" = "#4D4D4D",
    "antiCD20 (rituximab)" = "#0073C2",
    "antiCD20 (ofatumumab)" = "#EFC000"
  )) +
  labs(title = "", x = "", y = "Exhausted CD8 of all CD8 T cells (%)") +
  ggsignif::geom_signif(comparisons = list(c("BeforeTreatment", "Control"), c("BeforeTreatment", "antiCD20")), 
                        map_signif_level = TRUE, test = "t.test", color = "black", step_increase = 0.1) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) + 
  guides(fill = "none") +
  facet_wrap(~subtypes, scales = "free")






Idents(seurat.fna.cd8) <- seurat.fna.cd8$subtypes
cd8.subcelltype.markers <- FindAllMarkers(seurat.fna.cd8, assay = "RNA",
                                  only.pos = TRUE, 
                                  min.pct = 0.25, 
                                  logfc.threshold = 0.25)

cd8.subcelltype.markers <- cd8.subcelltype.markers[cd8.subcelltype.markers$p_val_adj < 0.05,]
cd8.subcelltype.markers.exhaust <- cd8.subcelltype.markers[cd8.subcelltype.markers$cluster %in% c("TRM-EX (c5)",
                                                                                                  "TRM-TEX (c6)",
                                                                                                  "TRM-TEX (c13)",
                                                                                                  "TEX (c15)"),]


seurat.fna.cd8.exhausted <- subset(seurat.fna.cd8, subset = subtypes %in% c("TRM-EX (c5)",
                                                                            "TRM-TEX (c6)",
                                                                            "TRM-TEX (c13)",
                                                                            "TEX (c15)"))

unique(seurat.fna.cd8.exhausted$subtypes)
unique(seurat.fna.cd8.exhausted$DMT)
unique(seurat.fna.cd8.exhausted$sample)


deg.exhaust.list.baseline <- RunConsensusPseudobulk(seurat_obj = seurat.fna.cd8.exhausted, 
                                        cell_type_col = "subtypes", 
                                        treatment_col = "DMT",
                                        sample_col = "sample", 
                                        baseline_group = "Control", 
                                        after_groups = c("BeforeTreatment"),
                                        p_val_cutoff = 0.05)

deg.exhaust.list.anticd20 <- RunConsensusPseudobulk(seurat_obj = seurat.fna.cd8.exhausted, 
                                                    cell_type_col = "subtypes", 
                                                    treatment_col = "DMT",
                                                    sample_col = "sample", 
                                                    baseline_group = "BeforeTreatment", 
                                                    after_groups = c("antiCD20"),
                                                    p_val_cutoff = 0.05)


write_csv(deg.exhaust.list.anticd20$significant_consensus, file = "results/cd8_exhaust/cd8_exhaust_pseudobulk_anticd20_vs_before_DEGs.csv")

library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(dplyr)

cds <- as.cell_data_set(seurat.fna.cd8)

# Transfer Seurat's UMAP embedding to avoid recomputing
reducedDims(cds)[["UMAP"]] <- Embeddings(seurat.fna.cd8, "umap")

# Transfer cluster info from Seurat
cds@clusters$UMAP$clusters <- seurat.fna.cd8$CD8_res_0.7

cds <- cluster_cells(cds, resolution = 1e-3)

# Or reuse Seurat clusters (recommended for consistency)
cds@clusters$UMAP$cluster_result$optim_res$membership <-
  setNames(seurat.fna.cd8$CD8_res_0.7, colnames(seurat.fna.cd8))

cds <- learn_graph(cds, use_partition = FALSE)
# use_partition = FALSE → fits one graph across all clusters
# use_partition = TRUE  → fits separate graphs per partition

plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster = TRUE,
           label_leaves = TRUE,
           label_branch_points = TRUE)

cds <- order_cells(cds)



# Visualize pseudotime
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 4)


write_rds(cds, file = "results/cd8_exhaust/cds_fna_cd8.rds")
cds <- read_rds("results/cd8_exhaust/cds_fna_cd8.rds")
# Test which genes change significantly along the trajectory
deg_results <- graph_test(cds, neighbor_graph = "principal_graph", cores = 10)

write_csv(deg_results, file = "results/cd8_exhaust/cds_fna_cd8_deg_results.csv")
deg_results <- read_csv("results/cd8_exhaust/cds_fna_cd8_deg_results.csv")
# Filter significant genes
sig_genes <- deg_results %>%
  filter(q_value < 0.05, morans_I > 0.1) %>%
  arrange(desc(morans_I))

sig_genes.exhaust <- sig_genes[rownames(sig_genes) %in% exhaustion.markers,]

rownames(sig_genes.exhaust)

rowData(cds)$gene_short_name <- rownames(cds)

plot.genes <- rownames(sig_genes.exhaust)

# 1. Extract pseudotime from CDS
pseudotime_df <- data.frame(
  cell        = colnames(cds),
  pseudotime  = pseudotime(cds),
  cluster     = cds@clusters$UMAP$clusters
)

# 2. Extract normalized expression from Seurat (log-normalized)
expr_mat <- GetAssayData(seurat.fna.cd8, assay = "RNA", slot = "data")[plot.genes, ]

# 3. Combine into long-format data frame
expr_df <- as.data.frame(t(as.matrix(expr_mat))) %>%
  tibble::rownames_to_column("cell") %>%
  pivot_longer(-cell, names_to = "gene", values_to = "expression") %>%
  left_join(pseudotime_df, by = "cell") %>%
  filter(is.finite(pseudotime))

ggplot(expr_df, aes(x = pseudotime, y = expression, color = pseudotime)) +
  scale_color_viridis_c(option = "plasma") +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
              color = "black", linewidth = 0.8) +
  facet_wrap(~ gene, scales = "free_y", ncol = 2) +
  theme_classic(base_size = 12)


# cytokines 

cytokines <- unique(c("IL10", "IL19", "IL24", "IL1B", "IL36A", "IL12A", "IL2", "IL21", "IL15", "IL5", "IL4", "IL17F", "IFNG", "IFNL1", "IFNL3", "IFNK", "IFNE", "IFNL4"))
cytokines <- unique(c("IFNL1", "IFNG", "IL10", "IL15", "IL19", "IL2", "IL21", "IL4"))


# 1. Extract pseudotime from CDS
pseudotime_df <- data.frame(
  cell        = colnames(cds),
  pseudotime  = pseudotime(cds),
  cluster     = cds@clusters$UMAP$clusters
)

# 2. Extract normalized expression from Seurat (log-normalized)
expr_mat <- GetAssayData(seurat.fna.cd8, assay = "RNA", slot = "data")[cytokines, ]

# 3. Combine into long-format data frame
expr_df <- as.data.frame(t(as.matrix(expr_mat))) %>%
  tibble::rownames_to_column("cell") %>%
  pivot_longer(-cell, names_to = "gene", values_to = "expression") %>%
  left_join(pseudotime_df, by = "cell") %>%
  filter(is.finite(pseudotime))

ggplot(expr_df, aes(x = pseudotime, y = expression, color = pseudotime)) +
  scale_color_viridis_c(option = "plasma") +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
              color = "black", linewidth = 0.8) +
  facet_wrap(~ gene, scales = "free_y", ncol = 2) +
  theme_classic(base_size = 12)




