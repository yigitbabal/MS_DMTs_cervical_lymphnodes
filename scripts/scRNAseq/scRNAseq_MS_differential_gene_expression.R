setwd("P:/grp_laakso/Project_repository/Yigit/MS_FNA_followup")

library(tidyverse)
library(Seurat)
library(data.table)

source("scripts/00_functions.R")

seurat.fna.bcell <- read_rds("seurat_objects/FNA/seurat_fna_bcell_subset.rds")
seurat.pbmc.bcell <- read_rds("seurat_objects/PBMC/seurat_pbmc_bcell_subset.rds")


# For FNA
seurat.fna <- read_rds("seurat_objects/FNA/seurat_fna_integrated_singleR_annotated_curated.rds")

seurat.fna$patient <- str_split(seurat.fna$sample, pattern = "_", simplify = T)[,1]
unique(seurat.fna$patient)


seurat.fna <- subset(seurat.fna, subset = Celltypes_curated_detailed != "Low quality/stressed T cells" & 
                       Celltypes_curated_detailed != "Stressed/ribo high T cells")

bcell.subset.id <- data.frame(cellid = seurat.fna.bcell$cellid,
                              Bsubclusters = seurat.fna.bcell$Bsubclusters)

seurat.fna$cellid <- colnames(seurat.fna)
seurat.fna$cellid <- gsub(pattern = "-", replacement = "_", x = seurat.fna$cellid)

seurat.fna$Celltypes_curated_detailed_withsub <- as.character(seurat.fna$Celltypes_curated_detailed)
seurat.fna@meta.data[seurat.fna@meta.data$cellid %in% bcell.subset.id$cellid, "Celltypes_curated_detailed_withsub"] <- bcell.subset.id$Bsubclusters

unique(seurat.fna$Celltypes_curated_detailed_withsub)

seurat.fna <- subset(seurat.fna, subset = Celltypes_curated_detailed_withsub != "Naive B cells" & 
                       Celltypes_curated_detailed_withsub != "Memory B cells" &
                       Celltypes_curated_detailed_withsub != "Pre B cells")

unique(seurat.fna$Celltypes_curated_detailed_withsub)

deg.fna.list <- RunConsensusPseudobulk(seurat_obj = seurat.fna, 
                                       cell_type_col = "Celltypes_curated_detailed_withsub", 
                                       treatment_col = "DMT",
                                       sample_col = "sample", 
                                       baseline_group = "BeforeTreatment", 
                                       after_groups = c("antiCD20", "DMF"),
                                       p_val_cutoff = 0.05)

write_rds(deg.fna.list, "results/DEGs/deg_fna_consensus_pseudobulk.rds")


deg.fna.wilcox.list <- RunSeuratWilcoxon(seurat_obj = seurat.fna, 
                                     cell_type_col = "Celltypes_curated_detailed_withsub", 
                                     treatment_col = "DMT",
                                     baseline_group = "BeforeTreatment", 
                                     after_groups = c("antiCD20", "DMF"), 
                                     min.pct = 0.25,
                                     p_val_cutoff = 0.05)

write_rds(deg.fna.wilcox.list, "results/DEGs/deg_fna_wilcox.rds")


# viz

library(ggplot2)
library(dplyr)
library(ggrepel)
library(scales) 

# ---------------------------------------------------------
# 1. Prepare Base Plot Data
# ---------------------------------------------------------
plot_df <- deg.fna.list$significant_consensus %>%
  mutate(
    sig_dir_DESeq2 = sign(logFC_DESeq2) * sig_DESeq2,
    sig_dir_edgeR  = sign(logFC_edgeR) * sig_edgeR,
    sig_dir_limma  = sign(logFC_limma) * sig_limma
  ) %>%
  filter(abs(sig_dir_DESeq2 + sig_dir_edgeR + sig_dir_limma) == consensus_score) %>%
  filter(comparison == "antiCD20_vs_BeforeTreatment") %>%
  mutate(
    avg_logFC = ((logFC_DESeq2 * sig_DESeq2) + 
                   (logFC_edgeR * sig_edgeR) + 
                   (logFC_limma * sig_limma)) / consensus_score,
    Consensus = as.factor(consensus_score),
    Direction = ifelse(avg_logFC > 0, "Up", "Down") 
  )

# ---------------------------------------------------------
# 2. Extract Top 10 Genes PER CELL TYPE
# ---------------------------------------------------------
top_genes_df <- plot_df %>%
  group_by(comparison, cell_type, Direction) %>%
  slice_max(order_by = abs(avg_logFC), n = 10) %>%
  ungroup()

# ---------------------------------------------------------
# 3. Calculate Gene Counts per Row/Facet
# ---------------------------------------------------------
counts_df <- plot_df %>%
  group_by(comparison, cell_type) %>%
  summarise(
    Up_Count = sum(Direction == "Up"),
    Down_Count = sum(Direction == "Down"),
    .groups = "drop"
  )

# ---------------------------------------------------------
# 4. CUSTOM COLOR MAPPING LOGIC
# ---------------------------------------------------------
# Find the absolute minimum and maximum log2FC in your actual data
min_fc <- min(plot_df$avg_logFC, na.rm = TRUE)
max_fc <- max(plot_df$avg_logFC, na.rm = TRUE)

# Define the anchor points for our colors: 
# [True Min] -> [-2] -> [0] -> [+2] -> [True Max]
# Note: The max() and min() functions inside ensure it doesn't break if your data doesn't naturally reach 2 or -2
color_breaks <- c(min_fc, max(min_fc, -2), 0, min(max_fc, 2), max_fc)

# Rescale those anchor points to a 0-to-1 scale (which ggplot requires for 'values')
color_values <- scales::rescale(color_breaks)

# ---------------------------------------------------------
# 5. Generate the Plot
# ---------------------------------------------------------
p <- ggplot(plot_df, aes(x = avg_logFC, y = cell_type)) +
  
  geom_jitter(aes(color = avg_logFC, size = Consensus), height = 0.2, alpha = 0.8) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  
  geom_text(data = counts_df, 
            aes(x = -Inf, y = cell_type, label = paste0(" \u2193 ", Down_Count)), 
            hjust = 0, vjust = 0.5, size = 3.5, color = "darkblue", fontface = "bold") +
  
  geom_text(data = counts_df, 
            aes(x = Inf, y = cell_type, label = paste0(Up_Count, " \u2191 ")), 
            hjust = 1, vjust = 0.5, size = 3.5, color = "darkred", fontface = "bold") +
  
  geom_text_repel(data = top_genes_df, 
                  aes(label = gene), 
                  size = 2.5,
                  color = "black",
                  fontface = "italic",
                  box.padding = 0.3,
                  point.padding = 0.2,
                  max.overlaps = Inf) + 
  
  facet_grid(. ~ comparison, scales = "free_y") +
  
  # NEW: Apply the custom gradient mapping
  # It assigns Dark Blue to everything from Min to -2, Grey to 0, and Dark Red to everything from 2 to Max.
  scale_color_gradientn(
    colors = c("darkblue", "darkblue", "grey85", "darkred", "darkred"),
    values = color_values,
    limits = c(min_fc, max_fc) # Ensures the legend spans your true data range
  ) +
  
  scale_size_manual(values = c("2" = 1.2, "3" = 2.5)) + 
  scale_x_continuous(expand = expansion(mult = 0.3)) +
  
  theme_bw() +
  labs(title = "dcLN - Differential Expression Shifts by Cell Type",
       subtitle = "Consensus Significant Genes (DESeq2, EdgeR and limma)",
       x = "Average Log2 Fold Change (Significant Tools Only)",
       y = "",
       color = "Avg Log2FC",
       size = "Tools Agreed") +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 11), 
    axis.text.y = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(p)


## go analysis
library(clusterProfiler)
library(org.Hs.eg.db) # CHANGE to org.Mm.eg.db if using mouse data
library(dplyr)
library(tibble)

# ---------------------------------------------------------
# 1. Prepare the Data (Using the same strict filters)
# ---------------------------------------------------------
# We start with your strictly filtered significant consensus genes
go_input_df <- deg.fna.list$significant_consensus %>%
  mutate(
    sig_dir_DESeq2 = sign(logFC_DESeq2) * sig_DESeq2,
    sig_dir_edgeR  = sign(logFC_edgeR) * sig_edgeR,
    sig_dir_limma  = sign(logFC_limma) * sig_limma
  ) %>%
  filter(abs(sig_dir_DESeq2 + sig_dir_edgeR + sig_dir_limma) == consensus_score) %>%
  filter(comparison == "antiCD20_vs_BeforeTreatment") %>%
  mutate(
    avg_logFC = ((logFC_DESeq2 * sig_DESeq2) + 
                   (logFC_edgeR * sig_edgeR) + 
                   (logFC_limma * sig_limma)) / consensus_score,
    Direction = ifelse(avg_logFC > 0, "Up", "Down") 
  )

# Get unique combinations of Cell Types and Comparisons to loop through
comparisons <- unique(go_input_df$comparison)
cell_types <- unique(go_input_df$cell_type)
directions <- c("Up", "Down")

all_go_results <- list()

# ---------------------------------------------------------
# 2. RUN CLUSTERPROFILER LOOP
# ---------------------------------------------------------
message("Starting Gene Ontology Enrichment Analysis...")

for (comp in comparisons) {
  for (ct in cell_types) {
    for (dir in directions) {
      
      # Extract the specific list of genes for this exact scenario
      target_genes <- go_input_df %>%
        filter(comparison == comp, cell_type == ct, Direction == dir) %>%
        pull(gene)
      
      # clusterProfiler usually needs at least 5-10 genes to find meaningful pathways
      if (length(target_genes) > 5) {
        
        message(paste("Running GO for:", comp, "|", ct, "|", dir, paste0("(", length(target_genes), " genes)")))
        
        # Run the enrichment test
        go_res <- enrichGO(
          gene          = target_genes,
          OrgDb         = org.Hs.eg.db, # CHANGE THIS if using mouse
          keyType       = "SYMBOL",     # Tells it we are providing Gene Symbols
          ont           = "BP",         # BP = Biological Process (usually most informative)
          pAdjustMethod = "BH",         # Benjamini-Hochberg FDR correction
          pvalueCutoff  = 0.05,
          qvalueCutoff  = 0.05,
          readable      = TRUE          # Ensures output gene lists are symbols, not Entrez IDs
        )
        
        # If the test found significant pathways, save them!
        if (!is.null(go_res) && nrow(go_res@result %>% filter(p.adjust < 0.05)) > 0) {
          
          # Convert the S4 object to a normal dataframe and add metadata
          res_df <- as.data.frame(go_res) %>%
            mutate(
              comparison = comp,
              cell_type = ct,
              Direction = dir
            )
          
          # Create a unique name for the list
          list_name <- paste(comp, ct, dir, sep = "_")
          all_go_results[[list_name]] <- res_df
          
        } else {
          message("  -> No significant GO terms found.")
        }
        
      } else {
        message(paste("Skipping:", comp, "|", ct, "|", dir, "- Not enough genes for GO analysis."))
      }
    }
  }
}

# ---------------------------------------------------------
# 3. COMBINE FINAL RESULTS
# ---------------------------------------------------------
final_go_table <- bind_rows(all_go_results)

message("GO Analysis Complete!")

# View the top pathways found
head(final_go_table)

write_rds(final_go_table, file = "results/DEGs/FNA_GO_enrichment_antiCD20_results.rds")
#write_rds(final_go_table, file = "results/DEGs/FNA_GO_enrichment_DMF_results.rds")

final_go_table <- read_rds("results/DEGs/FNA_GO_enrichment_antiCD20_results.rds")

final_go_table <- final_go_table[final_go_table$cell_type == "Naive B",]

cols = c("Down" = "darkblue", "Up" = "darkred")
p.fna.go <- final_go_table %>%
  ggplot(aes(x = -log10(p.adjust), y = reorder(Description, p.adjust, decreasing = T), fill = Direction)) +
  geom_bar(stat = "identity") + theme_classic() + scale_fill_manual(values = cols) +
  ggtitle(label = "dcLN - Naive B - GO terms") + ylab("") + xlab("-log10(adjusted p-value)")



## for PBMC
seurat.pbmc <- read_rds("seurat_objects/PBMC/seurat_pbmc_integrated_singleR_annotated_curated.rds")
seurat.pbmc$patient <- str_split(seurat.pbmc$sample, pattern = "-", simplify = T)[,1]
seurat.pbmc$patient <- gsub(pattern = "FU", replacement = "", x = seurat.pbmc$patient)
seurat.pbmc$cellid <- colnames(seurat.pbmc)
seurat.pbmc$cellid <- gsub(pattern = "-", replacement = "_", x = seurat.pbmc$cellid)

bcell.subset.id <- data.frame(cellid = seurat.pbmc.bcell$cellid,
                              Bsubclusters = seurat.pbmc.bcell$Bsubclusters)

seurat.pbmc$Celltypes_curated_detailed_withsub <- as.character(seurat.pbmc$Celltypes_curated_detailed)
seurat.pbmc@meta.data[seurat.pbmc@meta.data$cellid %in% bcell.subset.id$cellid, "Celltypes_curated_detailed_withsub"] <- bcell.subset.id$Bsubclusters

unique(seurat.pbmc$Celltypes_curated_detailed_withsub)

seurat.pbmc <- subset(seurat.pbmc, subset = Celltypes_curated_detailed_withsub != "Naive B cells" & 
                       Celltypes_curated_detailed_withsub != "Memory B cells")

unique(seurat.pbmc$Celltypes_curated_detailed_withsub)
DimPlot(seurat.pbmc, group.by = "Celltypes_curated_detailed_withsub", label = T)


deg.pbmc.list <- RunConsensusPseudobulk(seurat_obj = seurat.pbmc, 
                                       cell_type_col = "Celltypes_curated_detailed_withsub", 
                                       treatment_col = "DMT",
                                       sample_col = "sample", 
                                       baseline_group = "BeforeTreatment", 
                                       after_groups = c("antiCD20"),
                                       p_val_cutoff = 0.05)

write_rds(deg.pbmc.list, "results/DEGs/deg_pbmc_consensus_pseudobulk.rds")


deg.pbmc.wilcox.list <- RunSeuratWilcoxon(seurat_obj = seurat.pbmc, 
                                         cell_type_col = "Celltypes_curated_detailed_withsub", 
                                         treatment_col = "DMT",
                                         baseline_group = "BeforeTreatment", 
                                         after_groups = c("antiCD20"), 
                                         min.pct = 0.25,
                                         p_val_cutoff = 0.05)

write_rds(deg.pbmc.wilcox.list, "results/DEGs/deg_pbmc_wilcox.rds")

# ---------------------------------------------------------
# 1. Prepare Base Plot Data
# ---------------------------------------------------------
plot_df <- deg.pbmc.list$significant_consensus %>%
  mutate(
    sig_dir_DESeq2 = sign(logFC_DESeq2) * sig_DESeq2,
    sig_dir_edgeR  = sign(logFC_edgeR) * sig_edgeR,
    sig_dir_limma  = sign(logFC_limma) * sig_limma
  ) %>%
  filter(abs(sig_dir_DESeq2 + sig_dir_edgeR + sig_dir_limma) == consensus_score) %>%
  filter(comparison == "antiCD20_vs_BeforeTreatment") %>%
  mutate(
    avg_logFC = ((logFC_DESeq2 * sig_DESeq2) + 
                   (logFC_edgeR * sig_edgeR) + 
                   (logFC_limma * sig_limma)) / consensus_score,
    Consensus = as.factor(consensus_score),
    Direction = ifelse(avg_logFC > 0, "Up", "Down") 
  )

# ---------------------------------------------------------
# 2. Extract Top 10 Genes PER CELL TYPE
# ---------------------------------------------------------
top_genes_df <- plot_df %>%
  group_by(comparison, cell_type, Direction) %>%
  slice_max(order_by = abs(avg_logFC), n = 10) %>%
  ungroup()

# ---------------------------------------------------------
# 3. Calculate Gene Counts per Row/Facet
# ---------------------------------------------------------
counts_df <- plot_df %>%
  group_by(comparison, cell_type) %>%
  summarise(
    Up_Count = sum(Direction == "Up"),
    Down_Count = sum(Direction == "Down"),
    .groups = "drop"
  )

# ---------------------------------------------------------
# 4. CUSTOM COLOR MAPPING LOGIC
# ---------------------------------------------------------
# Find the absolute minimum and maximum log2FC in your actual data
min_fc <- min(plot_df$avg_logFC, na.rm = TRUE)
max_fc <- max(plot_df$avg_logFC, na.rm = TRUE)

# Define the anchor points for our colors: 
# [True Min] -> [-2] -> [0] -> [+2] -> [True Max]
# Note: The max() and min() functions inside ensure it doesn't break if your data doesn't naturally reach 2 or -2
color_breaks <- c(min_fc, max(min_fc, -2), 0, min(max_fc, 2), max_fc)

# Rescale those anchor points to a 0-to-1 scale (which ggplot requires for 'values')
color_values <- scales::rescale(color_breaks)

# ---------------------------------------------------------
# 5. Generate the Plot
# ---------------------------------------------------------
p2 <- ggplot(plot_df, aes(x = avg_logFC, y = cell_type)) +
  
  geom_jitter(aes(color = avg_logFC, size = Consensus), height = 0.2, alpha = 0.8) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  
  geom_text(data = counts_df, 
            aes(x = -Inf, y = cell_type, label = paste0(" \u2193 ", Down_Count)), 
            hjust = 0, vjust = 0.5, size = 3.5, color = "darkblue", fontface = "bold") +
  
  geom_text(data = counts_df, 
            aes(x = Inf, y = cell_type, label = paste0(Up_Count, " \u2191 ")), 
            hjust = 1, vjust = 0.5, size = 3.5, color = "darkred", fontface = "bold") +
  
  geom_text_repel(data = top_genes_df, 
                  aes(label = gene), 
                  size = 2.5,
                  color = "black",
                  fontface = "italic",
                  box.padding = 0.3,
                  point.padding = 0.2,
                  max.overlaps = Inf) + 
  
  facet_grid(. ~ comparison, scales = "free_y") +
  
  # NEW: Apply the custom gradient mapping
  # It assigns Dark Blue to everything from Min to -2, Grey to 0, and Dark Red to everything from 2 to Max.
  scale_color_gradientn(
    colors = c("darkblue", "darkblue", "grey85", "darkred", "darkred"),
    values = color_values,
    limits = c(min_fc, max_fc) # Ensures the legend spans your true data range
  ) +
  
  scale_size_manual(values = c("2" = 1.2, "3" = 2.5)) + 
  scale_x_continuous(expand = expansion(mult = 0.3)) +
  
  theme_bw() +
  labs(title = "PBMC - Differential Expression Shifts by Cell Type",
       subtitle = "Consensus Significant Genes (DESeq2, EdgeR and limma)",
       x = "Average Log2 Fold Change (Significant Tools Only)",
       y = "",
       color = "Avg Log2FC",
       size = "Tools Agreed") +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 11), 
    axis.text.y = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(p2)

## go analysis
library(clusterProfiler)
library(org.Hs.eg.db) # CHANGE to org.Mm.eg.db if using mouse data
library(dplyr)
library(tibble)

# ---------------------------------------------------------
# 1. Prepare the Data (Using the same strict filters)
# ---------------------------------------------------------
# We start with your strictly filtered significant consensus genes
go_input_df <- deg.pbmc.list$significant_consensus %>%
  mutate(
    sig_dir_DESeq2 = sign(logFC_DESeq2) * sig_DESeq2,
    sig_dir_edgeR  = sign(logFC_edgeR) * sig_edgeR,
    sig_dir_limma  = sign(logFC_limma) * sig_limma
  ) %>%
  filter(abs(sig_dir_DESeq2 + sig_dir_edgeR + sig_dir_limma) == consensus_score) %>%
  filter(comparison == "antiCD20_vs_BeforeTreatment") %>%
  mutate(
    avg_logFC = ((logFC_DESeq2 * sig_DESeq2) + 
                   (logFC_edgeR * sig_edgeR) + 
                   (logFC_limma * sig_limma)) / consensus_score,
    Direction = ifelse(avg_logFC > 0, "Up", "Down") 
  )

# Get unique combinations of Cell Types and Comparisons to loop through
comparisons <- unique(go_input_df$comparison)
cell_types <- unique(go_input_df$cell_type)
directions <- c("Up", "Down")

all_go_results <- list()

# ---------------------------------------------------------
# 2. RUN CLUSTERPROFILER LOOP
# ---------------------------------------------------------
message("Starting Gene Ontology Enrichment Analysis...")

for (comp in comparisons) {
  for (ct in cell_types) {
    for (dir in directions) {
      
      # Extract the specific list of genes for this exact scenario
      target_genes <- go_input_df %>%
        filter(comparison == comp, cell_type == ct, Direction == dir) %>%
        pull(gene)
      
      # clusterProfiler usually needs at least 5-10 genes to find meaningful pathways
      if (length(target_genes) > 5) {
        
        message(paste("Running GO for:", comp, "|", ct, "|", dir, paste0("(", length(target_genes), " genes)")))
        
        # Run the enrichment test
        go_res <- enrichGO(
          gene          = target_genes,
          OrgDb         = org.Hs.eg.db, # CHANGE THIS if using mouse
          keyType       = "SYMBOL",     # Tells it we are providing Gene Symbols
          ont           = "BP",         # BP = Biological Process (usually most informative)
          pAdjustMethod = "BH",         # Benjamini-Hochberg FDR correction
          pvalueCutoff  = 0.05,
          qvalueCutoff  = 0.05,
          readable      = TRUE          # Ensures output gene lists are symbols, not Entrez IDs
        )
        
        # If the test found significant pathways, save them!
        if (!is.null(go_res) && nrow(go_res@result %>% filter(p.adjust < 0.05)) > 0) {
          
          # Convert the S4 object to a normal dataframe and add metadata
          res_df <- as.data.frame(go_res) %>%
            mutate(
              comparison = comp,
              cell_type = ct,
              Direction = dir
            )
          
          # Create a unique name for the list
          list_name <- paste(comp, ct, dir, sep = "_")
          all_go_results[[list_name]] <- res_df
          
        } else {
          message("  -> No significant GO terms found.")
        }
        
      } else {
        message(paste("Skipping:", comp, "|", ct, "|", dir, "- Not enough genes for GO analysis."))
      }
    }
  }
}

# ---------------------------------------------------------
# 3. COMBINE FINAL RESULTS
# ---------------------------------------------------------
final_go_table.pbmc <- bind_rows(all_go_results)

write_rds(final_go_table.pbmc, file = "results/DEGs/PBMC_GO_enrichment_results.rds")

final_go_table.pbmc <- read_rds("results/DEGs/PBMC_GO_enrichment_results.rds")

final_go_table.pbmc <- final_go_table.pbmc[final_go_table.pbmc$cell_type == "Naive B",]
final_go_table.pbmc <- final_go_table.pbmc[1:12,]

cols = c("Down" = "darkblue", "Up" = "darkred")
p.pbmc.go <- final_go_table.pbmc %>%
  ggplot(aes(x = -log10(p.adjust), y = reorder(Description, p.adjust, decreasing = T), fill = Direction)) +
  geom_bar(stat = "identity") + theme_classic() + scale_fill_manual(values = cols) +
  ggtitle(label = "PBMC - Naive B - GO terms") + ylab("") + xlab("-log10(adjusted p-value)")




## venn diagram
library(dplyr)
library(ggVennDiagram)
library(ggplot2)

# ---------------------------------------------------------
# 1. Helper Function: Clean and Prep Data
# ---------------------------------------------------------
# We apply the exact same strict directional filter we used for plotting
prep_deg_data <- function(deg_list, dataset_name) {
  deg_list$significant_consensus %>%
    mutate(
      sig_dir_DESeq2 = sign(logFC_DESeq2) * sig_DESeq2,
      sig_dir_edgeR  = sign(logFC_edgeR) * sig_edgeR,
      sig_dir_limma  = sign(logFC_limma) * sig_limma
    ) %>%
    filter(abs(sig_dir_DESeq2 + sig_dir_edgeR + sig_dir_limma) == consensus_score) %>%
    mutate(
      avg_logFC = ((logFC_DESeq2 * sig_DESeq2) + 
                     (logFC_edgeR * sig_edgeR) + 
                     (logFC_limma * sig_limma)) / consensus_score,
      Direction = ifelse(avg_logFC > 0, "Up", "Down"),
      Dataset = dataset_name
    )
}

# Prep both datasets
fna_df <- prep_deg_data(deg.fna.list, "FNA")
pbmc_df <- prep_deg_data(deg.pbmc.list, "PBMC")

# ---------------------------------------------------------
# 2. Identify Shared Variables
# ---------------------------------------------------------
shared_cells <- intersect(unique(fna_df$cell_type), unique(pbmc_df$cell_type))
shared_comps <- intersect(unique(fna_df$comparison), unique(pbmc_df$comparison))
directions <- c("Up", "Down")

message(paste("Found", length(shared_cells), "shared cell types and", length(shared_comps), "shared comparisons."))

# ---------------------------------------------------------
# 3. Loop and Compare
# ---------------------------------------------------------
all_venn_plots <- list()
all_overlap_data <- list()

for (comp in shared_comps) {
  for (ct in shared_cells) {
    for (dir in directions) {
      
      # Extract gene lists for this specific scenario
      genes_fna <- fna_df %>% filter(comparison == comp, cell_type == ct, Direction == dir) %>% pull(gene)
      genes_pbmc <- pbmc_df %>% filter(comparison == comp, cell_type == ct, Direction == dir) %>% pull(gene)
      
      # Skip if BOTH datasets have 0 genes for this cell type/direction
      if (length(genes_fna) == 0 && length(genes_pbmc) == 0) next
      
      # --- A. Create the Venn Diagram Plot ---
      gene_list_for_venn <- list(
        FNA = genes_fna,
        PBMC = genes_pbmc
      )
      
      # Name the plot based on the variables
      plot_title <- paste(ct, "|", comp, "|", dir, "Regulated")
      plot_name <- paste(ct, comp, dir, sep = "_")
      
      v_plot <- ggVennDiagram(gene_list_for_venn, label_alpha = 0) +
        scale_fill_gradient(low = "white", high = ifelse(dir == "Up", "mistyrose", "aliceblue")) +
        scale_color_manual(values = c("black", "black")) +
        labs(title = plot_title, subtitle = "Significant Consensus Genes") +
        theme(plot.title = element_text(face = "bold", hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5))
      
      all_venn_plots[[plot_name]] <- v_plot
      
      # --- B. Classify the Genes for the Dataframe ---
      shared_genes <- intersect(genes_fna, genes_pbmc)
      fna_only <- setdiff(genes_fna, genes_pbmc)
      pbmc_only <- setdiff(genes_pbmc, genes_fna)
      
      # Build a temporary dataframe for this loop iteration
      temp_df <- data.frame(
        cell_type = ct,
        comparison = comp,
        Direction = dir,
        Gene = c(shared_genes, fna_only, pbmc_only),
        Overlap_Category = c(
          rep("Shared", length(shared_genes)),
          rep("FNA_Only", length(fna_only)),
          rep("PBMC_Only", length(pbmc_only))
        )
      )
      
      # Only add to our list if there are actually genes to report
      if (nrow(temp_df) > 0) {
        all_overlap_data[[plot_name]] <- temp_df
      }
    }
  }
}

# ---------------------------------------------------------
# 4. Finalize the Output Object
# ---------------------------------------------------------
# Bind the dataframes together into one master table
final_overlap_df <- bind_rows(all_overlap_data)

# Create the final list object you requested
cross_tissue_comparison <- list(
  Venn_Plots = all_venn_plots,
  Overlap_Data = final_overlap_df
)

message("Comparison complete! Use cross_tissue_comparison$Venn_Plots and cross_tissue_comparison$Overlap_Data")

cross_tissue_comparison$Venn_Plots$`Naive B_antiCD20_vs_BeforeTreatment_Up` +
cross_tissue_comparison$Venn_Plots$`Naive B_antiCD20_vs_BeforeTreatment_Down`

write_rds(cross_tissue_comparison, file = "results/DEGs/cross_tissue_DEG_comparison.rds")


# save seurat.fna and seurat.pbmc
write_rds(seurat.fna, "seurat_objects/FNA/seurat_fna_curated_with_subset.rds")
write_rds(seurat.pbmc, "seurat_objects/PBMC/seurat_pbmc_curated_with_subset.rds")
