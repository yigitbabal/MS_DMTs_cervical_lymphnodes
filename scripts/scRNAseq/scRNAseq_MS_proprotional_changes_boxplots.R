setwd("P:/grp_laakso/Project_repository/Yigit/MS_FNA_followup")

library(tidyverse)
library(Seurat)
library(data.table)
library(ggpubr)

source("scripts/00_functions.R")

seurat.fna <- read_rds("seurat_objects/FNA/seurat_fna_curated_with_subset.rds")
seurat.pbmc <- read_rds("seurat_objects/PBMC/seurat_pbmc_curated_with_subset.rds")


#########################################
####################### FNA ###########
#########################################
unique(seurat.fna$Celltypes_curated_detailed_withsub)

percentages.fna <- calculate_celltype_percentages(seurat_obj = seurat.fna, sample_col = "sample", celltype_col = "Celltypes_curated_detailed_withsub")

sample.meta <- fread("sample_meta_data.csv")

percentages.fna$group <- sample.meta$Group[match(percentages.fna$sample, sample.meta$ID)]
percentages.fna$DMT <- sample.meta$DMT[match(percentages.fna$sample, sample.meta$ID)]
percentages.fna$DMT_detail <- sample.meta$DMT_detail[match(percentages.fna$sample, sample.meta$ID)]
percentages.fna$DMT <- factor(percentages.fna$DMT, levels = c("Control","BeforeTreatment", "antiCD20", "DMF", "NTZ"))

fna.celltypes <- unique(percentages.fna$Celltypes_curated_detailed_withsub)

plot.percent.fna_with_wilcox <- lapply(fna.celltypes, function(x){
  data = percentages.fna %>% filter(Celltypes_curated_detailed_withsub == x)
  p <- ggplot(data, aes(x = DMT, y = Percentage, fill = DMT)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 3) +
    ggsignif::geom_signif(comparisons = list(c("BeforeTreatment", "Control"), c("BeforeTreatment", "antiCD20"), c("BeforeTreatment", "DMF")), 
                          map_signif_level = TRUE, test = "wilcox.test", color = "black", step_increase = 0.1) +
    labs(title = x) +
    xlab("") + ylab("Percentage of cell type within sample") +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) + 
    theme(legend.position = "none")
  return(p)
})

names(plot.percent.fna_with_wilcox) <- fna.celltypes

plot.percent.fna_with_wilcox$`Tfh-like`

write_rds(plot.percent.fna_with_wilcox, file = "figures/proportional_changes/plot.percent.fna_with_wilcox.rds")

plot.percent.fna <- lapply(fna.celltypes, function(x){
  data = percentages.fna %>% filter(Celltypes_curated_detailed_withsub == x)
  p <- ggplot(data, aes(x = DMT, y = Percentage, fill = DMT)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 3) +
    labs(title = x) +
    xlab("") + ylab("Percentage of cell type within sample") +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) + 
    theme(legend.position = "none")
  return(p)
})

names(plot.percent.fna) <- fna.celltypes

write_rds(plot.percent.fna, file = "figures/proportional_changes/plot.percent.fna.rds")

# save all plots as png
# Iterate over the names of your plot list
lapply(names(plot.percent.fna), function(cell_type) {
  
  safe_filename <- gsub("[/ ]", "_", cell_type)
  
  ggsave(filename = paste0("figures/proportional_changes/percentages_FNA/", paste0(safe_filename, ".png")),
         plot = plot.percent.fna[[cell_type]],
         width = 20, height = 20, 
         dpi = 300, units = "cm",
         bg = "white")
})

plot.percent.fna$`Naive B`



#########################################
####################### PBMC ###########
#########################################

unique(seurat.pbmc$Celltypes_curated_detailed_withsub)

percentages.pbmc <- calculate_celltype_percentages(seurat_obj = seurat.pbmc, sample_col = "sample", celltype_col = "Celltypes_curated_detailed_withsub")

sample.meta <- fread("sample_meta_data.csv")

percentages.pbmc$DMT <- sample.meta$DMT[match(percentages.pbmc$sample, sample.meta$PBMC_ID)]
percentages.pbmc$DMT <- factor(percentages.pbmc$DMT, levels = c("BeforeTreatment", "antiCD20", "NTZ"))

pbmc.celltypes <- unique(percentages.pbmc$Celltypes_curated_detailed_withsub)

plot.percent.fna_with_wilcox <- lapply(pbmc.celltypes, function(x){
  data = percentages.pbmc %>% filter(Celltypes_curated_detailed_withsub == x)
  p <- ggplot(data, aes(x = DMT, y = Percentage, fill = DMT)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 3) +
    ggsignif::geom_signif(comparisons = list(c("BeforeTreatment", "Control"), c("BeforeTreatment", "antiCD20"), c("BeforeTreatment", "DMF")), 
                          map_signif_level = TRUE, test = "wilcox.test", color = "black", step_increase = 0.1) +
    labs(title = x) +
    xlab("") + ylab("Percentage of cell type within sample") +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) + 
    theme(legend.position = "none")
  return(p)
})

names(plot.percent.fna_with_wilcox) <- pbmc.celltypes



write_rds(plot.percent.fna_with_wilcox, file = "figures/proportional_changes/plot.percent.pbmc_with_wilcox.rds")

plot.percent.pbmc <- lapply(pbmc.celltypes, function(x){
  data = percentages.pbmc %>% filter(Celltypes_curated_detailed_withsub == x)
  p <- ggplot(data, aes(x = DMT, y = Percentage, fill = DMT)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 3) +
    labs(title = x) +
    xlab("") + ylab("Percentage of cell type within sample") +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) + 
    theme(legend.position = "none")
  return(p)
})

names(plot.percent.pbmc) <- pbmc.celltypes

write_rds(plot.percent.pbmc, file = "figures/proportional_changes/plot.percent.pbmc.rds")

# save all plots as png
# Iterate over the names of your plot list
lapply(names(plot.percent.pbmc), function(cell_type) {
  
  safe_filename <- gsub("[/ ]", "_", cell_type)
  
  ggsave(filename = paste0("figures/proportional_changes/percentages_PBMC/", paste0(safe_filename, ".png")),
         plot = plot.percent.pbmc[[cell_type]],
         width = 20, height = 20, 
         dpi = 300, units = "cm",
         bg = "white")
})

plot.percent.pbmc$`Naive B`




#### B cell subset with total B cell percentage

b.cell.fna <- c("DN","GC B cells (dark zone)", "Plasmablasts",
                "Tbet+CD11c+","GC B cells (light zone)",
                "USM", "Pre-B", "Naive B", "SM", "IFN-stimulated B")


unique(percentages.fna$Celltypes_curated_detailed_withsub)

percentages.fna.bcell <- percentages.fna %>%
  filter(Celltypes_curated_detailed_withsub %in% b.cell.fna)

percentages.fna.bcell %>%
  group_by(sample) %>%
  summarise(total_bcell = sum(Cell_Count)) -> total.b.fna

percentages.fna.bcell$Total_Cells_in_Sample <- total.b.fna$total_bcell[match(percentages.fna.bcell$sample, total.b.fna$sample)]
percentages.fna.bcell$Percentage <- (percentages.fna.bcell$Cell_Count / percentages.fna.bcell$Total_Cells_in_Sample)*100

percentages.fna.bcell <- percentages.fna.bcell[, c(-6,-8)]
percentages.fna.bcell$tissue <- "FNA"
percentages.fna.bcell$tissue_DMT <- paste(percentages.fna.bcell$tissue, percentages.fna.bcell$DMT, sep = "_")

unique(percentages.pbmc$Celltypes_curated_detailed_withsub)

b.cell.pbmc <- c("SM", "USM", "Naive B", "Plasmablasts", "Tbet+CD11c+")

percentages.pbmc.bcell <- percentages.pbmc %>%
  filter(Celltypes_curated_detailed_withsub %in% b.cell.pbmc)

percentages.pbmc.bcell %>%
  group_by(sample) %>%
  summarise(total_bcell = sum(Cell_Count)) -> total.b.pbmc

percentages.pbmc.bcell$Total_Cells_in_Sample <- total.b.pbmc$total_bcell[match(percentages.pbmc.bcell$sample, total.b.pbmc$sample)]
percentages.pbmc.bcell$Percentage <- (percentages.pbmc.bcell$Cell_Count / percentages.pbmc.bcell$Total_Cells_in_Sample)*100
percentages.pbmc.bcell$tissue <- "PBMC"
percentages.pbmc.bcell$tissue_DMT <- paste(percentages.pbmc.bcell$tissue, percentages.pbmc.bcell$DMT, sep = "_")


percentage.bcell <- rbind(percentages.fna.bcell, percentages.pbmc.bcell)
percentage.bcell$tissue_DMT <- factor(percentage.bcell$tissue_DMT, levels = c("FNA_Control", "FNA_BeforeTreatment", "PBMC_BeforeTreatment", "FNA_antiCD20", "PBMC_antiCD20", "FNA_DMF"))

write_rds(percentage.bcell, file = "results/proportional_changes/percentage_fna_pbmc_bcell.rds")
percentage.bcell <- read_rds("results/proportional_changes/percentage_fna_pbmc_bcell.rds")

percentage.bcell$tissue_DMT <- paste(percentage.bcell$tissue, percentage.bcell$DMT, sep = "_")
percentage.bcell$tissue_DMT <- factor(percentage.bcell$tissue_DMT, levels = c("FNA_Control", "FNA_BeforeTreatment", "PBMC_BeforeTreatment", "FNA_antiCD20", "PBMC_antiCD20", "FNA_DMF"))

percentage.bcell$patient_id <- gsub(pattern = "_FU", replacement = "", percentage.bcell$sample)
percentage.bcell$patient_id <- gsub(pattern = "-PBMC", replacement = "", percentage.bcell$patient_id)
percentage.bcell$patient_id <- gsub(pattern = "-FU", replacement = "", percentage.bcell$patient_id)
percentage.bcell$patient_id <- gsub(pattern = "FU", replacement = "", percentage.bcell$patient_id)
percentage.bcell$patient_tissue <- paste(percentage.bcell$patient_id, percentage.bcell$tissue, sep = "_")

detail <- c(sample.meta$DMT_detail[match(percentage.bcell[percentage.bcell$tissue == "FNA",]$sample, sample.meta$ID)],
sample.meta$DMT_detail[match(percentage.bcell[percentage.bcell$tissue == "PBMC",]$sample, sample.meta$PBMC_ID)])
percentage.bcell$dmt_detail <- detail

percentage.bcell$relapse <- "Non-relapsing"
percentage.bcell[percentage.bcell$patient_id == "MS022", "relapse"] <- "Relapse-during-DMT" 


sample.meta <- fread("sample_meta_data.csv")


percentage.bcell %>%
  filter(DMT %in% c("Control", "BeforeTreatment", "antiCD20")) %>%
  ggplot(aes(x = Celltypes_curated_detailed_withsub, y = Percentage, fill = tissue_DMT, shape = relapse)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_dodge(width = 0.75)) +
  xlab("") + ylab("% of B cells") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  facet_wrap(~Celltypes_curated_detailed_withsub, scales = "free", ncol = 5)

plot_df <- percentage.bcell %>%
  filter(DMT %in% c("Control", "BeforeTreatment", "antiCD20"),
         !is.na(tissue_DMT))

tissue_levels <- sort(unique(plot_df$tissue_DMT))
n_lvls <- length(tissue_levels)

offset_lookup <- tibble(
  tissue_DMT = tissue_levels,
  x_offset   = seq(-0.75/2 + 0.75/(2*n_lvls), 0.75/2 - 0.75/(2*n_lvls), length.out = n_lvls)
)

plot_df <- plot_df %>%
  left_join(offset_lookup, by = "tissue_DMT") %>%
  mutate(x_pos = 1 + x_offset)

ggplot(plot_df, aes(x = x_pos, y = Percentage)) +
  geom_boxplot(aes(group = tissue_DMT, fill = tissue_DMT),
               width = 0.75 / n_lvls * 0.8, outlier.shape = NA) +
  geom_line(aes(group = patient_tissue), color = "grey50", alpha = 0.6, linewidth = 0.4) +
  geom_point(aes(color = dmt_detail, shape = relapse), size = 2) +
  scale_color_manual(values = c(
    "Control" = "#8C8C8C",
    "BeforeTreatment" = "#4D4D4D",
    "antiCD20 (rituximab)" = "#0073C2",
    "antiCD20 (ofatumumab)" = "#EFC000"
  )) +
  scale_x_continuous(breaks = NULL) +
  xlab("") + ylab("% of B cells") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  facet_wrap(~Celltypes_curated_detailed_withsub, scales = "free", ncol = 5)







#### fold chanhges
#FNA
percentages.fna.bcell %>%
  group_by(Celltypes_curated_detailed_withsub, DMT) %>%
  summarise(mean = mean(Percentage)) -> mean.percentage.fna

fc.comp <- list(c("BeforeTreatment", "Control"),
                c("antiCD20", "BeforeTreatment"),
                c("DMF", "BeforeTreatment"),
                c("NTZ", "BeforeTreatment"))

df_log2fc.fna <- mean.percentage.fna %>%
  # Pivot so that each DMT level is its own column
  pivot_wider(names_from = DMT, values_from = mean) %>%
  
  # Calculate log2FC for each pair: log2(Condition_1 / Condition_2)
  mutate(
    log2FC_BeforeTreatment_vs_Control = round(log2(BeforeTreatment / Control), digits = 2),
    log2FC_antiCD20_vs_BeforeTreatment = round(log2(antiCD20 / BeforeTreatment), digits = 2),
    log2FC_DMF_vs_BeforeTreatment      = round(log2(DMF / BeforeTreatment), digits = 2),
    log2FC_NTZ_vs_BeforeTreatment      = round(log2(NTZ / BeforeTreatment), digits = 2)
  )

# View the wide format result
print(df_log2fc.fna)

#PBMC
percentages.pbmc.bcell %>%
  group_by(Celltypes_curated_detailed_withsub, DMT) %>%
  summarise(mean = mean(Percentage)) -> mean.percentage.pbmc


df_log2fc.pbmc <- mean.percentage.pbmc %>%
  # Pivot so that each DMT level is its own column
  pivot_wider(names_from = DMT, values_from = mean) %>%
  
  # Calculate log2FC for each pair: log2(Condition_1 / Condition_2)
  mutate(
    log2FC_antiCD20_vs_BeforeTreatment = round(log2(antiCD20 / BeforeTreatment), digits = 2),
    log2FC_NTZ_vs_BeforeTreatment      = round(log2(NTZ / BeforeTreatment), digits = 2)
  )

# View the wide format result
print(df_log2fc.pbmc)




## relative percentage changes
#### Relative percentage changes
#FNA
percentages.fna.bcell %>%
  group_by(Celltypes_curated_detailed_withsub, DMT) %>%
  summarise(mean = mean(Percentage)) -> mean.percentage.fna

# quick check for zero baselines before dividing
mean.percentage.fna %>% filter(mean == 0)

df_pctchange.fna <- mean.percentage.fna %>%
  # Pivot so that each DMT level is its own column
  pivot_wider(names_from = DMT, values_from = mean) %>%
  
  # Relative percent change: (new - old) / old * 100
  mutate(
    pctChange_BeforeTreatment_vs_Control  = round((BeforeTreatment - Control) / Control * 100, digits = 2),
    pctChange_antiCD20_vs_BeforeTreatment = round((antiCD20 - BeforeTreatment) / BeforeTreatment * 100, digits = 2),
    pctChange_DMF_vs_BeforeTreatment      = round((DMF - BeforeTreatment) / BeforeTreatment * 100, digits = 2),
    pctChange_NTZ_vs_BeforeTreatment      = round((NTZ - BeforeTreatment) / BeforeTreatment * 100, digits = 2)
  )

# View the wide format result
print(df_pctchange.fna)

#PBMC
percentages.pbmc.bcell %>%
  group_by(Celltypes_curated_detailed_withsub, DMT) %>%
  summarise(mean = mean(Percentage)) -> mean.percentage.pbmc

# quick check for zero baselines before dividing
mean.percentage.pbmc %>% filter(mean == 0)

df_pctchange.pbmc <- mean.percentage.pbmc %>%
  # Pivot so that each DMT level is its own column
  pivot_wider(names_from = DMT, values_from = mean) %>%
  
  # Relative percent change: (new - old) / old * 100
  mutate(
    pctChange_antiCD20_vs_BeforeTreatment = round((antiCD20 - BeforeTreatment) / BeforeTreatment * 100, digits = 2),
    pctChange_NTZ_vs_BeforeTreatment      = round((NTZ - BeforeTreatment) / BeforeTreatment * 100, digits = 2)
  )

# View the wide format result
print(df_pctchange.pbmc)


fwrite(df_pctchange.fna, file = "figure_updates/df_pctchange_fna.csv", row.names = F)
fwrite(df_pctchange.pbmc, file = "figure_updates/df_pctchange_pbmc.csv", row.names = F)
