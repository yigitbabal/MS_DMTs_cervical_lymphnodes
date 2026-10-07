setwd("/mnt/h3048/grp_laakso/Project_repository/Yigit/MS_FNA_followup")

library(tidyverse)
library(data.table)


combined_TCR <- read_rds("results/TCR/FNA_combined_TCR.rds")
files <- list.files("results/TCR/tcrdist3_dawit/dbScanResults/", pattern = ".csv", full.names = T)


tcrdist.res <- lapply(files, function(file) {
  fread(file, colClasses = c(meta.clone.id = "character"))
})
tcrdist.res <- rbindlist(tcrdist.res)


sample.meta <- fread("sample_meta_data.csv")
sample.meta$ID
unique(tcrdist.res$Studysubject)

tcrdist.res$sample <- tcrdist.res$Studysubject
tcrdist.res$group = sample.meta$DMT[match(tcrdist.res$sample, sample.meta$ID)]

tcrdist.res$antigen_label <- paste(tcrdist.res$antigen.gene, tcrdist.res$antigen.species, sep = "-")


#expanded_all %>%
#  group_by(sample) %>%
#  summarise(total_tcr = n_distinct(cdr3_aa2), .groups = "drop") %>%
#  as.data.frame() -> total_clone_size
  
lapply(combined_TCR, function(x){
  x %>%
    group_by(sample) %>%
    summarise(total_tcr = n_distinct(CTstrict), .groups = "drop") %>%
    as.data.frame()
}) -> total_clone_size_list

total_clone_size_list <- rbindlist(total_clone_size_list)


plot_data <- tcrdist.res %>%
  filter(antigen.species == "EBV") %>%
  group_by(group, sample, antigen_label) %>%
  summarise(n_clones = n_distinct(TCR), .groups = "drop") %>%
  left_join(total_clone_size_list, by = "sample") %>%
  mutate(pct_clones = (n_clones / total_tcr) * 100)


tcrdist.res %>%
  group_by(group, sample, antigen_label) %>%
  summarise(n_clones = n_distinct(TCR), .groups = "drop") %>%
  as.data.frame()


#plot_data <- tcrdist.res %>%
#  group_by(group, sample, antigen_label) %>%
#  summarise(n_clones = n_distinct(cdr3_b_aa), .groups = "drop")

TOP_N <- 10

top_antigens <- plot_data %>%
  group_by(antigen_label) %>%
  summarise(total = sum(n_clones), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_max(order_by = total,n = TOP_N) %>%
  pull(antigen_label)

plot_data <- plot_data %>%
  mutate(antigen_label = if_else(antigen_label %in% top_antigens,
                                 antigen_label, "Other"),
         antigen_label = factor(antigen_label,
                                levels = c(top_antigens, "Other")))

sample_order <- sample.meta %>%
  arrange(DMT, ID) %>%
  filter(ID %in% unique(tcrdist.res$sample)) %>%
  pull(ID)

plot_data <- plot_data %>%
  mutate(sample = factor(sample, levels = sample_order),
         group  = factor(group,  levels = unique(sample.meta$DMT)))

plot_data$group <- factor(plot_data$group, levels = c("Control", "BeforeTreatment", "antiCD20", "DMF", "NTZ"))


# -----------------------------------------------------------
# STEP 5: Plot — let ggplot assign colours from results
# -----------------------------------------------------------
stacked_plot <- ggplot(
  plot_data,
  aes(x = sample, y = pct_clones, fill = antigen_label)
) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
  facet_grid(
    cols   = vars(group),
    scales = "free_x",
    space  = "free_x"
  ) +
  scale_fill_discrete(name = "Predicted epitope\n(gene · species)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                     breaks = scales::pretty_breaks(n = 6)) +
  labs(
    title = "Predicted specificity in expanded Memory CD8 clonotypes",
    x     = NULL,
    y     = "% of matched clones within sample"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background   = element_rect(fill = "grey92", colour = NA),
    strip.text         = element_text(face = "bold", size = 11),
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y        = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.title       = element_text(face = "bold", size = 10),
    legend.text        = element_text(size = 9),
    legend.key.size    = unit(0.45, "cm"),
    plot.title         = element_text(face = "bold", size = 12)
  )

stacked_plot



