setwd("P:/grp_laakso/Project_repository/Yigit/MS_FNA_followup")

library(tidyverse)
library(Seurat)
library(data.table)

source("scripts/00_functions.R")
## sccomp approach

library(sccomp)

seurat.fna <- read_rds("seurat_objects/FNA/seurat_fna_curated_with_subset.rds")
seurat.pbmc <- read_rds("seurat_objects/PBMC/seurat_pbmc_curated_with_subset.rds")

unique(seurat.fna$Celltypes_curated_detailed_withsub)


sccomp.res = 
  seurat.fna |>
  sccomp_estimate( 
    formula_composition = ~ 0 + DMT +(DMT | Protocol), 
    sample = "sample",
    cell_group = "Celltypes_curated_detailed_withsub",
    bimodal_mean_variability_association = TRUE,
    cores = 5, verbose = TRUE,
    max_sampling_iterations = 2000
  ) 

sccomp.test.control <- sccomp.res %>%
  sccomp_test(contrasts = c("DMTBeforeTreatment - DMTControl"), test_composition_above_logit_fold_change = 0)

sccomp.test.anticd20 <- sccomp.res %>%
  sccomp_test(contrasts = c("DMTantiCD20 - DMTBeforeTreatment"), test_composition_above_logit_fold_change = 0)

sccomp.test.dmf <- sccomp.res %>%
  sccomp_test(contrasts = c("DMTDMF - DMTBeforeTreatment"), test_composition_above_logit_fold_change = 0)

sccomp.test.ntz <- sccomp.res %>%
  sccomp_test(contrasts = c("DMTNTZ - DMTBeforeTreatment"), test_composition_above_logit_fold_change = 0)

sccomp.plot.control <- plot_1D_intervals(sccomp.test.control, test_composition_above_logit_fold_change = 0, significance_statistic = "FDR")
sccomp.plot.anticd20 <- plot_1D_intervals(sccomp.test.anticd20, test_composition_above_logit_fold_change = 0, significance_statistic = "FDR")
sccomp.plot.dmf <- plot_1D_intervals(sccomp.test.dmf, test_composition_above_logit_fold_change = 0, significance_statistic = "FDR")
sccomp.plot.ntz <- plot_1D_intervals(sccomp.test.ntz, test_composition_above_logit_fold_change = 0, significance_statistic = "FDR")

write_rds(list(model = sccomp.res,
               control = list(res = sccomp.test.control,
                               plot = sccomp.plot.control), 
              anticd20 = list(res = sccomp.test.anticd20,
                               plot = sccomp.plot.anticd20),
               DMF = list(res = sccomp.test.dmf,
                          plot = sccomp.plot.dmf),
               NTZ = list(res = sccomp.test.ntz,
                          plot = sccomp.plot.ntz)), "results/proportional_changes/sccomp_fna_results.rds")

### for pbmc

sccomp.res.pbmc = 
  seurat.pbmc |>
  sccomp_estimate( 
    formula_composition = ~ 0 + DMT, 
    sample = "sample",
    cell_group = "Celltypes_curated_detailed_withsub",
    bimodal_mean_variability_association = TRUE,
    cores = 5, verbose = TRUE,
    max_sampling_iterations = 2000
  ) 

sccomp.test.anticd20.pbmc <- sccomp.res.pbmc %>%
  sccomp_test(contrasts = c("DMTantiCD20 - DMTBeforeTreatment"), test_composition_above_logit_fold_change = 0)

sccomp.test.ntz.pbmc <- sccomp.res.pbmc %>%
  sccomp_test(contrasts = c("DMTNTZ - DMTBeforeTreatment"), test_composition_above_logit_fold_change = 0)

sccomp.plot.anticd20.pbmc <- plot_1D_intervals(sccomp.test.anticd20.pbmc, test_composition_above_logit_fold_change = 0, significance_statistic = "FDR")
sccomp.plot.ntz.pbmc <- plot_1D_intervals(sccomp.test.ntz.pbmc, test_composition_above_logit_fold_change = 0, significance_statistic = "FDR")

write_rds(list(model = sccomp.res.pbmc,
              anticd20 = list(res = sccomp.test.anticd20.pbmc,
                               plot = sccomp.plot.anticd20.pbmc),
               NTZ = list(res = sccomp.test.ntz.pbmc,
                          plot = sccomp.plot.ntz.pbmc)), "results/proportional_changes/sccomp_pbmc_results.rds")

