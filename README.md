# MS_DMTs_cervical_lymphnodes

The following summarizes the steps used to perform the analysis and produce the key results in the manuscript "Cervical lymph node B cell subset resists depletion and sustains lytic Epstein Barr-virus signature in patients with multiple sclerosis" by Babal and Sarkkinen et al., (submitted).

The study profiles fine-needle aspirates (FNA) of deep cervical lymph nodes, together with matched peripheral blood mononuclear cells (PBMC), from people with MS and controls. It uses single-cell RNA sequencing (scRNA-seq) with paired single-cell T cell receptor (scTCR-seq) and B cell receptor (scBCR-seq) sequencing.

The scripts below show the analysis workflow and reproduce the main results.
## Repository structure

```
scripts/
├── 00_functions.R                  # shared helper functions, sourced by the other scripts
├── scRNAseq/                       # preprocessing, annotation and downstream scRNA-seq analyses
├── scTCRseq/                       # TCR repertoire and epitope specificity analyses
└── scBCRseq/                       # BCR repertoire and clonality analyses
    └── germline_tree_construction/ # Immcantation-based clonal lineage (germline tree) pipeline
```

## Workflow

### 0. Helper functions
**`scripts/00_functions.R`** is sourced by most scripts. It defines:
- `counts_to_seurat()`: reads the 10x count matrices, maps gene annotation (org.Hs.eg.db), removes ncRNA/snRNA genes and creates per-sample Seurat objects
- `preQCSeurat()`: per-sample, MAD-based QC thresholds and QC plots
- `process_seurat()`: normalisation, selection of significant PCs with JackStraw, RPCA integration, multi-resolution clustering and UMAP
- `annotateCells()`: automatic annotation with SingleR (Human Primary Cell Atlas and Monaco Immune references)
- `calculate_celltype_percentages()`, `cluster_states()`: cell type composition helpers
- `RunConsensusPseudobulk()`: pseudobulk differential expression that combines DESeq2, edgeR and limma
- `RunSeuratWilcoxon()`: single-cell Wilcoxon differential expression
- `run_comparative_cellchat_lapply()`, `plot_split_LR_bubble()`: CellChat analysis by condition and visualisation
- `find_persistent_clones()`: finds TCR clones that are shared between paired pre-treatment and follow-up samples

### 1. scRNA-seq preprocessing, integration and annotation (`scripts/scRNAseq/`)
| Script | Description |
|---|---|
| `scRNAseq_MS_dcLN_seurat.R` | dcLN (FNA) samples: loads the data, runs per-sample QC filtering, merges the samples, runs RPCA integration and clustering, annotates with SingleR (including label transfer from our previous dcLN dataset), uses clustree for cluster stability, curates the annotation with marker genes, and subclusters where needed |
| `scRNAseq_MS_PBMC_seurat.R` | The same workflow for the matched PBMC samples |
| `scRNAseq_MS_dcLN_Bcell_subset.R` | Subsets, re-integrates and subclusters dcLN B cells; annotates B cell subsets (naive, USM, SM, DN, Tbet+CD11c+, IFN-stimulated B, GC dark and light zone, plasmablasts, pre-B) |
| `scRNAseq_MS_PBMC_Bcell_subset.R` | B cell subsetting and subclustering in PBMC |

### 2. Downstream scRNA-seq analyses (`scripts/scRNAseq/`)
| Script | Description |
|---|---|
| `scRNAseq_MS_differential_gene_expression.R` | Adds the B cell subset labels to the main objects (saved as the `*_curated_with_subset.rds` objects that the later scripts use). Runs differential expression between each DMT and before treatment with consensus pseudobulk (DESeq2, edgeR and limma) and single-cell Wilcoxon tests, then GO enrichment (clusterProfiler) and a comparison of DEGs between dcLN and PBMC |
| `scRNAseq_MS_proprotional_changes_sccomp.R` | Tests differences in cell type composition with sccomp (MS vs. controls, and each DMT vs. before treatment) in dcLN and PBMC |
| `scRNAseq_MS_proprotional_changes_boxplots.R` | Cell type proportions for each sample with Wilcoxon tests, B cell subset proportions, and fold and relative percentage changes after treatment |
| `scRNAseq_exhausted_cd8_analysis.R` | Subclusters CD8 T cells, identifies and characterises exhausted / tissue-resident memory CD8 subsets and their proportions, and runs pseudotime trajectory analysis (monocle3) |
| `scRNAseq_MS_cell_to_cell_communication.R` | Compares cell-cell communication between conditions with CellChat, focusing on T cell → B cell ligand-receptor interactions in dcLN and PBMC, including exhausted CD8 T cells |

### 3. scTCR-seq analysis (`scripts/scTCRseq/`)
| Script | Description |
|---|---|
| `scTCRseq_MS_TCR_analysis.R` | Merges Cell Ranger TCR contigs with the Seurat annotations (scRepertoire), then analyses clonality, diversity and clonal expansion for each cell type, extracts expanded CD8 clones for specificity prediction, and finds clones that persist between the pre-treatment and follow-up samples of the same patient |
| `scTCRseq_tcrdist3_epitope_targeting.R` | Predicts the specificity of expanded TCR clones: scans paired α/β TCRs against VDJdb with tcrdist3 (run through `reticulate`) |
| `scTCRseq_tcrdist3_epitope_targeting_barblot.R` | Summarises and plots the predicted epitope / antigen-species hits (e.g. EBV) for each sample and group |

### 4. scBCR-seq analysis (`scripts/scBCRseq/`)
| Script | Description |
|---|---|
| `scBCRseq_MS_dcLN_BCR_analysis.R` | Merges Cell Ranger BCR contigs with the B cell subset annotations (scRepertoire), analyses clonality and diversity for each B cell subset, and exports annotated B cell barcodes (`seurat_barcodes_annotated.csv`) for the germline pipeline |
| `germline_tree_construction/` | Clonal lineage pipeline (Immcantation: Change-O, SHazaM, dowser, IgPhyML), run step by step. Each `.sh` wrapper runs the matching step in an Immcantation Docker container |
| `scBCRseq_B_cell_clonality_analysis.R` | Run after the germline pipeline. Maps the Immcantation clone IDs back to the dcLN B cell Seurat object and analyses clonal sharing between B cell subsets and changes after treatment |

**Germline tree construction steps** (`scripts/scBCRseq/germline_tree_construction/`):
1. `01_merge_samples.sh`: merges the per-sample AIRR rearrangement files and adds sample prefixes and metadata
2. `01b_filter_by_seurat_barcodes.R/.sh`: keeps only the barcodes of QC-passed, annotated B cells from the Seurat object
3. `02_filter_heavy_chain.R/.sh`: keeps productive heavy chains and runs QC
4. `03_create_germlines.sh`: IgBLAST V(D)J assignment (`AssignGenes.py`), `MakeDb.py` and `CreateGermlines.py`
5. `04_distance_threshold.R/.sh`: distance-to-nearest analysis to choose the clonal clustering threshold
6. `05_define_clones.sh`: clonal clustering across samples (`DefineClones.py`)
7. `05b_reclone_germlines.sh`: rebuilds germlines for each clone (`CreateGermlines.py --cloned`)
8. `06_build_trees.R/.sh`: builds lineage trees with dowser (IgPhyML for large clones, maximum-likelihood `pml` for small clones)
9. `07_annotate_and_visualize.R/.sh`: annotates the trees with B cell subset labels and plots them

The shell wrappers read the `PROJECT_DIR`, `DOCKER_IMAGE`, `NPROC` and `CLONE_DIST` settings from a `project_config.sh` file.

## Main software
R (Seurat v5, SingleR/celldex, clustree, scCustomize, SCP, DESeq2, edgeR, limma, clusterProfiler, sccomp, CellChat, monocle3, scRepertoire, tidyverse), Python (tcrdist3), and the Immcantation suite (IgBLAST, Change-O, SHazaM, alakazam, dowser, IgPhyML). See the `library()` calls in each script for the full list.

## Contact
For more details or questions about the code, please contact the corresponding authors.
