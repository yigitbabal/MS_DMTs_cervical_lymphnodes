# Read10x, filter gene annotation and create Seurat object
## arguments
### dir.path: path of directory of 10x files. sub-folders should contain 10x files each samples
### mt.pattern: gene patterns for mitochondrial genes (e.g ^MT- for human, ^Mt- for mouse)
### remove ncRNA: logical to remove ncRNA genes from count matrix
### remove snRNA: logical to remove snRNA genes from count matrix
### mincells: Include features detected in at least this many cells. (see ?CreateSeuratObject)
### minfeature: Include cells where at least this many features are detected ((see ?CreateSeuratObject))

counts_to_seurat <- function(dir.path,
                             mt.pattern = "^MT-",
                             remove.ncRNA = FALSE,
                             remove.snRNA = FALSE,
                             mincells = 3,
                             minfeatures = 200){
  
  # 1. Safer directory listing (prevents grabbing the parent dir or non-sample folders)
  samples <- list.dirs(path = dir.path, full.names = TRUE, recursive = FALSE)
  samples.names <- basename(samples)
  
  counts <- lapply(samples, function(sample_path){
    print(paste("Processing sample:", sample_path))
    data <- Read10X(data.dir = sample_path)
    
    # 2. Safely handle multi-modal 10X data (e.g., Gene Expression + Antibody Capture)
    if (inherits(data, "list")) {
      data <- data[[1]] 
    }
    
    # 3. Consolidated logic (removed the duplicate if/else block)
    gene.anot <- data.frame(gene = rownames(data))
    gene.anot$id <- gene.anot$gene
    
    # Check for ENSEMBL IDs
    ensg_idx <- grepl(pattern = "^ENSG", x = gene.anot$gene)
    if(any(ensg_idx)){
      gene.anot$id[ensg_idx] <- mapIds(x = org.Hs.eg.db, 
                                       keys = gene.anot$gene[ensg_idx], 
                                       column = "SYMBOL", 
                                       keytype = "ENSEMBL")
    }
    
    # Map Entrez and Gene types
    gene.anot$entrez <- mapIds(x = org.Hs.eg.db, keys = gene.anot$id, keytype = "SYMBOL", column = "ENTREZID")
    gene.anot$genetype <- mapIds(x = org.Hs.eg.db, keys = gene.anot$id, keytype = "SYMBOL", column = "GENETYPE")
    
    # 4. Use your mt.pattern argument (previously, "^MT-" was hardcoded here)
    mito.genes <- gene.anot[grepl(pattern = mt.pattern, x = gene.anot$id), "id"]
    mito.genes <- mito.genes[!grepl(pattern = paste0(mt.pattern, "T"), x = mito.genes)] 
    
    if(remove.ncRNA){
      gene.anot <- gene.anot[!grepl(pattern = "ncRNA", x = gene.anot$genetype), ]
    }
    
    if(remove.snRNA){
      gene.anot <- gene.anot[!grepl(pattern = "snRNA", x = gene.anot$genetype), ]
    }
    
    gene.anot <- gene.anot[!is.na(gene.anot$id), ]
    gene.anot <- gene.anot[!grepl(pattern = "^ENSG", x = gene.anot$gene),]
    
    data <- data[rownames(data) %in% c(gene.anot$gene, mito.genes), ]
    return(data)
  })
  
  names(counts) <- samples.names
  
  
  
  seurat.list <- lapply(seq_along(counts), function(i){
    # 5. Moved the print statement BEFORE the return statement so it actually prints
    print(paste("Creating Seurat object for sample:", samples.names[i])) 
    
    colnames(counts[[i]]) <- paste0(samples.names[i], "_", colnames(counts[[i]]))
    
    s <- CreateSeuratObject(counts = counts[[i]], 
                            min.cells = mincells, 
                            min.features = minfeatures)
    s$sample <- samples.names[i]
    
    # 6. Use your mt.pattern argument here too!
    s[["percent_mito"]] <- PercentageFeatureSet(s, pattern = mt.pattern)
    return(s)
  })
  
  names(seurat.list) <- samples.names
  
  print("Seurat objects created for all samples.")
  return(seurat.list)
}




# QC plot
## arguments
### object: Seurat object to QC
### sample_col: column name in meta.data that contains sample labels (for grouping in plots)
### mito_pattern: gene pattern for mitochondrial genes (e.g. "^MT-" for human,
###               "^mt-" for mouse). Used to calculate percent_mito if not already present.
### k_mad: number of MADs to use for suggesting thresholds (default 3)
### show_guides: logical to overlay suggested thresholds on plots (default TRUE)
### show_ribo: logical to include percent_ribo plot (default FALSE). If TRUE, percent_ribo is calculated as percentage of ribosomal protein genes (RPS and RPL) if not already present.


preQCSeurat <- function(object,
                        sample_col   = "sample",
                        mito_pattern = "^MT-",      # use "^mt-" for mouse
                        k_mad        = 3,
                        show_guides  = TRUE,
                        show_ribo    = FALSE,
                        saveplot     = FALSE,
                        plot_path    = NULL) {
  
  stopifnot(sample_col %in% colnames(object@meta.data))
  
  seurat <- object
  
  sample_id <- unique(seurat[[sample_col]])
  
  # ---- QC features ----
  seurat[["percent_mito"]] <- PercentageFeatureSet(seurat, pattern = mito_pattern)
  
  small_subunit <- c("RPSA","RPS2","RPS3","RPS3A","RPS4X","RPS4Y","RPS5","RPS6","RPS7","RPS8","RPS9",
                     "RPS10","RPS11","RPS12","RPS13","RPS14","RPS15","RPS15A","RPS16","RPS17","RPS18",
                     "RPS19","RPS20","RPS21","RPS23","RPS24","RPS25","RPS26","RPS27","RPS27A","RPS28",
                     "RPS29","RPS30")
  large_subunit <- c("RPL3","RPL4","RPL5","RPL6","RPL7","RPL7A","RPL8","RPL9","RPL10","RPL10A","RPL11",
                     "RPL12","RPL13","RPL13A","RPL14","RPL15","RPL17","RPL18","RPL18A","RPL19","RPL21",
                     "RPL22","RPL23","RPL23A","RPL24","RPL26","RPL27","RPL27A","RPL28","RPL29","RPL30",
                     "RPL31","RPL32","RPL34","RPL35","RPL35A","RPL36","RPL36A","RPL37","RPL37A","RPL38",
                     "RPL39","RPL40","RPL41","RPLP0","RPLP1","RPLP2","RPLP3")
  ribosomal_protein_genes <- c(small_subunit, large_subunit)
  ribo.genes <- intersect(rownames(seurat), ribosomal_protein_genes)
  seurat[["percent_ribo"]] <- PercentageFeatureSet(seurat, features = ribo.genes)
  
  # ---- helper: suggest thresholds (per sample) ----
  .suggest_qc_thresholds <- function(seurat, sample_col, k_mad) {
    df <- FetchData(
      seurat,
      vars = c("nFeature_RNA","nCount_RNA","percent_mito","percent_ribo", sample_col)
    ) |>
      dplyr::rename(sample = !!sample_col)
    
    df |>
      dplyr::group_by(sample) |>
      dplyr::summarize(
        nFeature_low  = {
          x <- log10(nFeature_RNA + 1); med <- median(x, na.rm = TRUE); madv <- mad(x, na.rm = TRUE)
          pmax(0, round(10^(med - k_mad*madv) - 1))
        },
        nFeature_high = round(quantile(nFeature_RNA, 0.995, na.rm = TRUE)),
        nCount_low    = {
          x <- log10(nCount_RNA + 1); med <- median(x, na.rm = TRUE); madv <- mad(x, na.rm = TRUE)
          pmax(0, round(10^(med - k_mad*madv) - 1))
        },
        nCount_high   = round(quantile(nCount_RNA, 0.995, na.rm = TRUE)),
        mito_high_mad = {
          x <- percent_mito; med <- median(x, na.rm = TRUE); madv <- mad(x, na.rm = TRUE)
          med + k_mad*madv
        },
        mito_high_q95 = as.numeric(quantile(percent_mito, 0.95, na.rm = TRUE)),
        ribo_high_q95 = as.numeric(quantile(percent_ribo, 0.95, na.rm = TRUE)),
        .groups = "drop"
      ) |>
      dplyr::mutate(percent_mito_high = pmin(mito_high_mad, mito_high_q95)) |>
      dplyr::select(
        sample,
        nFeature_low, nFeature_high,
        nCount_low, nCount_high,
        percent_mito_high, ribo_high_q95
      )
  }
  
  # ---- base plots ----
  feats <- c("nFeature_RNA","nCount_RNA","percent_mito")
  if (show_ribo) feats <- c(feats, "percent_ribo")
  
  # list of ggplots rather than combined patchwork
  p_list <- VlnPlot(
    seurat,
    features = feats,
    group.by = sample_col,
    ncol = length(feats),
    raster = FALSE,
    pt.size = 0,
    combine = FALSE
  )
  names(p_list) <- feats
  
  p_scatter <- FeatureScatter(
    seurat,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA",
    group.by = sample_col
  )
  
  # ---- overlay guides ----
  thresholds <- NULL
  if (show_guides) {
    thresholds <- .suggest_qc_thresholds(seurat, sample_col, k_mad)
    
    ## nFeature_RNA: low & high per sample
    if ("nFeature_RNA" %in% names(p_list)) {
      p_list[["nFeature_RNA"]] <- p_list[["nFeature_RNA"]] +
        ggplot2::geom_crossbar(
          data = thresholds,
          ggplot2::aes(x = sample, y = nFeature_low, ymin = nFeature_low, ymax = nFeature_low),
          inherit.aes = FALSE,
          width = 0.7,
          linetype = "dashed"
        ) +
        ggplot2::geom_crossbar(
          data = thresholds,
          ggplot2::aes(x = sample, y = nFeature_high, ymin = nFeature_high, ymax = nFeature_high),
          inherit.aes = FALSE,
          width = 0.7,
          linetype = "dashed"
        )
    }
    
    ## nCount_RNA: low & high per sample
    if ("nCount_RNA" %in% names(p_list)) {
      p_list[["nCount_RNA"]] <- p_list[["nCount_RNA"]] +
        ggplot2::geom_crossbar(
          data = thresholds,
          ggplot2::aes(x = sample, y = nCount_low, ymin = nCount_low, ymax = nCount_low),
          inherit.aes = FALSE,
          width = 0.7,
          linetype = "dashed"
        ) +
        ggplot2::geom_crossbar(
          data = thresholds,
          ggplot2::aes(x = sample, y = nCount_high, ymin = nCount_high, ymax = nCount_high),
          inherit.aes = FALSE,
          width = 0.7,
          linetype = "dashed"
        )
    }
    
    ## mito guide only on percent_mito plot
    if ("percent_mito" %in% names(p_list)) {
      p_list[["percent_mito"]] <- p_list[["percent_mito"]] +
        ggplot2::stat_summary(
          data = thresholds,
          ggplot2::aes(x = sample, y = percent_mito_high),
          fun = "identity",
          geom = "crossbar",
          width = 0.7,
          inherit.aes = FALSE
        )
    }
    
    ## ribo guide only on percent_ribo plot
    if (show_ribo && "percent_ribo" %in% names(p_list)) {
      p_list[["percent_ribo"]] <- p_list[["percent_ribo"]] +
        ggplot2::stat_summary(
          data = thresholds,
          ggplot2::aes(x = sample, y = ribo_high_q95),
          fun = "identity",
          geom = "crossbar",
          width = 0.7,
          inherit.aes = FALSE
        )
    }
    
    # Global min/max guides on scatter
    x_min <- min(thresholds$nCount_low,  na.rm = TRUE)
    x_max <- max(thresholds$nCount_high, na.rm = TRUE)
    y_min <- min(thresholds$nFeature_low,  na.rm = TRUE)
    y_max <- max(thresholds$nFeature_high, na.rm = TRUE)
    
    p_scatter <- p_scatter +
      ggplot2::geom_vline(xintercept = x_min, linetype = "dashed") +
      ggplot2::geom_vline(xintercept = x_max, linetype = "dashed") +
      ggplot2::geom_hline(yintercept = y_min, linetype = "dashed") +
      ggplot2::geom_hline(yintercept = y_max, linetype = "dashed") +
      ggplot2::labs(subtitle = sprintf(
        "Suggested guides: nCount [%s–%s], nFeature [%s–%s]",
        format(x_min, big.mark = ","), format(x_max, big.mark = ","),
        format(y_min, big.mark = ","), format(y_max, big.mark = ",")
      ))
  }
  
  # combine violins into single patchwork layout
  p_violin <- patchwork::wrap_plots(p_list, ncol = length(p_list), guides = "collect")
  
  
  if (show_guides) {
    plot_name = "_preQC.pdf"
  } else {
    plot_name = "_postQC.pdf"
  }
  
  if (saveplot && !is.null(plot_path)) {
    ggsave(paste(plot_path, paste(sample_id, plot_name), sep = "/"), p_violin / p_scatter, width = 300, height = 200, units = "mm")
  }
  
  # ---- output ----
  list(
    plot = p_violin / p_scatter,
    thresholds = thresholds
  )
}




# Basic seurat workflow

process_seurat <- function(obj, 
                                          integrate = FALSE, 
                                          use_jackstraw = TRUE,
                                          res_range = c(0.2, 0.5, 0.8, 1.2),
                                          js_threshold = 0.05,
                                          project_name = "MultiResProject") {
  
  # 1. Initial Merge
  merged_obj <- obj
  
  # 2. Standard Normalization and Scaling
  merged_obj <- NormalizeData(merged_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  merged_obj <- FindVariableFeatures(merged_obj, selection.method = "vst", nfeatures = 2000)
  merged_obj <- ScaleData(merged_obj)
  merged_obj <- RunPCA(merged_obj, features = VariableFeatures(object = merged_obj))
  
  # 3. JackStraw Procedure (Statistically determine PC significance)
  # Note: JackStraw is computationally expensive on large datasets
  dims_to_use <- 1:20 # Default fallback
  
  if (use_jackstraw) {
    message("Running JackStraw (this may take a moment)...")
    merged_obj <- JackStraw(merged_obj, num.replicate = 100, dims = 50)
    merged_obj <- ScoreJackStraw(merged_obj, dims = 1:50)
    
    # Extract p-values for each PC
    js_stats <- merged_obj[["pca"]]@jackstraw@overall.p.values
    
    # Select PCs where p-value < threshold
    significant_pcs <- js_stats[js_stats[, "Score"] < js_threshold, "PC"]
    
    if (length(significant_pcs) > 0) {
      dims_to_use <- significant_pcs
      message(paste("JackStraw selected the first", length(significant_pcs), "PCs."))
    }
  }
  
  # 4. Clustering & UMAP (Integrated vs Standard)
  if (integrate) {
    message("Performing Integration...")
    merged_obj <- IntegrateLayers(
      object = merged_obj, 
      method = RPCAIntegration, 
      orig.reduction = "pca", 
      new.reduction = "integrated.rpca"
    )
    reduction_target <- "integrated.rpca"
  } else {
    reduction_target <- "pca"
  }
  
  # 5. Neighbors and UMAP
  merged_obj <- FindNeighbors(merged_obj, reduction = reduction_target, dims = dims_to_use)
  merged_obj <- RunUMAP(merged_obj, reduction = reduction_target, dims = dims_to_use)
  
  # 6. Multi-Resolution Clustering
  message(paste("Calculating clusters for resolutions:", paste(res_range, collapse=", ")))
  merged_obj <- FindClusters(merged_obj, resolution = res_range)
  
  return(merged_obj)
}




# This function returns annotations form singleR and scCATCH.
suppressPackageStartupMessages(library(SingleR))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(celldex))

annotateCells <- function(srtObj, cluster_col = "seurat_clusters", cluster_level = TRUE) {
  library(SingleR)
  library(SummarizedExperiment)
  library(celldex)
  
  # 1. Ensure references exist (Downloading only if necessary)
  if (!exists("hpca.se")) hpca.se <- HumanPrimaryCellAtlasData()
  if (!exists("immune.monaco")) immune.monaco <- MonacoImmuneData()
  
  # 2. Prepare Data
  # JoinLayers ensures Seurat v5 data is accessible as a single matrix
  srtObj <- JoinLayers(srtObj)
  normData <- GetAssayData(srtObj, assay = "RNA", slot = "data")
  
  # Define clusters if using cluster-level annotation
  clusters_to_use <- NULL
  if (cluster_level) {
    message("Mode: Cluster-level annotation (Memory Efficient)")
    clusters_to_use <- srtObj[[cluster_col, drop = TRUE]]
  } else {
    message("Mode: Cell-level annotation (High Memory Usage)")
  }
  
  # 3. Internal helper to run SingleR and handle mapping
  run_singler_logic <- function(ref_data, label_type) {
    message(paste("Annotating:", label_type))
    
    # Run SingleR
    results <- SingleR(
      test = normData, 
      ref = ref_data, 
      labels = ref_data[[label_type]],
      clusters = clusters_to_use, # NULL if cell-level, Vector if cluster-level
      assay.type.test = 1
    )
    
    # Mapping logic
    if (cluster_level) {
      # Map the cluster-level prediction back to every cell in that cluster
      # results is indexed by cluster name
      final_labels <- results$pruned.labels[match(clusters_to_use, rownames(results))]
    } else {
      # results is indexed by cell name
      final_labels <- results$pruned.labels
    }
    
    final_labels[is.na(final_labels)] <- "unannotated"
    return(final_labels)
  }
  
  # 4. Execute Annotations
  output <- list(
    hcmain  = run_singler_logic(hpca.se, "label.main"),
    hcfine  = run_singler_logic(hpca.se, "label.fine"),
    immMain = run_singler_logic(immune.monaco, "label.main"),
    immFine = run_singler_logic(immune.monaco, "label.fine")
  )
  
  return(output)
}

cluster_states <- function(x, y, z){
  
  
  meta.data.seurat <- x@meta.data
  
  meta.data.seurat %>%
    dplyr::select(!!sym(z)) %>%
    dplyr::group_by(!!sym(z)) %>%
    summarise(number = n()) %>%
    as.data.frame() -> total.number
  
  meta.data.seurat %>%
    dplyr::select(!!sym(z), !!sym(y)) %>%
    dplyr::group_by(!!sym(z), !!sym(y)) %>%
    summarise(number = n()) %>%
    as.data.frame() -> cell.number
  
  cell.number$total <- total.number$number[match(cell.number[,z], total.number[,z])]
  cell.number$percent <- (cell.number$number / cell.number$total) * 100
  
  plot <- cell.number %>%
    ggplot(aes(y = !!sym(z), x = round(percent, digits = 2), fill = !!sym(y))) +
    geom_bar(stat = "identity") + 
    theme_bw() + ylab("") + xlab("percent (%)")
  
  return(plot)
}


calculate_celltype_percentages <- function(seurat_obj, sample_col, celltype_col) {
  
  # 1. Extract metadata
  meta <- seurat_obj@meta.data
  
  # 2. Safety check
  if (!sample_col %in% colnames(meta)) {
    stop(paste("Column", sample_col, "not found in Seurat metadata."))
  }
  if (!celltype_col %in% colnames(meta)) {
    stop(paste("Column", celltype_col, "not found in Seurat metadata."))
  }
  
  # Create symbols for tidy evaluation
  s_col <- sym(sample_col)
  c_col <- sym(celltype_col)
  
  # 3. Calculate with explicit package routing
  percentage_df <- meta %>%
    # Group and count (avoids the count() masking issue)
    dplyr::group_by(!!s_col, !!c_col) %>%
    dplyr::summarise(Cell_Count = n(), .groups = 'drop') %>%
    # Force all missing combinations to exist with 0 counts
    tidyr::complete(!!s_col, !!c_col, fill = list(Cell_Count = 0)) %>%
    # Regroup by sample to calculate the total cells and proportions
    dplyr::group_by(!!s_col) %>%
    dplyr::mutate(
      Total_Cells_in_Sample = sum(Cell_Count),
      Percentage = (Cell_Count / Total_Cells_in_Sample) * 100
    ) %>%
    dplyr::ungroup()
  
  return(percentage_df)
}


### DEG analysis
# Load required libraries
library(Seurat)
library(DESeq2)
library(edgeR)
library(limma)
library(dplyr)
library(tibble)
library(Matrix)

RunConsensusPseudobulk <- function(seurat_obj, 
                                   cell_type_col = "cell_type",
                                   treatment_col = "treatment",
                                   sample_col = "sample_id",
                                   baseline_group = "before",
                                   after_groups = c("after_1", "after_2", "after_3"),
                                   p_val_cutoff = 0.05 
) {
  
  message("Aggregating single-cell counts to pseudobulk...")
  
  # ---------------------------------------------------------
  # 1. SETUP
  # ---------------------------------------------------------
  seurat_obj$safe_ct <- gsub("_| ", "-", seurat_obj[[cell_type_col]][[1]])
  seurat_obj$safe_sample <- gsub("_| ", "-", seurat_obj[[sample_col]][[1]])
  
  sample_meta <- seurat_obj@meta.data %>%
    dplyr::select(safe_sample, all_of(treatment_col)) %>%
    distinct() %>%
    remove_rownames() %>%
    column_to_rownames(var = "safe_sample")
  
  agg_counts <- AggregateExpression(seurat_obj, 
                                    group.by = c("safe_ct", "safe_sample"),
                                    assays = "RNA", 
                                    return.seurat = FALSE)$RNA
  
  safe_cell_types <- unique(seurat_obj$safe_ct)
  
  all_results_list <- list()
  consensus_sig_list <- list()
  
  stats_formula <- as.formula(paste("~", treatment_col))
  
  # ---------------------------------------------------------
  # 2. ANALYSIS LOOP
  # ---------------------------------------------------------
  for (ct_safe in safe_cell_types) {
    
    orig_ct <- unique(seurat_obj[[cell_type_col]][[1]][seurat_obj$safe_ct == ct_safe])[1]
    
    ct_cols <- grep(paste0("^", ct_safe, "_"), colnames(agg_counts), value = TRUE)
    
    if (length(ct_cols) < 2) next 
    
    ct_counts_full <- agg_counts[, ct_cols, drop = FALSE]
    colnames(ct_counts_full) <- gsub(paste0("^", ct_safe, "_"), "", colnames(ct_counts_full))
    
    valid_samples <- intersect(colnames(ct_counts_full), rownames(sample_meta))
    ct_counts_full <- ct_counts_full[, valid_samples, drop = FALSE]
    ct_meta_full <- sample_meta[valid_samples, , drop = FALSE]
    
    for (trt in after_groups) {
      
      if (trt %in% ct_meta_full[[treatment_col]] && baseline_group %in% ct_meta_full[[treatment_col]]) {
        
        message(paste("Running consensus DGE for:", orig_ct, "|", trt, "vs", baseline_group))
        
        # --- Subset Pseudobulk Data ---
        samples_to_keep <- rownames(ct_meta_full)[ct_meta_full[[treatment_col]] %in% c(trt, baseline_group)]
        sub_counts <- ct_counts_full[, samples_to_keep, drop = FALSE]
        
        # We still remove genes with absolute 0 counts across all samples in this subset
        # to prevent the statistical models from crashing
        sub_counts <- sub_counts[rowSums(sub_counts) > 0, , drop = FALSE] 
        
        if(nrow(sub_counts) < 2 || ncol(sub_counts) < 2) {
          message("  -> Skipping: Not enough genes or samples left after filtering.")
          next
        }
        
        sub_meta <- ct_meta_full[samples_to_keep, , drop = FALSE]
        sub_meta[[treatment_col]] <- factor(sub_meta[[treatment_col]], levels = c(baseline_group, trt))
        
        # --- A. DESeq2 ---
        dds <- DESeqDataSetFromMatrix(countData = sub_counts, colData = sub_meta, design = stats_formula)
        dds <- DESeq(dds, quiet = TRUE)
        res_deseq <- results(dds, name = paste0(treatment_col, "_", trt, "_vs_", baseline_group)) %>%
          as.data.frame() %>% rownames_to_column("gene") %>%
          dplyr::select(gene, logFC_DESeq2 = log2FoldChange, padj_DESeq2 = padj)
        
        # --- B. edgeR ---
        y_edger <- DGEList(counts = sub_counts, group = sub_meta[[treatment_col]])
        y_edger <- calcNormFactors(y_edger)
        design <- model.matrix(stats_formula, data = sub_meta)
        y_edger <- estimateDisp(y_edger, design)
        fit_edger <- glmQLFit(y_edger, design)
        qlf <- glmQLFTest(fit_edger, coef = 2)
        res_edger <- topTags(qlf, n = Inf)$table %>%
          rownames_to_column("gene") %>%
          dplyr::select(gene, logFC_edgeR = logFC, padj_edgeR = FDR)
        
        # --- C. limma-voom ---
        y_limma <- DGEList(counts = sub_counts)
        y_limma <- calcNormFactors(y_limma)
        v <- voom(y_limma, design, plot = FALSE) 
        fit_limma <- lmFit(v, design)
        fit_limma <- eBayes(fit_limma)
        res_limma <- topTable(fit_limma, coef = 2, number = Inf) %>%
          rownames_to_column("gene") %>%
          dplyr::select(gene, logFC_limma = logFC, padj_limma = adj.P.Val)
        
        # --- Merge and Score ---
        merged_res <- res_deseq %>%
          full_join(res_edger, by = "gene") %>%
          full_join(res_limma, by = "gene") %>%
          mutate(
            cell_type = orig_ct, 
            comparison = paste0(trt, "_vs_", baseline_group),
            sig_DESeq2 = ifelse(!is.na(padj_DESeq2) & padj_DESeq2 < p_val_cutoff, 1, 0),
            sig_edgeR = ifelse(!is.na(padj_edgeR) & padj_edgeR < p_val_cutoff, 1, 0),
            sig_limma = ifelse(!is.na(padj_limma) & padj_limma < p_val_cutoff, 1, 0),
            consensus_score = sig_DESeq2 + sig_edgeR + sig_limma
          )
        
        list_name <- paste(orig_ct, trt, "vs", baseline_group, sep = "_")
        all_results_list[[list_name]] <- merged_res
        
        consensus_sig_list[[list_name]] <- merged_res %>% filter(consensus_score >= 2)
      }
    }
  }
  
  # ---------------------------------------------------------
  # 3. RETURN OUTPUTS
  # ---------------------------------------------------------
  final_consensus_sig_df <- bind_rows(consensus_sig_list)
  
  message("Done!")
  
  return(list(
    all_results = all_results_list,
    significant_consensus = final_consensus_sig_df
  ))
}

# Load required libraries
library(Seurat)
library(dplyr)
library(tibble)

RunSeuratWilcoxon <- function(seurat_obj, 
                              cell_type_col = "cell_type",
                              treatment_col = "treatment",
                              baseline_group = "before",
                              after_groups = c("after_1", "after_2", "after_3"),
                              min.pct = 0,          # Set to 0 to compare all genes against pseudobulk
                              logfc.threshold = 0,  # Set to 0 to compare all genes against pseudobulk
                              p_val_cutoff = 0.05 
) {
  
  message("Running standard Seurat single-cell Wilcoxon tests...")
  
  # ---------------------------------------------------------
  # 1. SETUP: Create combined identity class
  # ---------------------------------------------------------
  # Create a temporary column combining cell type and treatment (e.g., "Th17 Cell_antiCD20")
  seurat_obj$celltype_treatment <- paste(seurat_obj[[cell_type_col]][[1]], 
                                         seurat_obj[[treatment_col]][[1]], 
                                         sep = "_")
  
  # Set this new column as the active identity for FindMarkers
  Idents(seurat_obj) <- "celltype_treatment"
  
  cell_types <- unique(seurat_obj[[cell_type_col]][[1]])
  
  all_results_list <- list()
  sig_results_list <- list()
  
  # ---------------------------------------------------------
  # 2. ANALYSIS LOOP
  # ---------------------------------------------------------
  for (ct in cell_types) {
    for (trt in after_groups) {
      
      # Construct the exact group names we are looking for
      ident_1 <- paste(ct, trt, sep = "_")            # e.g., "Th17 Cell_antiCD20" (Numerator/After)
      ident_2 <- paste(ct, baseline_group, sep = "_") # e.g., "Th17 Cell_BeforeTreatment" (Denominator/Before)
      
      # Check if both groups actually exist in the data
      if (ident_1 %in% Idents(seurat_obj) && ident_2 %in% Idents(seurat_obj)) {
        
        # Check if there are enough cells to run the test (Seurat prefers >3)
        cells_1 <- sum(Idents(seurat_obj) == ident_1)
        cells_2 <- sum(Idents(seurat_obj) == ident_2)
        
        if (cells_1 > 3 && cells_2 > 3) {
          
          message(paste("Running Wilcoxon for:", ct, "|", trt, "vs", baseline_group))
          
          # Run standard Seurat DGE
          res <- FindMarkers(seurat_obj, 
                             ident.1 = ident_1, 
                             ident.2 = ident_2,
                             test.use = "wilcox",
                             min.pct = min.pct,
                             logfc.threshold = logfc.threshold,
                             verbose = FALSE)
          
          # Tidy up results to match our previous format
          res <- res %>%
            rownames_to_column("gene") %>%
            mutate(
              cell_type = ct,
              comparison = paste0(trt, "_vs_", baseline_group),
              # Flag significance based on Seurat's adjusted p-value (Bonferroni)
              sig_wilcox = ifelse(!is.na(p_val_adj) & p_val_adj < p_val_cutoff, 1, 0)
            )
          
          list_name <- paste(ct, trt, "vs", baseline_group, sep = "_")
          all_results_list[[list_name]] <- res
          
          # Filter for significant genes
          sig_results_list[[list_name]] <- res %>% filter(sig_wilcox == 1)
          
        } else {
          message(paste("  -> Skipping:", ct, trt, "vs", baseline_group, "- Not enough cells."))
        }
      } 
    }
  }
  
  # ---------------------------------------------------------
  # 3. RETURN OUTPUTS
  # ---------------------------------------------------------
  final_sig_df <- bind_rows(sig_results_list)
  
  message("Done!")
  
  return(list(
    all_results = all_results_list,
    significant_results = final_sig_df
  ))
}



## CellChat
library(Seurat)
library(CellChat)

run_comparative_cellchat_lapply <- function(seurat_obj, celltype_col, condition_col, min_cells = 10) {
  
  # 2. Split the Seurat object by the condition (e.g., Pre/Post treatment)
  message("Splitting Seurat object by ", condition_col, "...")
  seurat.list <- SplitObject(seurat_obj, split.by = condition_col)
  
  # 3. Create CellChat objects directly from the Seurat objects
  message("Creating CellChat objects...")
  cellChat.list <- lapply(seurat.list, function(x) {
    createCellChat(object = x, group.by = celltype_col, assay = "RNA")
  })
  
  # 4. Clean up memory (very useful for large MS datasets!)
  rm(seurat.list)
  gc()
  
  # 5. Set Database
  CellChatDB <- CellChatDB.human
  
  # 6. Run the core CellChat pipeline using lapply
  message("Running CellChat pipeline on all conditions (this may take a while)...")
  
  # We can use lapply here. It will automatically preserve the list names (e.g., your condition names)
  cellChat.list <- lapply(cellChat.list, function(cc) {
    
    cc@DB <- CellChatDB
    cc <- subsetData(cc) # This step is necessary even if using the whole database
    cc <- identifyOverExpressedGenes(cc)
    cc <- identifyOverExpressedInteractions(cc)
    
    # Compute probabilities
    cc <- computeCommunProb(cc, type = "triMean")
    cc <- filterCommunication(cc, min.cells = min_cells)
    
    # Pathways and Network Aggregation
    cc <- computeCommunProbPathway(cc)
    cc <- aggregateNet(cc)
    
    # Compute centrality for downstream network analysis
    cc <- netAnalysis_computeCentrality(cc)
    
    return(cc)
  })
  
  message("--- Pipeline Complete! ---")
  return(cellChat.list)
}




# cellchat split bubble plot

library(Seurat)
library(tidyverse)
library(ggnewscale)
library(scales)

plot_split_LR_bubble <- function(LR.data, seurat_obj, celltype_col, condition_col) {
  
  # 1. Create a combined ID column (e.g., "Memory B cells_Pre_Tx")
  seurat_obj$combo_id <- paste(seurat_obj[[celltype_col, drop=TRUE]], 
                               seurat_obj[[condition_col, drop=TRUE]], sep = "_")
  
  # 2. Extract unique genes and split complexes dynamically
  all_features <- unique(c(LR.data$ligand, LR.data$receptor))
  gene_list <- strsplit(all_features, "_")
  names(gene_list) <- all_features
  all_individual_genes <- unique(unlist(gene_list))
  
  # Check against Seurat object to prevent missing-gene errors
  all_genes <- intersect(all_individual_genes, rownames(seurat_obj))
  if(length(all_genes) == 0) stop("None of the genes were found in the Seurat object!")
  
  # 3. Get Expression Data using AverageExpression
  message("Extracting expression data...")
  exp_list <- AverageExpression(seurat_obj, features = all_genes, group.by = "combo_id")
  # Extract the matrix natively (prevents spaces from turning into dots)
  exp_mat <- if("RNA" %in% names(exp_list)) exp_list$RNA else exp_list[[1]]
  
  exp_df_list <- lapply(names(gene_list), function(feat) {
    genes <- intersect(gene_list[[feat]], rownames(exp_mat))
    
    # Air-tight fallback
    if (length(genes) == 0) {
      return(data.frame(feature = feat, combo_id = colnames(exp_mat), avg_exp = 0, stringsAsFactors = FALSE))
    }
    
    sub_exp <- exp_mat[genes, , drop = FALSE]
    avg_vals <- colMeans(sub_exp, na.rm = TRUE)
    
    combo_names <- if(!is.null(names(avg_vals))) names(avg_vals) else colnames(exp_mat)
    
    data.frame(
      feature = feat, 
      combo_id = combo_names, 
      avg_exp = as.numeric(avg_vals), 
      stringsAsFactors = FALSE
    )
  })
  exp_stats <- dplyr::bind_rows(exp_df_list)
  
  # 4. Get Percentage Data using DotPlot
  message("Extracting percentage data...")
  dp <- DotPlot(seurat_obj, features = all_genes, group.by = "combo_id")
  dp_data <- dp$data
  
  pct_df_list <- lapply(names(gene_list), function(feat) {
    genes <- intersect(gene_list[[feat]], rownames(seurat_obj))
    
    if(length(genes) == 0) {
      return(data.frame(feature = feat, combo_id = as.character(unique(dp_data$id)), pct_exp = 0, stringsAsFactors = FALSE))
    }
    
    sub_dp <- dp_data[dp_data$features.plot %in% genes, ]
    if(nrow(sub_dp) == 0) {
      return(data.frame(feature = feat, combo_id = as.character(unique(dp_data$id)), pct_exp = 0, stringsAsFactors = FALSE))
    }
    
    agg <- aggregate(pct.exp ~ id, data = sub_dp, FUN = mean, na.rm = TRUE)
    
    data.frame(
      feature = feat,
      combo_id = as.character(agg$id),
      pct_exp = as.numeric(agg$pct.exp),
      stringsAsFactors = FALSE
    )
  })
  pct_stats <- dplyr::bind_rows(pct_df_list)
  
  # 5. Merge Expression and Percentages
  message("Joining datasets...")
  feature_stats <- dplyr::full_join(exp_stats, pct_stats, by = c("feature", "combo_id"))
  
  # 6. Prepare and Join exactly mapped names to LR.data
  # EXPLICITLY forcing dplyr::rename to completely block plyr from crashing the script
  ligand_mapped <- feature_stats %>% 
    dplyr::rename(ligand = feature, source_combo = combo_id, ligand_exp = avg_exp, ligand_perc = pct_exp)
  
  receptor_mapped <- feature_stats %>% 
    dplyr::rename(receptor = feature, target_combo = combo_id, receptor_exp = avg_exp, receptor_perc = pct_exp)
  
  LR.data <- LR.data %>%
    dplyr::mutate(
      source_combo = paste(source, dataset, sep = "_"),
      target_combo = paste(target, dataset, sep = "_")
    ) %>%
    dplyr::left_join(ligand_mapped, by = c("ligand", "source_combo")) %>%
    dplyr::left_join(receptor_mapped, by = c("receptor", "target_combo"))
  
  # 7. Geometry Prep
  message("Preparing geometric plot...")
  LR.data <- LR.data %>%
    dplyr::mutate(
      # Replace any NAs from failed matches with 0 before scaling
      ligand_exp = tidyr::replace_na(ligand_exp, 0),
      receptor_exp = tidyr::replace_na(receptor_exp, 0),
      ligand_perc = tidyr::replace_na(ligand_perc, 0),
      receptor_perc = tidyr::replace_na(receptor_perc, 0),
      
      ligand_radius = sqrt(ligand_perc / pi) / 10,
      receptor_radius = sqrt(receptor_perc / pi) / 10,
      ligand_exp_scaled = scales::rescale(ligand_exp, to = c(0, 1)),
      receptor_exp_scaled = scales::rescale(receptor_exp, to = c(0, 1)),
      row_id = dplyr::row_number()
    )
  
  x_levels <- unique(LR.data$source)
  y_levels <- unique(LR.data$interaction_name_2)
  
  make_semicircle <- function(center, radius, direction = "left", n = 50) { 
    if (direction == "left") {
      theta <- seq(pi/2, 3*pi/2, length.out = n)
    } else {
      theta <- seq(3*pi/2, pi/2 + 2*pi, length.out = n)
    }
    tibble::tibble(x = center[1] + radius * cos(theta), y = center[2] + radius * sin(theta))
  }
  
  left_polys <- LR.data %>%
    dplyr::rowwise() %>%
    dplyr::mutate(poly = list(
      make_semicircle(c(match(source, x_levels), match(interaction_name_2, y_levels)), ligand_radius, "left") %>% 
        dplyr::mutate(fill = ligand_exp_scaled, semi = "left", row_id = row_id, target = target, dataset = dataset)
    )) %>% dplyr::ungroup() %>% dplyr::select(poly) %>% tidyr::unnest(poly)
  
  right_polys <- LR.data %>%
    dplyr::rowwise() %>%
    dplyr::mutate(poly = list(
      make_semicircle(c(match(source, x_levels), match(interaction_name_2, y_levels)), receptor_radius, "right") %>% 
        dplyr::mutate(fill = receptor_exp_scaled, semi = "right", row_id = row_id, target = target, dataset = dataset)
    )) %>% dplyr::ungroup() %>% dplyr::select(poly) %>% tidyr::unnest(poly)
  
  # 8. Generate the Plot
  message("Generating final ggplot...")
  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = left_polys, ggplot2::aes(x = x, y = y, group = interaction(row_id, semi), fill = fill), color = "black", linewidth = 0.15) +
    ggplot2::scale_fill_gradient(name = "Ligand (scaled exp)", low = "#fee0d2", high = "darkred", limits = c(0, 1)) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_polygon(data = right_polys, ggplot2::aes(x = x, y = y, group = interaction(row_id, semi), fill = fill), color = "black", linewidth = 0.15) +
    ggplot2::scale_fill_gradient(name = "Receptor (scaled exp)", low = "#deebf7", high = "darkblue", limits = c(0, 1)) +
    ggplot2::scale_x_continuous(breaks = seq_along(x_levels), labels = x_levels) +
    ggplot2::scale_y_continuous(breaks = seq_along(y_levels), labels = y_levels) +
    ggplot2::facet_wrap(target~dataset, ncol = 5, scales = "free") +
    ggplot2::labs(x = "Source", y = "Interaction") +
    ggplot2::theme_minimal(base_size = 13) + ggplot2::xlab("") + ggplot2::ylab("") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  return(p)
}


### persistent clones
find_persistent_clones <- function(tcr_list, pre_sample_name, post_sample_name, cloneCall = "CTstrict") {
  
  # 1. Extract the dataframes for the paired samples
  pre_data <- tcr_list[[pre_sample_name]]
  post_data <- tcr_list[[post_sample_name]]
  
  # Safety check: Ensure both samples exist in the list and have cells
  if (is.null(pre_data) || is.null(post_data) || nrow(pre_data) == 0 || nrow(post_data) == 0) {
    message(sprintf("Skipping %s vs %s: Missing data or 0 cells.", pre_sample_name, post_sample_name))
    return(NULL)
  }
  
  # 2. Identify unique clonotypes
  pre_clones <- unique(pre_data[[cloneCall]])
  post_clones <- unique(post_data[[cloneCall]])
  
  # 3. Find the intersection (shared clones) and remove NAs
  shared_clones <- intersect(pre_clones, post_clones)
  shared_clones <- shared_clones[!is.na(shared_clones)]
  
  # 4. Extract the specific cell BARCODES that belong to these shared clones
  pre_barcodes <- pre_data$barcode[pre_data[[cloneCall]] %in% shared_clones]
  post_barcodes <- post_data$barcode[post_data[[cloneCall]] %in% shared_clones]
  all_persistent_barcodes <- c(pre_barcodes, post_barcodes)
  
  # 5. Calculate stats for reporting
  pre_total <- nrow(pre_data)
  post_total <- nrow(post_data)
  pre_shared_count <- length(pre_barcodes)
  post_shared_count <- length(post_barcodes)
  
  # Print a quick summary to the console
  cat(sprintf("\n--- Results for %s vs %s ---\n", pre_sample_name, post_sample_name))
  cat(sprintf("Shared distinct clonotypes: %d\n", length(shared_clones)))
  cat(sprintf("Pre-treatment: %d / %d cells (%.2f%%) are persistent.\n", 
              pre_shared_count, pre_total, (pre_shared_count / pre_total) * 100))
  cat(sprintf("Post-treatment: %d / %d cells (%.2f%%) are persistent.\n", 
              post_shared_count, post_total, (post_shared_count / post_total) * 100))
  
  # 6. Return all useful data as a list
  return(list(
    num_shared_clones = length(shared_clones),
    shared_clonotypes = shared_clones,
    persistent_barcodes = all_persistent_barcodes
  ))
}
