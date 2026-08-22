# ============================================================
# HMPV scRNA-seq preprocessing pipeline
#
# Description:
#   This script performs:
#   1. Raw 10X matrix loading
#   2. Seurat object creation
#   3. Multi-sample merging
#   4. Preliminary batch-effect evaluation
#   5. Cell-level quality control
#   6. Doublet detection using scDblFinder
#   7. Normalization, PCA, graph-based clustering, and UMAP
#
# Recommended software environment:
#   Seurat >= 5.0
#   SeuratObject >= 5.0
#   SingleCellExperiment
#   scDblFinder
#   clustree
#   kBET
#
# Important:
#   Raw count matrices and RDS files should not be uploaded to GitHub.
#   Use .gitignore to exclude large data files and generated outputs.
# ============================================================


# ============================================================
# 1. Environment setup
# ============================================================

.libPaths("/home/lzh/miniconda3/envs/seurat5/lib/R/library")

set.seed(1234)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(dplyr)
  library(ggplot2)
  library(clustree)
})

message("Seurat version: ", packageVersion("Seurat"))
message("SeuratObject version: ", packageVersion("SeuratObject"))
message("SingleCellExperiment version: ",
        packageVersion("SingleCellExperiment"))
message("scDblFinder version: ", packageVersion("scDblFinder"))


# ============================================================
# 2. Project paths
# ============================================================

project_dir <- "/mnt/sda1/lzh/test/cjl/hmpv"

pre_dir <- file.path(project_dir, "pre")
pre_figure_dir <- file.path(pre_dir, "figures")
pre_table_dir <- file.path(pre_dir, "tables")
pre_log_dir <- file.path(pre_dir, "logs")

required_directories <- c(
  pre_dir,
  pre_figure_dir,
  pre_table_dir,
  pre_log_dir
)

invisible(
  lapply(
    required_directories,
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  )
)


# ============================================================
# 3. Utility functions
# ============================================================

save_pdf <- function(plot_object, file, width, height) {
  ggsave(
    filename = file,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    device = "pdf",
    limitsize = FALSE
  )
}


check_directory <- function(path, sample_id) {
  if (!dir.exists(path)) {
    stop(
      "Input directory does not exist for sample ",
      sample_id,
      ": ",
      path
    )
  }

  message("Input directory verified for sample ", sample_id)
}


read_10x_counts <- function(path, sample_id) {
  check_directory(path, sample_id)

  message("Reading 10X matrix for sample: ", sample_id)
  message("Directory: ", path)

  counts <- Read10X(data.dir = path)

  # Read10X() may return either:
  # 1. A sparse count matrix
  # 2. A list containing multiple modalities
  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      stop(
        "The 10X directory for sample ",
        sample_id,
        " contains multiple modalities, but ",
        "'Gene Expression' was not found."
      )
    }
  }

  if (!inherits(counts, "Matrix")) {
    stop(
      "The count matrix for sample ",
      sample_id,
      " is not a sparse Matrix object."
    )
  }

  if (ncol(counts) == 0 || nrow(counts) == 0) {
    stop(
      "The count matrix for sample ",
      sample_id,
      " is empty."
    )
  }

  seurat_object <- CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.cells = 0,
    min.features = 0
  )

  # Keep the original sample identity explicitly.
  seurat_object$sample_id <- sample_id

  return(seurat_object)
}


write_object_summary <- function(object, file) {
  summary_lines <- c(
    capture.output(print(object)),
    "",
    "Assays:",
    capture.output(print(Assays(object))),
    "",
    "Reductions:",
    capture.output(print(Reductions(object))),
    "",
    "Metadata columns:",
    capture.output(print(colnames(object@meta.data))),
    "",
    "Cell counts by sample:",
    capture.output(print(table(object$sample_id)))
  )

  writeLines(summary_lines, con = file)
}


# ============================================================
# 4. Define the input samples
# ============================================================

sample_manifest <- data.frame(
  sample_id = c(
    "hMPVsevere01",
    "hMPVsevere02",
    "hMPVsevere03",
    "hMPVsevere04",
    "hMPVsevere05",
    "Control02",
    "Control03",
    "Control04",
    "hMPVmild01",
    "hMPVmild02",
    "hMPVmild03",
    "hMPVmild04",
    "RSVmild01",
    "RSVsevere01",
    "RSVsevere02"
  ),
  input_dir = c(
    "/mnt/sda1/lzh/data/inhouse/HMPV_RSV/HMPV-ZZ-01/20250429/Matrix/HMPV-ZZ-01_EmptyDrops_CR_matrix",
    "/mnt/sda1/lzh/data/inhouse/HMPV_RSV/HMPV-ZZ-02/20250429/Matrix/HMPV-ZZ-02_EmptyDrops_CR_matrix",
    "/mnt/sda1/lzh/data/inhouse/HMPV_RSV/HMPV-ZZ-03",
    "/mnt/sda1/lzh/data/inhouse/HMPV_RSV/HMPV-ZZ-04",
    "/mnt/sda1/lzh/data/inhouse/HMPV_RSV/HMPV-ZZ-05",
    "/mnt/sda1/lzh/data/inhouse/HMPV_RSV/DZ-02/20250429/Matrix",
    "/mnt/sda1/lzh/data/inhouse/HMPV_RSV/DZ-03/20250429/Matrix",
    "/mnt/sda1/lzh/data/inhouse/HMPV_RSV/DZ-04/20250429/Matrix",
    "/mnt/sda1/lzh/data/inhouse/HMPV_RSV/HMPV-QZ-01/20250429/Matrix",
    "/mnt/sda1/lzh/data/inhouse/HMPV_RSV/HMPV-QZ-02/20250429/Matrix",
    "/mnt/sda1/lzh/data/inhouse/HMPV_RSV/HMPV-QZ-03/20250429/Matrix",
    "/mnt/sda1/lzh/data/inhouse/HMPV_RSV/HMPV-QZ-04/20250429/Matrix",
    "/mnt/sda1/lzh/data/inhouse/HMPV_RSV/RSV-QZ-01/20250429/Matrix",
    "/mnt/sda1/lzh/data/inhouse/HMPV_RSV/RSV-ZZ-01/20250429/Matrix",
    "/mnt/sda1/lzh/data/inhouse/HMPV_RSV/RSV-ZZ-02"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  sample_manifest,
  file = file.path(pre_table_dir, "sample_manifest.csv"),
  row.names = FALSE
)

if (anyDuplicated(sample_manifest$sample_id) > 0) {
  stop("Duplicated sample_id values were found in sample_manifest.")
}


# ============================================================
# 5. Read all samples and create Seurat objects
# ============================================================

seurat_objects <- lapply(
  seq_len(nrow(sample_manifest)),
  function(i) {
    read_10x_counts(
      path = sample_manifest$input_dir[i],
      sample_id = sample_manifest$sample_id[i]
    )
  }
)

names(seurat_objects) <- sample_manifest$sample_id

sample_summary_before_merge <- data.frame(
  sample_id = names(seurat_objects),
  genes = vapply(
    seurat_objects,
    nrow,
    numeric(1)
  ),
  cells = vapply(
    seurat_objects,
    ncol,
    numeric(1)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  sample_summary_before_merge,
  file = file.path(
    pre_table_dir,
    "sample_summary_before_merge.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 6. Merge all samples
# ============================================================

merged_seurat <- merge(
  x = seurat_objects[[1]],
  y = seurat_objects[-1],
  add.cell.ids = names(seurat_objects),
  project = "HMPV_RSV_PBMC"
)

# Explicitly define sample_id from orig.ident after merging.
merged_seurat$sample_id <- as.character(
  merged_seurat$orig.ident
)

message("Merged Seurat object:")
print(merged_seurat)

write_object_summary(
  merged_seurat,
  file.path(
    pre_log_dir,
    "01_merged_object_summary.txt"
  )
)

saveRDS(
  merged_seurat,
  file.path(
    pre_dir,
    "rawData_merged.rds"
  )
)


# ============================================================
# 7. Prepare RNA layers for Seurat v5
# ============================================================

DefaultAssay(merged_seurat) <- "RNA"

message("RNA layers before joining:")
print(Layers(merged_seurat[["RNA"]]))

if (length(Layers(merged_seurat[["RNA"]])) > 1) {
  merged_seurat[["RNA"]] <- JoinLayers(
    merged_seurat[["RNA"]]
  )
}

message("RNA layers after joining:")
print(Layers(merged_seurat[["RNA"]]))


# ============================================================
# 8. Preliminary normalization and PCA
# ============================================================

merged_seurat <- NormalizeData(
  merged_seurat,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

merged_seurat <- FindVariableFeatures(
  merged_seurat,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

merged_seurat <- ScaleData(
  merged_seurat,
  features = VariableFeatures(merged_seurat),
  verbose = FALSE
)

merged_seurat <- RunPCA(
  merged_seurat,
  features = VariableFeatures(merged_seurat),
  npcs = 50,
  verbose = FALSE
)


# ============================================================
# 9. Qualitative batch-effect evaluation
# ============================================================

p_pca_sample <- DimPlot(
  merged_seurat,
  reduction = "pca",
  group.by = "orig.ident",
  raster = TRUE
) +
  ggtitle("PCA colored by sample") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank()
  )

save_pdf(
  p_pca_sample,
  file.path(
    pre_figure_dir,
    "01_PCA_by_sample.pdf"
  ),
  width = 7,
  height = 6
)


# ============================================================
# 10. Quantitative batch-effect evaluation using kBET
# ============================================================

kbet_result <- NULL

if (requireNamespace("kBET", quietly = TRUE)) {
  message("Running kBET. This step may take a long time.")

  pca_data <- Embeddings(
    merged_seurat,
    reduction = "pca"
  )[, 1:30, drop = FALSE]

  batch_labels <- as.character(
    merged_seurat$orig.ident
  )

  set.seed(1234)

  kbet_result <- tryCatch(
    {
      kBET::kBET(
        df = pca_data,
        batch = batch_labels,
        plot = FALSE
      )
    },
    error = function(e) {
      warning("kBET failed: ", conditionMessage(e))
      NULL
    }
  )

  if (!is.null(kbet_result)) {
    saveRDS(
      kbet_result,
      file.path(
        pre_dir,
        "kBET_result.rds"
      )
    )

    writeLines(
      capture.output(print(kbet_result)),
      con = file.path(
        pre_log_dir,
        "kBET_result.txt"
      )
    )

    message(
      "kBET completed. Please inspect the full kBET output ",
      "rather than relying on a single summary value."
    )
  }
} else {
  warning(
    "The kBET package is not installed. ",
    "Quantitative batch-effect evaluation was skipped."
  )
}


# ============================================================
# 11. Calculate mitochondrial percentages
# ============================================================

# Human mitochondrial genes are generally annotated as MT-.
# If the dataset uses a different naming convention, modify the pattern.
mitochondrial_genes <- grep(
  pattern = "^MT-",
  x = rownames(merged_seurat),
  value = TRUE
)

if (length(mitochondrial_genes) == 0) {
  warning(
    "No mitochondrial genes matching '^MT-' were found. ",
    "Please check the gene naming convention."
  )
} else {
  merged_seurat[["percent.mt"]] <- PercentageFeatureSet(
    merged_seurat,
    features = mitochondrial_genes
  )
}


# ============================================================
# 12. Visualize QC metrics before filtering
# ============================================================

qc_features <- c(
  "nFeature_RNA",
  "nCount_RNA",
  "percent.mt"
)

qc_features <- qc_features[
  qc_features %in% colnames(merged_seurat@meta.data)
]

p_qc_violin <- VlnPlot(
  merged_seurat,
  features = qc_features,
  pt.size = 0.1,
  group.by = "orig.ident",
  ncol = length(qc_features)
) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

save_pdf(
  p_qc_violin,
  file.path(
    pre_figure_dir,
    "02_QC_metrics_before_filtering.pdf"
  ),
  width = 13,
  height = 5
)

p_feature_scatter <- FeatureScatter(
  merged_seurat,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA",
  group.by = "orig.ident"
) +
  geom_smooth(
    method = "lm",
    color = "red",
    se = FALSE
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  )

save_pdf(
  p_feature_scatter,
  file.path(
    pre_figure_dir,
    "03_nCount_nFeature_scatter_before_filtering.pdf"
  ),
  width = 7,
  height = 6
)


# ============================================================
# 13. Cell-level quality control
# ============================================================

min_features <- 200
max_features <- 6000
max_mt_percent <- 20

cells_before_qc <- ncol(merged_seurat)

if (!"percent.mt" %in% colnames(merged_seurat@meta.data)) {
  stop(
    "percent.mt is not available. QC filtering cannot be performed."
  )
}

merged_seurat <- subset(
  merged_seurat,
  subset =
    nFeature_RNA > min_features &
    nFeature_RNA < max_features &
    percent.mt < max_mt_percent
)

cells_after_qc <- ncol(merged_seurat)

qc_cell_summary <- data.frame(
  stage = c(
    "Before QC",
    "After QC"
  ),
  cells = c(
    cells_before_qc,
    cells_after_qc
  )
)

write.csv(
  qc_cell_summary,
  file = file.path(
    pre_table_dir,
    "cell_counts_before_after_QC.csv"
  ),
  row.names = FALSE
)

message(
  "Cells before QC: ",
  cells_before_qc
)

message(
  "Cells after QC: ",
  cells_after_qc
)


# ============================================================
# 14. Recalculate normalization and PCA after QC
# ============================================================

# Re-run normalization and dimensionality reduction after QC.
# This prevents cells removed during QC from influencing the
# downstream feature selection and PCA space.

merged_seurat <- NormalizeData(
  merged_seurat,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

merged_seurat <- FindVariableFeatures(
  merged_seurat,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

merged_seurat <- ScaleData(
  merged_seurat,
  features = VariableFeatures(merged_seurat),
  verbose = FALSE
)

merged_seurat <- RunPCA(
  merged_seurat,
  features = VariableFeatures(merged_seurat),
  npcs = 50,
  verbose = FALSE
)


# ============================================================
# 15. Save the QC-filtered object before doublet removal
# ============================================================

saveRDS(
  merged_seurat,
  file.path(
    pre_dir,
    "QC_filtered_before_doublet_removal.rds"
  )
)


# ============================================================
# 16. Doublet detection using scDblFinder
# ============================================================

message("Converting Seurat object to SingleCellExperiment.")

sce <- as.SingleCellExperiment(
  merged_seurat,
  assay = "RNA"
)

# Confirm that cell barcodes are consistent between objects.
if (!identical(colnames(sce), colnames(merged_seurat))) {
  stop(
    "Cell barcode order differs between the Seurat and ",
    "SingleCellExperiment objects."
  )
}

colData(sce)$sample_id <- as.character(
  merged_seurat$orig.ident
)

set.seed(1234)

sce <- scDblFinder(
  sce,
  samples = "sample_id"
)

required_doublet_columns <- c(
  "scDblFinder.score",
  "scDblFinder.class"
)

missing_doublet_columns <- setdiff(
  required_doublet_columns,
  colnames(colData(sce))
)

if (length(missing_doublet_columns) > 0) {
  stop(
    "scDblFinder did not generate the following columns: ",
    paste(missing_doublet_columns, collapse = ", ")
  )
}

doublet_metadata <- as.data.frame(
  colData(sce)
)[
  ,
  required_doublet_columns,
  drop = FALSE
]

rownames(doublet_metadata) <- colnames(sce)

# Match doublet results explicitly by cell barcode.
merged_seurat$scDblFinder.score <- doublet_metadata[
  colnames(merged_seurat),
  "scDblFinder.score"
]

merged_seurat$scDblFinder.class <- doublet_metadata[
  colnames(merged_seurat),
  "scDblFinder.class"
]

doublet_summary <- merged_seurat@meta.data %>%
  as.data.frame() %>%
  count(
    orig.ident,
    scDblFinder.class,
    name = "cells"
  ) %>%
  group_by(orig.ident) %>%
  mutate(
    proportion = cells / sum(cells)
  ) %>%
  ungroup()

write.csv(
  doublet_summary,
  file = file.path(
    pre_table_dir,
    "doublet_summary_by_sample.csv"
  ),
  row.names = FALSE
)

print(
  table(
    merged_seurat$orig.ident,
    merged_seurat$scDblFinder.class,
    useNA = "ifany"
  )
)


# ============================================================
# 17. Remove predicted doublets
# ============================================================

cells_before_doublet_removal <- ncol(merged_seurat)

merged_seurat <- subset(
  merged_seurat,
  subset = scDblFinder.class == "singlet"
)

cells_after_doublet_removal <- ncol(merged_seurat)

doublet_filter_summary <- data.frame(
  stage = c(
    "Before doublet removal",
    "After doublet removal"
  ),
  cells = c(
    cells_before_doublet_removal,
    cells_after_doublet_removal
  )
)

write.csv(
  doublet_filter_summary,
  file = file.path(
    pre_table_dir,
    "cell_counts_before_after_doublet_removal.csv"
  ),
  row.names = FALSE
)

message(
  "Cells before doublet removal: ",
  cells_before_doublet_removal
)

message(
  "Cells after doublet removal: ",
  cells_after_doublet_removal
)


# ============================================================
# 18. Save the final QC-filtered object
# ============================================================

clean_data_file <- file.path(
  pre_dir,
  "cleanData.rds"
)

saveRDS(
  merged_seurat,
  clean_data_file
)

write_object_summary(
  merged_seurat,
  file.path(
    pre_log_dir,
    "02_cleanData_object_summary.txt"
  )
)


# ============================================================
# 19. Final normalization and PCA after doublet removal
# ============================================================

merged_seurat <- NormalizeData(
  merged_seurat,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

merged_seurat <- FindVariableFeatures(
  merged_seurat,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

merged_seurat <- ScaleData(
  merged_seurat,
  features = VariableFeatures(merged_seurat),
  verbose = FALSE
)

merged_seurat <- RunPCA(
  merged_seurat,
  features = VariableFeatures(merged_seurat),
  npcs = 50,
  verbose = FALSE
)


# ============================================================
# 20. Construct the neighbor graph
# ============================================================

merged_seurat <- FindNeighbors(
  merged_seurat,
  reduction = "pca",
  dims = 1:50,
  graph.name = c(
    "myKNN_nn",
    "myKNN_snn"
  ),
  verbose = FALSE
)


# ============================================================
# 21. Multi-resolution clustering
# ============================================================

clustering_resolutions <- seq(
  from = 0.2,
  to = 1.5,
  by = 0.1
)

merged_seurat <- FindClusters(
  merged_seurat,
  graph.name = "myKNN_snn",
  resolution = clustering_resolutions,
  algorithm = 1,
  random.seed = 1234,
  verbose = FALSE
)

clustering_columns <- grep(
  "^myKNN_snn_res\\.",
  colnames(merged_seurat@meta.data),
  value = TRUE
)

message("Generated clustering columns:")
print(clustering_columns)


# ============================================================
# 22. Clustering tree
# ============================================================

p_clustree <- clustree(
  merged_seurat@meta.data,
  prefix = "myKNN_snn_res.",
  node_color = "purple"
) +
  ggtitle("Clustering tree by resolution") +
  theme_bw()

save_pdf(
  p_clustree,
  file.path(
    pre_figure_dir,
    "04_clustering_tree.pdf"
  ),
  width = 10,
  height = 8
)


# ============================================================
# 23. Select the final clustering resolution
# ============================================================

# Change this value if another resolution is selected after
# inspecting the clustering tree and marker genes.
final_resolution <- 0.6

final_cluster_column <- paste0(
  "myKNN_snn_res.",
  final_resolution
)

if (!final_cluster_column %in% colnames(merged_seurat@meta.data)) {
  stop(
    "The selected final clustering column was not found: ",
    final_cluster_column
  )
}

merged_seurat$seurat_clusters <- merged_seurat@meta.data[
  [final_cluster_column]
]

Idents(merged_seurat) <- "seurat_clusters"

message("Final cluster sizes:")
print(table(merged_seurat$seurat_clusters))


# ============================================================
# 24. Run UMAP using the PCA coordinates
# ============================================================

merged_seurat <- RunUMAP(
  merged_seurat,
  reduction = "pca",
  dims = 1:50,
  reduction.name = "umap.pca",
  reduction.key = "PCUMAP_",
  seed.use = 1234,
  verbose = FALSE
)


# ============================================================
# 25. Generate final clustering plots
# ============================================================

p_cluster <- DimPlot(
  merged_seurat,
  reduction = "umap.pca",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  label.size = 5,
  raster = TRUE
) +
  ggtitle(
    paste0(
      "Final clustering, resolution = ",
      final_resolution
    )
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank()
  )

save_pdf(
  p_cluster,
  file.path(
    pre_figure_dir,
    "05_final_clustering_UMAP.pdf"
  ),
  width = 8,
  height = 6
)

p_cluster_sample <- DimPlot(
  merged_seurat,
  reduction = "umap.pca",
  group.by = "orig.ident",
  raster = TRUE
) +
  ggtitle("Final UMAP colored by sample") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank()
  )

save_pdf(
  p_cluster_sample,
  file.path(
    pre_figure_dir,
    "06_final_UMAP_by_sample.pdf"
  ),
  width = 8,
  height = 6
)


# ============================================================
# 26. Save the final preprocessing object
# ============================================================

final_object_file <- file.path(
  pre_dir,
  "cluster.rds"
)

saveRDS(
  merged_seurat,
  final_object_file
)

write_object_summary(
  merged_seurat,
  file.path(
    pre_log_dir,
    "03_cluster_object_summary.txt"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    pre_log_dir,
    "preprocessing_sessionInfo.txt"
  )
)

message("Preprocessing pipeline completed successfully.")
message("Final object saved to: ", final_object_file)
