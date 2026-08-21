# ============================================================
# Figure 4: Global PBMC preprocessing, annotation, and
# composition analysis
#
# Project:
#   HMPV single-cell RNA-seq analysis
#
# Author:
#   lzh
#
# Date:
#   2026-08-21
#
# Description:
#   This script performs the following analyses:
#   1. Load the QC-filtered Seurat object
#   2. Validate sample metadata
#   3. Perform normalization, variable feature selection,
#      scaling, PCA, Harmony integration, UMAP, and clustering
#   4. Annotate major immune cell types
#   5. Generate marker plots and UMAP plots
#   6. Calculate cell-type composition by group and sample
#   7. Compare scRNA-seq neutrophil proportions with CBC results
#   8. Compare NK-cell proportions between RSV mild and severe groups
#
# Note:
#   This script is intended for Seurat v5.
# ============================================================


# ============================================================
# 1. Environment setup
# ============================================================

.libPaths("/home/lzh/miniconda3/envs/seurat5/lib/R/library")

set.seed(1234)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(harmony)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)
  library(scales)
  library(RColorBrewer)
  library(ggrepel)
  library(ggrastr)
})

message("Seurat version: ", packageVersion("Seurat"))
message("SeuratObject version: ", packageVersion("SeuratObject"))
message("ggplot2 version: ", packageVersion("ggplot2"))
message("patchwork version: ", packageVersion("patchwork"))
message("Matrix version: ", packageVersion("Matrix"))


# ============================================================
# 2. Project paths
# ============================================================

project_dir <- "/mnt/sda1/lzh/test/cjl/hmpv"

pre_dir <- file.path(project_dir, "pre")
metadata_dir <- file.path(project_dir, "metadata")
common_dir <- file.path(project_dir, "common")

fig4_dir <- file.path(project_dir, "Figure4_global_PBMC")
fig4_object_dir <- file.path(fig4_dir, "objects")
fig4_figure_dir <- file.path(fig4_dir, "figures")
fig4_table_dir <- file.path(fig4_dir, "tables")
fig4_log_dir <- file.path(fig4_dir, "logs")

required_directories <- c(
  pre_dir,
  metadata_dir,
  common_dir,
  fig4_object_dir,
  fig4_figure_dir,
  fig4_table_dir,
  fig4_log_dir
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
# 3. General plotting themes and colors
# ============================================================

theme_global <- theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    aspect.ratio = 1
  )

group_levels <- c(
  "Control",
  "HMPVmild",
  "HMPVsevere",
  "RSVmild",
  "RSVsevere"
)

virus_levels <- c(
  "Control",
  "HMPV",
  "RSV"
)

severity_levels <- c(
  "Control",
  "mild",
  "severe"
)

major_celltype_levels <- c(
  "B cell",
  "CD4+ T cell",
  "CD8+ T cell",
  "HSC-G-CSF",
  "Monocyte",
  "Neutrophil",
  "NK cell",
  "Unassigned"
)

group_colors <- c(
  "Control" = "#808080",
  "HMPVmild" = "#B2DF8A",
  "HMPVsevere" = "#78B478",
  "RSVmild" = "#74ADD1",
  "RSVsevere" = "#3278C8"
)

celltype_colors <- c(
  "B cell" = "#FDBF6F",
  "CD4+ T cell" = "#A6CEE3",
  "CD8+ T cell" = "#74ADD1",
  "HSC-G-CSF" = "#33A02C",
  "Monocyte" = "#DECBE4",
  "Neutrophil" = "#BFA4CB",
  "NK cell" = "#1F78B4",
  "Unassigned" = "#FB9A99"
)


# ============================================================
# 4. Utility functions
# ============================================================

check_required_columns <- function(data, required_columns, object_name) {
  missing_columns <- setdiff(required_columns, colnames(data))

  if (length(missing_columns) > 0) {
    stop(
      object_name,
      " is missing the following columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  invisible(TRUE)
}


save_plot_pdf <- function(plot_object, file, width, height) {
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


get_qc_column <- function(metadata, candidates) {
  existing_columns <- candidates[candidates %in% colnames(metadata)]

  if (length(existing_columns) == 0) {
    return(NULL)
  }

  existing_columns[1]
}


# ============================================================
# 5. Locate and load the QC-filtered Seurat object
# ============================================================

rds_files <- list.files(
  path = pre_dir,
  pattern = "\\.(rds|RDS)$",
  full.names = TRUE
)

if (length(rds_files) == 0) {
  stop("No RDS file was found in: ", pre_dir)
}

if (length(rds_files) > 1) {
  stop(
    "More than one RDS file was found in ", pre_dir,
    ". Please specify the input file manually."
  )
}

input_rds <- rds_files[1]

message("Input RDS file: ", input_rds)

obj0 <- readRDS(input_rds)

if (inherits(obj0, "Seurat")) {
  pbmc <- obj0
} else if (
  is.list(obj0) &&
  length(obj0) > 0 &&
  all(vapply(obj0, inherits, logical(1), what = "Seurat"))
) {
  pbmc <- merge(
    x = obj0[[1]],
    y = obj0[-1],
    add.cell.ids = names(obj0)
  )
} else {
  stop(
    "The RDS file does not contain a Seurat object or a list of Seurat objects."
  )
}

message("Loaded object:")
print(pbmc)

writeLines(
  capture.output(print(pbmc)),
  con = file.path(fig4_log_dir, "Fig4_input_object_summary.txt")
)


# ============================================================
# 6. Basic object validation
# ============================================================

if (!"RNA" %in% Assays(pbmc)) {
  stop("The Seurat object does not contain an RNA assay.")
}

DefaultAssay(pbmc) <- "RNA"

message("Assays:")
print(Assays(pbmc))

message("Default assay: ", DefaultAssay(pbmc))

message("Available reductions:")
print(Reductions(pbmc))

message("Available graphs:")
print(names(pbmc@graphs))

message("Object dimensions:")
print(dim(pbmc))


# ============================================================
# 7. Clean and standardize cell-level metadata
# ============================================================

metadata <- pbmc@meta.data

# Remove obsolete custom KNN columns.
metadata <- metadata[
  ,
  !str_detect(colnames(metadata), "^myKNN_"),
  drop = FALSE
]

pbmc@meta.data <- metadata

# Generate sample_id if it is not already available.
if (!"sample_id" %in% colnames(pbmc@meta.data)) {
  if ("orig.ident" %in% colnames(pbmc@meta.data)) {
    pbmc$sample_id <- as.character(pbmc$orig.ident)
  } else {
    pbmc$sample_id <- sub(
      pattern = "_.*$",
      replacement = "",
      x = colnames(pbmc)
    )
  }
}

pbmc$sample_id <- as.character(pbmc$sample_id)

if (anyNA(pbmc$sample_id) || any(pbmc$sample_id == "")) {
  stop("Some cells have missing or empty sample_id values.")
}

sample_detected <- pbmc@meta.data %>%
  as.data.frame() %>%
  count(sample_id, name = "cells") %>%
  arrange(sample_id)

write.csv(
  sample_detected,
  file = file.path(fig4_table_dir, "Fig4_detected_samples.csv"),
  row.names = FALSE
)

message("Detected samples:")
print(sample_detected)


# ============================================================
# 8. Load and validate sample metadata
# ============================================================

sample_metadata_file <- file.path(
  metadata_dir,
  "sample_metadata.csv"
)

if (!file.exists(sample_metadata_file)) {
  stop("Sample metadata file was not found: ", sample_metadata_file)
}

sample_meta <- read.csv(
  sample_metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

check_required_columns(
  data = sample_meta,
  required_columns = c(
    "sample_id",
    "group",
    "virus",
    "severity"
  ),
  object_name = "sample_metadata.csv"
)

sample_meta <- sample_meta %>%
  mutate(
    sample_id = as.character(sample_id),
    group = as.character(group),
    virus = as.character(virus),
    severity = as.character(severity)
  )

missing_in_metadata <- setdiff(
  unique(pbmc$sample_id),
  unique(sample_meta$sample_id)
)

if (length(missing_in_metadata) > 0) {
  stop(
    "The following sample_id values in the Seurat object are missing ",
    "from sample_metadata.csv: ",
    paste(missing_in_metadata, collapse = ", ")
  )
}

unused_metadata_samples <- setdiff(
  unique(sample_meta$sample_id),
  unique(pbmc$sample_id)
)

if (length(unused_metadata_samples) > 0) {
  warning(
    "The following metadata samples are not present in the Seurat object: ",
    paste(unused_metadata_samples, collapse = ", ")
  )
}

sample_meta_used <- sample_meta %>%
  filter(sample_id %in% unique(pbmc$sample_id))

sample_meta_used <- sample_meta_used[
  match(unique(pbmc$sample_id), sample_meta_used$sample_id),
  ,
  drop = FALSE
]

sample_index <- match(
  pbmc$sample_id,
  sample_meta_used$sample_id
)

pbmc$group <- sample_meta_used$group[sample_index]
pbmc$virus <- sample_meta_used$virus[sample_index]
pbmc$severity <- sample_meta_used$severity[sample_index]

pbmc$group <- factor(
  pbmc$group,
  levels = group_levels
)

pbmc$virus <- factor(
  pbmc$virus,
  levels = virus_levels
)

pbmc$severity <- factor(
  pbmc$severity,
  levels = severity_levels
)

if (anyNA(pbmc$group)) {
  warning("Some cells have NA values in group.")
}

message("Cells by group:")
print(table(pbmc$group, useNA = "ifany"))

message("Cells by virus:")
print(table(pbmc$virus, useNA = "ifany"))

message("Cells by severity:")
print(table(pbmc$severity, useNA = "ifany"))


# ============================================================
# 9. Save metadata summary
# ============================================================

sample_summary <- pbmc@meta.data %>%
  as.data.frame() %>%
  count(
    sample_id,
    group,
    virus,
    severity,
    name = "cells"
  ) %>%
  arrange(group, sample_id)

write.csv(
  sample_summary,
  file = file.path(fig4_table_dir, "Fig4_cells_by_sample_group.csv"),
  row.names = FALSE
)

message("Total number of cells: ", ncol(pbmc))


# ============================================================
# 10. QC summary
# ============================================================

metadata_df <- pbmc@meta.data %>%
  as.data.frame()

ncount_column <- get_qc_column(
  metadata_df,
  c("nCount_RNA", "nCount_SCT")
)

nfeature_column <- get_qc_column(
  metadata_df,
  c("nFeature_RNA", "nFeature_SCT")
)

mitochondrial_column <- get_qc_column(
  metadata_df,
  c(
    "percent.mt",
    "percent.MT",
    "mito.percent",
    "percent_mito"
  )
)

qc_summary <- metadata_df %>%
  group_by(sample_id, group) %>%
  summarise(
    cells = n(),
    median_nCount = if (!is.null(ncount_column)) {
      median(.data[[ncount_column]], na.rm = TRUE)
    } else {
      NA_real_
    },
    median_nFeature = if (!is.null(nfeature_column)) {
      median(.data[[nfeature_column]], na.rm = TRUE)
    } else {
      NA_real_
    },
    median_percent_mt = if (!is.null(mitochondrial_column)) {
      median(.data[[mitochondrial_column]], na.rm = TRUE)
    } else {
      NA_real_
    },
    .groups = "drop"
  )

write.csv(
  qc_summary,
  file = file.path(fig4_table_dir, "Fig4_QC_summary_by_sample.csv"),
  row.names = FALSE
)


# ============================================================
# 11. QC and sample-size plots
# ============================================================

p_cells <- ggplot(
  sample_summary,
  aes(
    x = sample_id,
    y = cells,
    fill = group
  )
) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  labs(
    x = NULL,
    y = "Number of cells",
    title = "Number of cells per sample"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.title = element_blank()
  )

save_plot_pdf(
  p_cells,
  file.path(fig4_figure_dir, "Fig4_cells_per_sample.pdf"),
  width = 10,
  height = 5
)

if (!is.null(ncount_column) && !is.null(nfeature_column)) {
  qc_features <- c(ncount_column, nfeature_column)

  p_qc <- VlnPlot(
    pbmc,
    features = qc_features,
    group.by = "sample_id",
    pt.size = 0,
    ncol = 1
  ) &
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    )

  save_plot_pdf(
    p_qc,
    file.path(fig4_figure_dir, "Fig4_QC_nCount_nFeature_by_sample.pdf"),
    width = 12,
    height = 8
  )
}

if (!is.null(mitochondrial_column)) {
  p_mito <- VlnPlot(
    pbmc,
    features = mitochondrial_column,
    group.by = "sample_id",
    pt.size = 0
  ) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    )

  save_plot_pdf(
    p_mito,
    file.path(fig4_figure_dir, "Fig4_QC_mitochondrial_content_by_sample.pdf"),
    width = 12,
    height = 5
  )
}


# ============================================================
# 12. Save object with sample metadata
# ============================================================

metadata_added_file <- file.path(
  fig4_object_dir,
  "Fig4_01_metadata_added.rds"
)

saveRDS(pbmc, metadata_added_file)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    fig4_log_dir,
    "Fig4_01_metadata_added_sessionInfo.txt"
  )
)


# ============================================================
# 13. Prepare Seurat v5 RNA layers
# ============================================================

DefaultAssay(pbmc) <- "RNA"

rna_layers <- Layers(pbmc[["RNA"]])
message("RNA layers before preparation:")
print(rna_layers)

if (length(rna_layers) > 1) {
  pbmc[["RNA"]] <- JoinLayers(pbmc[["RNA"]])
}

pbmc[["RNA"]] <- split(
  pbmc[["RNA"]],
  f = pbmc$sample_id
)

message("RNA layers after sample-based splitting:")
print(Layers(pbmc[["RNA"]]))


# ============================================================
# 14. Normalization, variable features, scaling, and PCA
# ============================================================

pbmc <- NormalizeData(
  object = pbmc,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

pbmc <- FindVariableFeatures(
  object = pbmc,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

pbmc <- ScaleData(
  object = pbmc,
  features = VariableFeatures(pbmc),
  verbose = FALSE
)

pbmc <- RunPCA(
  object = pbmc,
  features = VariableFeatures(pbmc),
  npcs = 50,
  verbose = FALSE
)

p_elbow <- ElbowPlot(
  pbmc,
  ndims = 50
) +
  ggtitle("PCA elbow plot")

save_plot_pdf(
  p_elbow,
  file.path(fig4_figure_dir, "Fig4_PCA_elbow_plot.pdf"),
  width = 6,
  height = 5
)


# ============================================================
# 15. Harmony integration
# ============================================================

pbmc <- IntegrateLayers(
  object = pbmc,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = TRUE
)

if (!"harmony" %in% Reductions(pbmc)) {
  stop("Harmony integration failed: harmony reduction was not found.")
}

message("Harmony embedding dimensions:")
print(dim(Embeddings(pbmc, reduction = "harmony")))


# ============================================================
# 16. Neighbor graph, UMAP, and clustering
# ============================================================

pbmc <- FindNeighbors(
  object = pbmc,
  reduction = "harmony",
  dims = 1:50,
  graph.name = c(
    "harmony_nn",
    "harmony_snn"
  ),
  verbose = FALSE
)

pbmc <- RunUMAP(
  object = pbmc,
  reduction = "harmony",
  dims = 1:50,
  reduction.name = "umap.harmony",
  reduction.key = "harmonyUMAP_",
  seed.use = 1234,
  verbose = FALSE
)

pbmc <- FindClusters(
  object = pbmc,
  graph.name = "harmony_snn",
  resolution = c(0.2, 0.3, 0.4, 0.5, 0.6),
  algorithm = 1,
  random.seed = 1234,
  verbose = FALSE
)

resolution_columns <- grep(
  "^harmony_snn_res\\.",
  colnames(pbmc@meta.data),
  value = TRUE
)

message("Available clustering resolutions:")
print(resolution_columns)


# ============================================================
# 17. Compare clustering resolutions
# ============================================================

resolution_plots <- lapply(
  c("0.2", "0.3", "0.4", "0.5", "0.6"),
  function(resolution) {
    cluster_column <- paste0("harmony_snn_res.", resolution)

    DimPlot(
      pbmc,
      reduction = "umap.harmony",
      group.by = cluster_column,
      label = TRUE,
      repel = TRUE,
      raster = TRUE
    ) +
      ggtitle(paste0("Resolution = ", resolution)) +
      theme_global
  }
)

p_resolution_comparison <- wrap_plots(
  resolution_plots,
  ncol = 2
)

save_plot_pdf(
  p_resolution_comparison,
  file.path(fig4_figure_dir, "Fig4_clustering_resolution_comparison.pdf"),
  width = 12,
  height = 10
)


# ============================================================
# 18. Select the final clustering resolution
# ============================================================

final_cluster_column <- "harmony_snn_res.0.6"

if (!final_cluster_column %in% colnames(pbmc@meta.data)) {
  stop(
    "The selected clustering column was not found: ",
    final_cluster_column
  )
}

pbmc$seurat_clusters <- pbmc@meta.data[[final_cluster_column]]

Idents(pbmc) <- "seurat_clusters"

message("Final cluster sizes:")
print(table(pbmc$seurat_clusters))


# ============================================================
# 19. Global cluster UMAP
# ============================================================

p_cluster_harmony <- DimPlot(
  pbmc,
  reduction = "umap.harmony",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  raster = FALSE,
  pt.size = 0.8
) +
  theme_global +
  ggtitle("Global PBMC clusters")

save_plot_pdf(
  p_cluster_harmony,
  file.path(fig4_figure_dir, "Fig4_Harmony_UMAP_cluster.pdf"),
  width = 8,
  height = 6
)


# ============================================================
# 20. Marker genes for global cluster annotation
# ============================================================

global_marker_genes <- c(
  # T cells
  "CD3D",
  "CD3E",
  "TRAC",
  "TRBC1",
  "TRBC2",

  # CD4 T cells
  "CD4",
  "IL7R",
  "LTB",
  "CCR7",

  # CD8 T cells
  "CD8A",
  "CD8B",
  "GZMK",
  "GZMH",

  # NK cells
  "NKG7",
  "GNLY",
  "PRF1",
  "KLRD1",
  "KLRF1",

  # B cells
  "MS4A1",
  "CD79A",
  "CD79B",
  "CD74",
  "HLA-DRA",

  # Plasma cells
  "MZB1",
  "JCHAIN",
  "SDC1",
  "CD38",
  "IGHG1",
  "IGKC",

  # Monocytes
  "LYZ",
  "LST1",
  "FCN1",
  "CTSS",
  "FCGR3A",

  # Neutrophils
  "FCGR3B",
  "CSF3R",
  "CXCR2",
  "CEACAM8",
  "S100A8",
  "S100A9",
  "MPO",
  "ELANE",
  "AZU1",

  # Platelets and megakaryocytes
  "PPBP",
  "PF4",
  "NRGN",
  "GP9",
  "ITGA2B",
  "GP1BA",

  # Erythroid cells
  "HBB",
  "HBA1",
  "HBA2",
  "AHSP",
  "ALAS2",
  "GYPA",

  # HSC/progenitor-like cells
  "KIT",
  "CD34"
)

global_marker_genes <- intersect(
  global_marker_genes,
  rownames(pbmc)
)

p_dot_cluster <- DotPlot(
  pbmc,
  features = global_marker_genes,
  group.by = "seurat_clusters",
  assay = "RNA"
) +
  RotatedAxis() +
  scale_color_distiller(
    palette = "RdYlBu"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  ) +
  labs(
    x = NULL,
    y = NULL
  )

save_plot_pdf(
  p_dot_cluster,
  file.path(
    fig4_figure_dir,
    "Fig4_Harmony_marker_DotPlot_cluster.pdf"
  ),
  width = 14,
  height = 7
)


# ============================================================
# 21. Find cluster marker genes
# ============================================================

pbmc_marker <- pbmc

DefaultAssay(pbmc_marker) <- "RNA"

pbmc_marker[["RNA"]] <- JoinLayers(
  pbmc_marker[["RNA"]]
)

Idents(pbmc_marker) <- "seurat_clusters"

cluster_markers <- FindAllMarkers(
  object = pbmc_marker,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.3,
  logfc.threshold = 0.25
)

cluster_markers <- cluster_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  ungroup()

write.csv(
  cluster_markers,
  file = file.path(
    fig4_table_dir,
    "Fig4_Harmony_all_cluster_markers.csv"
  ),
  row.names = FALSE
)

cluster_top20 <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 20,
    with_ties = FALSE
  ) %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  ungroup()

write.csv(
  cluster_top20,
  file = file.path(
    fig4_table_dir,
    "Fig4_Harmony_cluster_top20_markers.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 22. Manual major cell-type annotation
# ============================================================

cluster_celltype <- c(
  "0" = "Neutrophil",
  "1" = "CD4+ T cell",
  "2" = "B cell",
  "3" = "Monocyte",
  "4" = "NK cell",
  "5" = "Neutrophil",
  "6" = "CD8+ T cell",
  "7" = "NK cell",
  "8" = "CD4+ T cell",
  "9" = "CD4+ T cell",
  "10" = "B cell",
  "11" = "Monocyte",
  "12" = "CD8+ T cell",
  "13" = "Neutrophil",
  "14" = "CD4+ T cell",
  "15" = "Monocyte",
  "16" = "Monocyte",
  "17" = "B cell",
  "18" = "Unassigned",
  "19" = "B cell",
  "20" = "B cell",
  "21" = "Unassigned",
  "22" = "Unassigned",
  "23" = "Unassigned",
  "24" = "Neutrophil",
  "25" = "HSC-G-CSF",
  "26" = "Unassigned"
)

all_clusters <- sort(
  unique(as.character(pbmc$seurat_clusters))
)

missing_annotations <- setdiff(
  all_clusters,
  names(cluster_celltype)
)

if (length(missing_annotations) > 0) {
  stop(
    "The following clusters do not have annotations: ",
    paste(missing_annotations, collapse = ", ")
  )
}

pbmc$major_celltype <- unname(
  cluster_celltype[as.character(pbmc$seurat_clusters)]
)

pbmc$major_celltype <- factor(
  pbmc$major_celltype,
  levels = major_celltype_levels
)

Idents(pbmc) <- "major_celltype"

message("Major cell-type composition:")
print(table(pbmc$major_celltype, useNA = "ifany"))

if (anyNA(pbmc$major_celltype)) {
  warning("Some cells have NA values in major_celltype.")
}


# ============================================================
# 23. Major cell-type UMAPs
# ============================================================

p_major_umap <- DimPlot(
  pbmc,
  reduction = "umap.harmony",
  group.by = "major_celltype",
  raster = FALSE,
  label = TRUE,
  repel = TRUE,
  pt.size = 0.8,
  cols = celltype_colors
) +
  theme_global +
  theme(legend.title = element_blank()) +
  ggtitle("Major cell types")

p_group_umap <- DimPlot(
  pbmc,
  reduction = "umap.harmony",
  group.by = "group",
  raster = FALSE,
  repel = TRUE,
  pt.size = 0.8,
  cols = group_colors
) +
  theme_global +
  theme(legend.title = element_blank()) +
  ggtitle("Clinical groups")

save_plot_pdf(
  p_major_umap,
  file.path(fig4_figure_dir, "Fig4_major_celltype_UMAP.pdf"),
  width = 7,
  height = 6
)

save_plot_pdf(
  p_group_umap,
  file.path(fig4_figure_dir, "Fig4_group_UMAP.pdf"),
  width = 7,
  height = 6
)


# ============================================================
# 24. Major cell-type marker DotPlot
# ============================================================

major_marker_genes <- c(
  "MS4A1",
  "CD79A",
  "IGHG1",
  "CD3D",
  "IL7R",
  "LTB",
  "CD8A",
  "CD8B",
  "GZMK",
  "KIT",
  "CD34",
  "LYZ",
  "FCN1",
  "CTSS",
  "FCGR3B",
  "CSF3R",
  "CXCR2",
  "GNLY",
  "KLRD1",
  "PRF1"
)

major_marker_genes <- intersect(
  major_marker_genes,
  rownames(pbmc)
)

p_major_dotplot <- DotPlot(
  pbmc,
  features = major_marker_genes,
  group.by = "major_celltype",
  dot.min = 0,
  dot.scale = 6,
  scale = TRUE
) +
  scale_color_distiller(palette = "RdYlBu") +
  scale_y_discrete(limits = major_celltype_levels) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    ),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  ) +
  labs(
    x = NULL,
    y = NULL
  )

save_plot_pdf(
  p_major_dotplot,
  file.path(
    fig4_figure_dir,
    "Fig4_major_celltype_marker_DotPlot.pdf"
  ),
  width = 8,
  height = 5
)


# ============================================================
# 25. Marker FeaturePlots
# ============================================================

pbmc_plot <- pbmc

DefaultAssay(pbmc_plot) <- "RNA"

if (length(Layers(pbmc_plot[["RNA"]])) > 1) {
  pbmc_plot[["RNA"]] <- JoinLayers(
    pbmc_plot[["RNA"]]
  )
}

pdf(
  file = file.path(
    fig4_figure_dir,
    "Fig4_marker_UMAP.pdf"
  ),
  width = 6,
  height = 6
)

for (gene in major_marker_genes) {
  plot_result <- tryCatch(
    {
      FeaturePlot(
        pbmc_plot,
        features = gene,
        reduction = "umap.harmony",
        raster = TRUE,
        order = TRUE
      ) +
        theme_global +
        scale_color_distiller(palette = "RdYlBu") +
        theme(
          panel.border = element_rect(
            color = "black",
            fill = NA,
            linewidth = 0.5
          ),
          plot.title = element_text(hjust = 0.5)
        )
    },
    error = function(e) {
      warning("Failed to plot gene: ", gene)
      NULL
    }
  )

  if (!is.null(plot_result)) {
    print(plot_result)
  }
}

dev.off()


# ============================================================
# 26. Average-expression heatmap
# ============================================================

Idents(pbmc_plot) <- "major_celltype"

average_expression <- AverageExpression(
  pbmc_plot,
  assays = "RNA",
  features = major_marker_genes,
  group.by = "major_celltype",
  return.seurat = TRUE
)

DefaultAssay(average_expression) <- "RNA"

average_expression <- ScaleData(
  average_expression,
  features = major_marker_genes,
  verbose = FALSE
)

available_celltypes <- intersect(
  major_celltype_levels,
  colnames(average_expression)
)

average_expression <- average_expression[
  ,
  available_celltypes
]

heatmap_palette <- colorRampPalette(
  brewer.pal(11, "RdYlBu")
)(256)

p_heatmap_mean <- DoHeatmap(
  average_expression,
  features = major_marker_genes,
  slot = "scale.data",
  group.colors = celltype_colors[available_celltypes],
  size = 4,
  draw.lines = FALSE
) +
  scale_fill_gradientn(
    colours = rev(heatmap_palette)
  ) +
  labs(
    x = NULL,
    y = NULL
  )

save_plot_pdf(
  p_heatmap_mean,
  file.path(
    fig4_figure_dir,
    "Fig4_major_celltype_marker_mean_DoHeatmap.pdf"
  ),
  width = 7,
  height = 7
)


# ============================================================
# 27. Save the annotated global object
# ============================================================

annotated_object_file <- file.path(
  fig4_object_dir,
  "Fig4_03_major_celltype_harmony.rds"
)

saveRDS(
  pbmc,
  annotated_object_file
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    fig4_log_dir,
    "Fig4_03_major_celltype_harmony_sessionInfo.txt"
  )
)


# ============================================================
# 28. Cell-type composition by clinical group
# ============================================================

cell_count_group <- pbmc@meta.data %>%
  as.data.frame() %>%
  count(
    group,
    major_celltype,
    name = "cell_number"
  ) %>%
  complete(
    group = group_levels,
    major_celltype = major_celltype_levels,
    fill = list(cell_number = 0)
  )

cell_prop_group <- cell_count_group %>%
  group_by(group) %>%
  mutate(
    total_cells = sum(cell_number),
    proportion = cell_number / total_cells
  ) %>%
  ungroup()

write.csv(
  cell_prop_group,
  file = file.path(
    fig4_table_dir,
    "Fig4_cell_proportion_by_group_long.csv"
  ),
  row.names = FALSE
)

p_cell_prop_group <- ggplot(
  cell_prop_group,
  aes(
    x = group,
    y = proportion,
    fill = major_celltype
  )
) +
  geom_col(
    width = 0.85,
    color = "white",
    linewidth = 0
  ) +
  scale_fill_manual(
    values = celltype_colors,
    drop = FALSE
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.03, 0.05))
  ) +
  labs(
    x = NULL,
    y = "Cell proportion",
    fill = NULL
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(color = "black"),
    legend.title = element_blank()
  )

save_plot_pdf(
  p_cell_prop_group,
  file.path(
    fig4_figure_dir,
    "Fig4_celltype_proportion_by_group.pdf"
  ),
  width = 6,
  height = 5.5
)


# ============================================================
# 29. Cell-type composition by sample
# ============================================================

cell_count_sample <- pbmc@meta.data %>%
  as.data.frame() %>%
  count(
    sample_id,
    group,
    major_celltype,
    name = "cell_number"
  ) %>%
  complete(
    sample_id,
    group,
    major_celltype = major_celltype_levels,
    fill = list(cell_number = 0)
  )

cell_count_sample_wide <- cell_count_sample %>%
  pivot_wider(
    names_from = major_celltype,
    values_from = cell_number,
    values_fill = 0
  ) %>%
  relocate(
    group,
    .after = last_col()
  )

write.csv(
  cell_count_sample_wide,
  file = file.path(
    fig4_table_dir,
    "Fig4_cell_number_by_sample.csv"
  ),
  row.names = FALSE
)

cell_prop_sample <- cell_count_sample %>%
  group_by(sample_id) %>%
  mutate(
    total_cells = sum(cell_number),
    proportion = cell_number / total_cells
  ) %>%
  ungroup()

cell_prop_sample_wide <- cell_prop_sample %>%
  select(
    sample_id,
    group,
    major_celltype,
    proportion
  ) %>%
  pivot_wider(
    names_from = major_celltype,
    values_from = proportion,
    values_fill = 0
  ) %>%
  relocate(
    group,
    .after = last_col()
  )

cell_prop_sample_percent <- cell_prop_sample_wide

celltype_columns <- setdiff(
  colnames(cell_prop_sample_percent),
  c("sample_id", "group")
)

cell_prop_sample_percent[
  ,
  celltype_columns
] <- cell_prop_sample_percent[
  ,
  celltype_columns
] * 100

write.csv(
  cell_prop_sample_percent,
  file = file.path(
    fig4_table_dir,
    "Fig4_cell_proportion_percent_by_sample.csv"
  ),
  row.names = FALSE
)



# ============================================================
# 30. Save final session information
# ============================================================

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    fig4_log_dir,
    "Fig4_global_PBMC_sessionInfo.txt"
  )
)

message("Figure 4 analysis completed successfully.")
