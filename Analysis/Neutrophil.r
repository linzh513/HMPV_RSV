# Figure 5: Neutrophil subcluster analysis, differential expression,
# pathway enrichment, pseudotime analysis, and reviewer-requested analyses

.libPaths("/home/lzh/miniconda3/envs/seurat5/lib/R/library")
set.seed(1234)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(harmony)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(scales)
  library(data.table)
  library(RColorBrewer)
  library(ggrepel)
  library(ggrastr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

project_dir <- Sys.getenv("HMPV_PROJECT_DIR", unset = "/mnt/sda1/lzh/test/cjl/hmpv")
fig4_object_dir <- file.path(project_dir, "Figure4_global_PBMC", "objects")
fig5_dir <- file.path(project_dir, "Figure5_neutrophil")
fig5_object_dir <- file.path(fig5_dir, "objects")
fig5_figure_dir <- file.path(fig5_dir, "figures")
fig5_table_dir <- file.path(fig5_dir, "tables")
fig5_log_dir <- file.path(fig5_dir, "logs")

required_directories <- c(fig5_object_dir, fig5_figure_dir, fig5_table_dir, fig5_log_dir)
invisible(lapply(required_directories, dir.create, recursive = TRUE, showWarnings = FALSE))

theme_global <- theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1)

group_levels <- c("Control", "HMPVmild", "HMPVsevere", "RSVmild", "RSVsevere")
virus_levels <- c("Control", "HMPV", "RSV")
severity_levels <- c("Control", "mild", "severe")

group_colors <- c(
  "Control" = "#808080",
  "HMPVmild" = "#B2DF8A",
  "HMPVsevere" = "#78B478",
  "RSVmild" = "#74ADD1",
  "RSVsevere" = "#3278C8"
)

virus_colors <- c("HMPV" = "#78B478", "RSV" = "#3278C8")

neutrophil_levels <- c(
  "IL1R2+ Immature Neutrophils",
  "IL1R2+ Mature Neutrophils",
  "Immature Neutrophils",
  "Mature Neutrophils",
  "Transitional Neutrophils"
)

neutrophil_colors <- c(
  "IL1R2+ Immature Neutrophils" = "#DE9DB5",
  "IL1R2+ Mature Neutrophils" = "#85B293",
  "Immature Neutrophils" = "#E7EFC8",
  "Mature Neutrophils" = "#6FABDC",
  "Transitional Neutrophils" = "#AB9AC4"
)

save_plot_pdf <- function(plot_object, file, width, height) {
  ggsave(filename = file, plot = plot_object, width = width, height = height, units = "in", device = "pdf", limitsize = FALSE)
}

mean_gene_expression <- function(expression_matrix, genes, min_genes = 3) {
  genes_use <- intersect(genes, rownames(expression_matrix))
  if (length(genes_use) < min_genes) {
    return(rep(NA_real_, ncol(expression_matrix)))
  }
  as.numeric(Matrix::colMeans(expression_matrix[genes_use, , drop = FALSE]))
}

pseudobulk_program_score <- function(expression_matrix, genes, min_genes = 3) {
  genes_use <- intersect(genes, rownames(expression_matrix))
  if (length(genes_use) < min_genes) {
    return(rep(NA_real_, ncol(expression_matrix)))
  }

  expression_subset <- expression_matrix[genes_use, , drop = FALSE]
  gene_sd <- apply(expression_subset, 1, sd, na.rm = TRUE)
  expression_subset <- expression_subset[is.finite(gene_sd) & gene_sd > 0, , drop = FALSE]

  if (nrow(expression_subset) < min_genes) {
    return(rep(NA_real_, ncol(expression_matrix)))
  }

  z_matrix <- t(scale(t(expression_subset)))
  colMeans(z_matrix, na.rm = TRUE)
}



# 1. Load Figure 4 object -------------------------------------------------
input_object_file <- file.path(fig4_object_dir, "Fig4_03_major_celltype_harmony.rds")

if (!file.exists(input_object_file)) {
  stop("The Figure 4 object was not found: ", input_object_file)
}

pbmc <- readRDS(input_object_file)

if (!inherits(pbmc, "Seurat")) {
  stop("The input object is not a Seurat object.")
}

required_metadata <- c("sample_id", "group", "virus", "severity", "major_celltype")
missing_metadata <- setdiff(required_metadata, colnames(pbmc@meta.data))

if (length(missing_metadata) > 0) {
  stop("Missing metadata columns: ", paste(missing_metadata, collapse = ", "))
}

DefaultAssay(pbmc) <- "RNA"



# 2. Extract neutrophils -------------------------------------------------
neut <- subset(pbmc, subset = major_celltype == "Neutrophil")

if (ncol(neut) == 0) {
  stop("No neutrophil cells were found.")
}

neut$group <- factor(as.character(neut$group), levels = group_levels)
neut$virus <- factor(as.character(neut$virus), levels = virus_levels)
neut$severity <- factor(as.character(neut$severity), levels = severity_levels)

saveRDS(neut, file.path(fig5_object_dir, "Fig5_01_neutrophil_subset.rds"))

neut_sample_summary <- neut@meta.data %>%
  as.data.frame() %>%
  count(sample_id, group, virus, severity, name = "neutrophil_cells") %>%
  arrange(group, sample_id)

write.csv(neut_sample_summary, file.path(fig5_table_dir, "Fig5_neutrophil_cells_by_sample.csv"), row.names = FALSE)



# 3. Normalize and integrate neutrophils --------------------------------
DefaultAssay(neut) <- "RNA"

if (length(Layers(neut[["RNA"]])) > 1) {
  neut[["RNA"]] <- JoinLayers(neut[["RNA"]])
}

neut[["RNA"]] <- split(neut[["RNA"]], f = neut$sample_id)

neut <- NormalizeData(neut, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
neut <- FindVariableFeatures(neut, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
neut <- ScaleData(neut, features = VariableFeatures(neut), verbose = FALSE)
neut <- RunPCA(neut, features = VariableFeatures(neut), npcs = 30, verbose = FALSE)

p_elbow <- ElbowPlot(neut, ndims = 30) + ggtitle("Neutrophil PCA elbow plot")
save_plot_pdf(p_elbow, file.path(fig5_figure_dir, "Fig5_neutrophil_ElbowPlot.pdf"), width = 6, height = 5)

dims <- 1:30

neut <- IntegrateLayers(
  object = neut,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = TRUE
)

if (!"harmony" %in% Reductions(neut)) {
  stop("Harmony integration failed.")
}

neut <- FindNeighbors(neut, reduction = "harmony", dims = dims, graph.name = c("neut_harmony_nn", "neut_harmony_snn"), verbose = FALSE)

neut <- RunUMAP(
  neut,
  reduction = "harmony",
  dims = dims,
  reduction.name = "umap.neut",
  reduction.key = "neutUMAP_",
  seed.use = 1234,
  verbose = FALSE
)



# 4. UMAP by sample and group --------------------------------------------
p_neut_sample <- DimPlot(neut, reduction = "umap.neut", group.by = "sample_id", raster = FALSE, pt.size = 0.8) +
  theme_global +
  ggtitle("Neutrophils by sample")

p_neut_group <- DimPlot(neut, reduction = "umap.neut", group.by = "group", raster = FALSE, pt.size = 0.8, cols = group_colors) +
  theme_global +
  ggtitle("Neutrophils by group")

save_plot_pdf(p_neut_sample, file.path(fig5_figure_dir, "Fig5_neutrophil_Harmony_UMAP_by_sample.pdf"), width = 8, height = 6)
save_plot_pdf(p_neut_group, file.path(fig5_figure_dir, "Fig5_neutrophil_Harmony_UMAP_by_group.pdf"), width = 8, height = 6)



# 5. Multi-resolution clustering ----------------------------------------
resolution_list <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)

neut <- FindClusters(
  neut,
  graph.name = "neut_harmony_snn",
  resolution = resolution_list,
  algorithm = 1,
  random.seed = 1234,
  verbose = FALSE
)

resolution_columns <- grep("^neut_harmony_snn_res\\.", colnames(neut@meta.data), value = TRUE)

resolution_summary <- data.frame(
  resolution_column = resolution_columns,
  cluster_number = vapply(resolution_columns, function(x) length(unique(neut@meta.data[[x]])), integer(1))
) %>%
  mutate(distance_to_9 = abs(cluster_number - 9)) %>%
  arrange(distance_to_9, resolution_column)

write.csv(resolution_summary, file.path(fig5_table_dir, "Fig5_clustering_resolution_summary.csv"), row.names = FALSE)

selected_cluster_column <- resolution_summary$resolution_column[1]

if (is.na(selected_cluster_column)) {
  stop("No clustering resolution was selected.")
}

neut$neut_cluster <- factor(
  as.character(neut@meta.data[[selected_cluster_column]]),
  levels = sort(unique(as.character(neut@meta.data[[selected_cluster_column]])))
)

Idents(neut) <- "neut_cluster"

p_resolution <- lapply(resolution_columns, function(cluster_column) {
  DimPlot(neut, reduction = "umap.neut", group.by = cluster_column, label = TRUE, repel = TRUE, raster = TRUE) +
    theme_global +
    ggtitle(cluster_column)
})

p_resolution_comparison <- wrap_plots(p_resolution, ncol = 2)
save_plot_pdf(p_resolution_comparison, file.path(fig5_figure_dir, "Fig5_neutrophil_clustering_resolution_comparison.pdf"), width = 12, height = 12)

p_cluster <- DimPlot(neut, reduction = "umap.neut", group.by = "neut_cluster", label = TRUE, repel = TRUE, raster = FALSE, pt.size = 0.8) +
  theme_global + theme(legend.title = element_blank()) + ggtitle("Neutrophil clusters")

save_plot_pdf(p_cluster, file.path(fig5_figure_dir, "Fig5A_neutrophil_cluster_UMAP.pdf"), width = 8, height = 6)



# 6. Cluster marker analysis ---------------------------------------------
neut_marker <- neut
DefaultAssay(neut_marker) <- "RNA"

if (length(Layers(neut_marker[["RNA"]])) > 1) {
  neut_marker[["RNA"]] <- JoinLayers(neut_marker[["RNA"]])
}

Idents(neut_marker) <- "neut_cluster"

neut_cluster_markers <- FindAllMarkers(
  neut_marker,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.3,
  logfc.threshold = 0.25
) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  ungroup()

write.csv(neut_cluster_markers, file.path(fig5_table_dir, "Fig5_neutrophil_all_cluster_markers.csv"), row.names = FALSE)

neut_cluster_top30 <- neut_cluster_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 30, with_ties = FALSE) %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  ungroup()

write.csv(neut_cluster_top30, file.path(fig5_table_dir, "Fig5_neutrophil_cluster_top30_markers.csv"), row.names = FALSE)

# 7. Neutrophil marker plots ---------------------------------------------

neut_feature_markers <- c(
  "CSF3R", "FCGR3B", "CEACAM8",
  "MPO", "ELANE", "AZU1", "OLFM4", "CAMP", "LTF", "KIT", "ITGA4", "CXCR4",
  "S100A8", "S100A9", "MMP8", "MMP9", "ITGAM",
  "MME", "CD101", "CXCR2", "FPR1", "SELL", "IL1R2"
)

neut_feature_markers <- intersect(neut_feature_markers, rownames(neut))

p_neut_dotplot <- DotPlot(neut, features = neut_feature_markers, group.by = "neut_cluster") +
  RotatedAxis() +
  scale_color_distiller(palette = "RdYlBu") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

save_plot_pdf(p_neut_dotplot, file.path(fig5_figure_dir, "Fig5_neutrophil_marker_DotPlot.pdf"), width = 10, height = 5)

neut_plot <- neut
DefaultAssay(neut_plot) <- "RNA"

if (length(Layers(neut_plot[["RNA"]])) > 1) {
  neut_plot[["RNA"]] <- JoinLayers(neut_plot[["RNA"]])
}

pdf(file.path(fig5_figure_dir, "Fig5B_neutrophil_marker_FeaturePlot.pdf"), width = 6, height = 6)
for (gene in neut_feature_markers) {
  plot_result <- tryCatch(
    {
      FeaturePlot(neut_plot, features = gene, reduction = "umap.neut", raster = TRUE, order = TRUE) +
        theme_global +
        scale_color_distiller(palette = "RdYlBu") +
        theme(
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          plot.title = element_text(hjust = 0.5)
        )
    },
    error = function(e) NULL
  )

  if (!is.null(plot_result)) {
    print(plot_result)
  }
}
dev.off()



# 8. Manual neutrophil-state annotation ---------------------------------
if (!all(as.character(0:8) %in% levels(neut$neut_cluster))) {
  stop("The expected nine neutrophil clusters were not found.")
}

cluster_neutrophil <- c(
  "0" = "Mature Neutrophils",
  "1" = "Transitional Neutrophils",
  "2" = "IL1R2+ Mature Neutrophils",
  "3" = "Transitional Neutrophils",
  "4" = "IL1R2+ Mature Neutrophils",
  "5" = "IL1R2+ Immature Neutrophils",
  "6" = "Transitional Neutrophils",
  "7" = "Transitional Neutrophils",
  "8" = "Immature Neutrophils"
)

neut$neutrophil <- unname(cluster_neutrophil[as.character(neut$neut_cluster)])
neut$neutrophil <- factor(neut$neutrophil, levels = neutrophil_levels)

neut$group <- factor(neut$group, levels = group_levels)
neut$virus <- factor(neut$virus, levels = virus_levels)
neut$severity <- factor(neut$severity, levels = severity_levels)

if (anyNA(neut$neutrophil)) {
  stop("Some cells have missing neutrophil-state annotations.")
}

saveRDS(neut, file.path(fig5_object_dir, "Fig5_03_neutrophil_annotation_final.rds"))



# 9. Neutrophil subtype composition --------------------------------------
cell_prop_neut_group <- neut@meta.data %>%
  as.data.frame() %>%
  count(group, neutrophil, name = "cell_number") %>%
  complete(group = group_levels, neutrophil = neutrophil_levels, fill = list(cell_number = 0)) %>%
  group_by(group) %>%
  mutate(
    total_cells = sum(cell_number),
    proportion = cell_number / total_cells,
    percentage = proportion * 100
  ) %>%
  ungroup()

write.csv(cell_prop_neut_group, file.path(fig5_table_dir, "Fig5D_neutrophil_proportion_by_group.csv"), row.names = FALSE)

p_fig5D <- ggplot(cell_prop_neut_group, aes(x = group, y = proportion, fill = neutrophil)) +
  geom_col(width = 0.85, color = "white", linewidth = 0) +
  scale_fill_manual(values = neutrophil_colors, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), breaks = seq(0, 1, by = 0.2), expand = expansion(mult = c(0.01, 0.04))) +
  labs(x = NULL, y = "Cell proportion", fill = NULL) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )

save_plot_pdf(p_fig5D, file.path(fig5_figure_dir, "Fig5D_neutrophil_proportion_by_group.pdf"), width = 6, height = 5.5)

p_neut_umap <- DimPlot(neut, reduction = "umap.neut", group.by = "neutrophil", raster = FALSE, label = TRUE, repel = TRUE, pt.size = 0.8, cols = neutrophil_colors) +
  theme_global +
  theme(legend.title = element_blank()) +
  ggtitle("Neutrophil states")

save_plot_pdf(p_neut_umap, file.path(fig5_figure_dir, "Fig5C_neutrophil_celltype_UMAP.pdf"), width = 7, height = 6)

cell_prop_neut_sample <- neut@meta.data %>%
  as.data.frame() %>%
  count(sample_id, group, neutrophil, name = "cell_number") %>%
  complete(nesting(sample_id, group), neutrophil = neutrophil_levels, fill = list(cell_number = 0)) %>%
  group_by(sample_id, group) %>%
  mutate(
    total_neutrophil = sum(cell_number),
    proportion = cell_number / total_neutrophil,
    percentage = proportion * 100
  ) %>%
  ungroup()

write.csv(cell_prop_neut_sample, file.path(fig5_table_dir, "Fig5D_neutrophil_proportion_by_sample.csv"), row.names = FALSE)



# 10. HMPV versus RSV differential expression ----------------------------
neut_de <- subset(neut, subset = virus %in% c("HMPV", "RSV"))

if (length(Layers(neut_de[["RNA"]])) > 1) {
  neut_de[["RNA"]] <- JoinLayers(neut_de[["RNA"]])
}

Idents(neut_de) <- "virus"

deg_HMPV_RSV <- FindMarkers(
  neut_de,
  ident.1 = "HMPV",
  ident.2 = "RSV",
  assay = "RNA",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.2,
  only.pos = FALSE,
  verbose = FALSE
)

deg_HMPV_RSV$gene <- rownames(deg_HMPV_RSV)

if (!"avg_log2FC" %in% colnames(deg_HMPV_RSV) && "avg_logFC" %in% colnames(deg_HMPV_RSV)) {
  deg_HMPV_RSV$avg_log2FC <- deg_HMPV_RSV$avg_logFC
}

if (!"avg_log2FC" %in% colnames(deg_HMPV_RSV)) {
  stop("No log-fold-change column was found.")
}

deg_HMPV_RSV <- deg_HMPV_RSV %>%
  mutate(
    regulation = case_when(
      p_val_adj < 0.0001 & avg_log2FC >= 0.25 ~ "Up",
      p_val_adj < 0.0001 & avg_log2FC <= -0.25 ~ "Down",
      TRUE ~ "NS"
    )
  ) %>%
  arrange(desc(avg_log2FC))

write.csv(deg_HMPV_RSV, file.path(fig5_table_dir, "Fig5E_DEG_HMPV_vs_RSV_all_neutrophils.csv"), row.names = FALSE)
write.csv(filter(deg_HMPV_RSV, regulation != "NS"), file.path(fig5_table_dir, "Fig5E_DEG_HMPV_vs_RSV_significant.csv"), row.names = FALSE)

HMPV_up_genes <- filter(deg_HMPV_RSV, regulation == "Up")
RSV_up_genes <- filter(deg_HMPV_RSV, regulation == "Down")

write.csv(HMPV_up_genes, file.path(fig5_table_dir, "Fig5E_HMPV_up_genes.csv"), row.names = FALSE)
write.csv(RSV_up_genes, file.path(fig5_table_dir, "Fig5E_RSV_up_genes.csv"), row.names = FALSE)



# 11. KEGG enrichment -----------------------------------------------------
up_gene_symbol <- unique(filter(deg_HMPV_RSV, regulation == "Up")$gene)
down_gene_symbol <- unique(filter(deg_HMPV_RSV, regulation == "Down")$gene)

up_gene_map <- bitr(up_gene_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
down_gene_map <- bitr(down_gene_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
background_gene_map <- bitr(rownames(neut_de), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

up_entrez <- unique(up_gene_map$ENTREZID)
down_entrez <- unique(down_gene_map$ENTREZID)
background_entrez <- unique(background_gene_map$ENTREZID)

kegg_up <- enrichKEGG(gene = up_entrez, organism = "hsa", universe = background_entrez, pvalueCutoff = 0.05, qvalueCutoff = 0.2, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500)
kegg_down <- enrichKEGG(gene = down_entrez, organism = "hsa", universe = background_entrez, pvalueCutoff = 0.05, qvalueCutoff = 0.2, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500)

kegg_up_df <- as.data.frame(kegg_up)
kegg_down_df <- as.data.frame(kegg_down)

if (nrow(kegg_up_df) > 0) kegg_up_df$Direction <- "Up"
if (nrow(kegg_down_df) > 0) kegg_down_df$Direction <- "Down"

write.csv(kegg_up_df, file.path(fig5_table_dir, "Fig5E_KEGG_HMPV_up.csv"), row.names = FALSE)
write.csv(kegg_down_df, file.path(fig5_table_dir, "Fig5E_KEGG_HMPV_down.csv"), row.names = FALSE)

kegg_all <- bind_rows(kegg_up_df, kegg_down_df) %>%
  filter(!is.na(p.adjust)) %>%
  mutate(
    FDR = p.adjust,
    signed_log10FDR = if_else(Direction == "Up", -log10(pmax(FDR, 1e-50)), log10(pmax(FDR, 1e-50)))
  )

write.csv(kegg_all, file.path(fig5_table_dir, "Fig5E_KEGG_HMPV_vs_RSV_all.csv"), row.names = FALSE)



# 12. Selected KEGG pathways ---------------------------------------------
pathway_select <- tibble(
  Description = c(
    "RIG-I-like receptor signaling pathway",
    "NOD-like receptor signaling pathway",
    "Toll-like receptor signaling pathway",
    "C-type lectin receptor signaling pathway",
    "TNF signaling pathway",
    "IL-17 signaling pathway",
    "Chemokine signaling pathway",
    "MAPK signaling pathway",
    "Fc gamma R-mediated phagocytosis",
    "Phagosome",
    "Platelet activation",
    "Leukocyte transendothelial migration"
  ),
  pathway_class = c(
    "Antiviral sensing",
    "Innate sensing",
    "Innate sensing",
    "Innate sensing",
    "Inflammation",
    "Inflammation",
    "Chemotaxis",
    "Activation",
    "Phagocytosis",
    "Phagocytosis",
    "Activation",
    "Migration"
  ),
  pathway_order = seq_len(12)
)

plotdata <- kegg_all %>%
  inner_join(pathway_select, by = "Description") %>%
  group_by(Description, Direction) %>%
  slice_min(order_by = FDR, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(pathway_order) %>%
  mutate(
    Direction = factor(Direction, levels = c("Down", "Up")),
    Description_plot = factor(Description, levels = rev(unique(Description))),
    logFDR = -log10(pmax(FDR, 1e-50)),
    signed_logFDR = if_else(Direction == "Up", logFDR, -logFDR),
    label = paste0("FDR = ", format(FDR, scientific = TRUE, digits = 2), ", count = ", Count)
  )

write.csv(plotdata, file.path(fig5_table_dir, "Fig5E_KEGG_selected_plotdata.csv"), row.names = FALSE)

if (nrow(plotdata) > 0) {
  x_max <- ceiling(max(abs(plotdata$signed_logFDR), na.rm = TRUE) + 1)

  p_fig5E <- ggplot(plotdata, aes(y = Description_plot)) +
    geom_segment(aes(x = 0, xend = signed_logFDR, yend = Description_plot, color = Direction), linewidth = 1) +
    geom_point(aes(x = signed_logFDR), color = "black", fill = "white", size = 8, shape = 21, stroke = 0.7) +
    geom_vline(xintercept = 0, linetype = 4, linewidth = 0.6) +
    geom_text(data = filter(plotdata, Direction == "Down"), aes(x = 0.1, label = label), hjust = 0, size = 2.5, inherit.aes = TRUE) +
    geom_text(data = filter(plotdata, Direction == "Up"), aes(x = -0.1, label = label), hjust = 1, size = 2.5, inherit.aes = TRUE) +
    scale_color_manual(values = c("Down" = "#A0D8EA", "Up" = "#EFA39F")) +
    scale_x_continuous(limits = c(-x_max, x_max)) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 42)) +
    labs(x = expression("Signed " * -log[10] * "(FDR)"), y = NULL) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(color = "black"),
      legend.position = "none"
    )

  save_plot_pdf(p_fig5E, file.path(fig5_figure_dir, "Fig5E_KEGG_selected_butterfly.pdf"), width = 7, height = 5)
}



# 13. Save object and session information --------------------------------
saveRDS(neut, file.path(fig5_object_dir, "Fig5_05_neutrophil_analysis_final.rds"))
writeLines(capture.output(sessionInfo()), file.path(fig5_log_dir, "Fig5_sessionInfo.txt"))

message("Figure 5 neutrophil analysis completed successfully.")


# 14. Monocle2 pseudotime analysis ---------------------------------------
suppressPackageStartupMessages({
  library(monocle)
  library(Biobase)
})

set.seed(1234)

pseudotime_object_file <- file.path(fig5_object_dir, "Fig5_06_neutrophil_30percent_for_pseudotime.rds")

if ("neut_cluster" %in% colnames(neut@meta.data)) {
  sampling_cluster <- "neut_cluster"
} else {
  sampling_cluster <- "seurat_clusters"
}

sampling_df <- neut@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell")

selected_cells <- sampling_df %>%
  group_by(.data[[sampling_cluster]]) %>%
  slice_sample(prop = 0.3) %>%
  pull(cell)

write.csv(data.frame(cell = selected_cells), file.path(fig5_table_dir, "Fig5FG_selected_30percent_cells.csv"), row.names = FALSE)

pseudotime_data <- subset(neut, cells = selected_cells)
DefaultAssay(pseudotime_data) <- "RNA"

if (length(Layers(pseudotime_data[["RNA"]])) > 1) {
  pseudotime_data[["RNA"]] <- JoinLayers(pseudotime_data[["RNA"]])
}

pseudotime_data <- NormalizeData(
  pseudotime_data,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

pseudotime_data <- FindVariableFeatures(
  pseudotime_data,
  assay = "RNA",
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

saveRDS(pseudotime_data, pseudotime_object_file)

counts_matrix <- GetAssayData(
  pseudotime_data,
  assay = "RNA",
  layer = "counts"
)

gene_annotation <- data.frame(
  gene_short_name = rownames(counts_matrix),
  row.names = rownames(counts_matrix)
)

feature_data <- new("AnnotatedDataFrame", data = gene_annotation)
phenotype_data <- new("AnnotatedDataFrame", data = pseudotime_data@meta.data)

monocle_cds <- newCellDataSet(
  counts_matrix,
  phenoData = phenotype_data,
  featureData = feature_data,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit = 1
)

monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- detectGenes(monocle_cds)

ordering_genes <- VariableFeatures(pseudotime_data)
ordering_genes <- intersect(ordering_genes, rownames(monocle_cds))

ordering_genes <- ordering_genes[
  fData(monocle_cds)[ordering_genes, "num_cells_expressed"] >= 10
]

write.csv(data.frame(gene = ordering_genes), file.path(fig5_table_dir, "Fig5FG_Monocle2_ordering_genes.csv"), row.names = FALSE)

monocle_cds <- setOrderingFilter(monocle_cds, ordering_genes)

set.seed(1234)

monocle_cds <- reduceDimension(
  monocle_cds,
  max_components = 2,
  reduction_method = "DDRTree",
  norm_method = "log",
  verbose = TRUE
)

monocle_cds <- orderCells(monocle_cds)

available_states <- sort(unique(as.numeric(as.character(pData(monocle_cds)$State))))

if (1 %in% available_states) {
  monocle_cds <- orderCells(monocle_cds, root_state = 1)
} else {
  warning("State 1 was not found. The original Monocle2 root state was retained.")
}

saveRDS(monocle_cds, file.path(fig5_object_dir, "Fig5_09_neutrophil_Monocle2_DDRTree_final.rds"))

trajectory_metadata <- pData(monocle_cds) %>%
  as.data.frame() %>%
  rownames_to_column("cell")

trajectory_coordinates <- t(reducedDimS(monocle_cds)) %>% as.data.frame()

colnames(trajectory_coordinates)[1:2] <- c("Component1", "Component2")

trajectory_coordinates <- trajectory_coordinates %>% rownames_to_column("cell")

trajectory_metadata <- trajectory_metadata %>% left_join(trajectory_coordinates, by = "cell")

write.csv(trajectory_metadata, file.path(fig5_table_dir, "Fig5FG_neutrophil_pseudotime_metadata.csv"), row.names = FALSE)

state_by_neutrophil <- trajectory_metadata %>%
  count(State, neutrophil, name = "cells") %>%
  group_by(State) %>%
  mutate(proportion = cells / sum(cells)) %>%
  ungroup()

state_by_group <- trajectory_metadata %>%
  count(State, group, name = "cells") %>%
  group_by(State) %>%
  mutate(proportion = cells / sum(cells)) %>%
  ungroup()

write.csv(state_by_neutrophil, file.path(fig5_table_dir, "Fig5FG_State_by_neutrophil.csv"), row.names = FALSE)
write.csv(state_by_group, file.path(fig5_table_dir, "Fig5FG_State_by_group.csv"), row.names = FALSE)

sample_pseudotime_summary <- trajectory_metadata %>%
  group_by(sample_id, group) %>%
  summarise(
    cells = n(),
    median_pseudotime = median(Pseudotime, na.rm = TRUE),
    mean_pseudotime = mean(Pseudotime, na.rm = TRUE),
    Q1 = quantile(Pseudotime, 0.25, na.rm = TRUE),
    Q3 = quantile(Pseudotime, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(sample_pseudotime_summary, file.path(fig5_table_dir, "Fig5FG_sample_pseudotime_summary.csv"), row.names = FALSE)

plot_state <- plot_cell_trajectory(
  monocle_cds,
  color_by = "State",
  cell_size = 0.7,
  show_branch_points = TRUE
) +
  theme_classic()

save_plot_pdf(
  plot_state,
  file.path(fig5_figure_dir, "Fig5F_neutrophil_trajectory_State.pdf"),
  width = 7,
  height = 5.5
)

plot_pseudotime <- plot_cell_trajectory(
  monocle_cds,
  color_by = "Pseudotime",
  cell_size = 0.7,
  show_branch_points = FALSE
) +
  scale_color_gradientn(colors = c("#102F46", "#1F5576", "#347CAC", "#63ACE1")) +
  theme_classic()

save_plot_pdf(plot_pseudotime, file.path(fig5_figure_dir, "Fig5G_neutrophil_trajectory_Pseudotime.pdf"), width = 7, height = 5.5)

plot_group <- plot_cell_trajectory(
  monocle_cds,
  color_by = "group",
  cell_size = 0.7,
  show_branch_points = FALSE
) +
  scale_color_manual(values = group_colors) +
  theme_classic()

save_plot_pdf(plot_group, file.path(fig5_figure_dir, "Fig5G_neutrophil_trajectory_Group.pdf"), width = 7, height = 5.5)



# 15. Reviewer-requested IL1R2 and interferon analysis -------------------
ISG_core <- c(
  "ISG15", "IFIT1", "IFIT2", "IFIT3",
  "IFITM1", "IFITM2", "IFITM3",
  "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL",
  "RSAD2", "IFI6", "IFI27", "IFI35", "IFI44", "IFI44L",
  "USP18", "HERC5", "BST2", "XAF1", "DDX58", "IFIH1", "IRF7", "STAT1"
)

IFNA_genes <- ISG_core

IFNG_genes <- c(
  "STAT1", "IRF1", "GBP1", "GBP2", "GBP4", "GBP5",
  "CXCL9", "CXCL10", "TAP1", "TAP2",
  "PSMB8", "PSMB9", "HLA-B", "HLA-C",
  "ISG15", "IFIT1", "IFIT3", "MX1", "OAS1", "BST2"
)

Inflammatory_genes <- c(
  "IL1B", "TNF", "NFKBIA", "NFKBIZ", "CXCL8", "CXCL2",
  "CCL3", "CCL4", "PTGS2", "FOS", "JUN", "JUNB",
  "S100A8", "S100A9", "S100A12",
  "ICAM1", "PLAUR", "C5AR1", "FCGR3B", "TREM1"
)

IFN_combined_genes <- unique(c(IFNA_genes, IFNG_genes))

DefaultAssay(neut) <- "RNA"

rna_data <- tryCatch(
  GetAssayData(neut, assay = "RNA", layer = "data"),
  error = function(e) NULL
)

if (is.null(rna_data)) {
  neut[["RNA"]] <- JoinLayers(neut[["RNA"]])
  rna_data <- GetAssayData(neut, assay = "RNA", layer = "data")
}

IFNA_present <- intersect(IFNA_genes, rownames(rna_data))
IFNG_present <- intersect(IFNG_genes, rownames(rna_data))
IFN_combined_present <- intersect(IFN_combined_genes, rownames(rna_data))
Inflammatory_present <- intersect(Inflammatory_genes, rownames(rna_data))
ISG_core_present <- intersect(ISG_core, rownames(rna_data))

if (!"IL1R2" %in% rownames(rna_data)) {
  stop("IL1R2 was not detected in the RNA assay.")
}

neut$IFNA_score <- mean_gene_expression(rna_data, IFNA_present)
neut$IFNG_score <- mean_gene_expression(rna_data, IFNG_present)
neut$IFN_combined_score <- mean_gene_expression(rna_data, IFN_combined_present)
neut$Inflammatory_score <- mean_gene_expression(rna_data, Inflammatory_present)
neut$IL1R2_expression <- as.numeric(rna_data["IL1R2", ])

ISG_detected_number <- Matrix::colSums(
  rna_data[ISG_core_present, , drop = FALSE] > 0
)

neut$ISG_detected_number <- as.numeric(ISG_detected_number)
neut$IL1R2_positive <- neut$IL1R2_expression > 0
neut$ISG_positive <- neut$ISG_detected_number >= 3

neut$IL1R2_ISG_status <- case_when(
  neut$IL1R2_positive & neut$ISG_positive ~ "IL1R2+ / ISG+",
  neut$IL1R2_positive & !neut$ISG_positive ~ "IL1R2+ / ISG-",
  !neut$IL1R2_positive & neut$ISG_positive ~ "IL1R2- / ISG+",
  TRUE ~ "IL1R2- / ISG-"
)

neut$IL1R2_ISG_status <- factor(
  neut$IL1R2_ISG_status,
  levels = c(
    "IL1R2- / ISG-",
    "IL1R2+ / ISG-",
    "IL1R2- / ISG+",
    "IL1R2+ / ISG+"
  )
)

cell_metadata <- neut@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell")

gene_set_check <- bind_rows(
  data.frame(programme = "IFNA", gene = IFNA_genes),
  data.frame(programme = "IFNG", gene = IFNG_genes),
  data.frame(programme = "Inflammatory", gene = Inflammatory_genes),
  data.frame(programme = "ISG_core", gene = ISG_core)
) %>%
  mutate(detected = gene %in% rownames(rna_data))

write.csv(gene_set_check, file.path(fig5_table_dir, "Reviewer15_gene_set_detection.csv"), row.names = FALSE)



# 15.1 Sample-level neutrophil subtype composition -----------------------
sample_info <- cell_metadata %>%
  distinct(sample_id, group, virus, severity)

sample_subtype_composition <- cell_metadata %>%
  count(sample_id, neutrophil, name = "cells") %>%
  complete(
    nesting(sample_id),
    neutrophil = neutrophil_levels,
    fill = list(cells = 0)
  ) %>%
  left_join(sample_info, by = "sample_id") %>%
  group_by(sample_id) %>%
  mutate(
    total_neutrophils = sum(cells),
    proportion = cells / total_neutrophils,
    percent = proportion * 100
  ) %>%
  ungroup()

sample_subtype_composition$group <- factor(
  sample_subtype_composition$group,
  levels = group_levels
)

sample_subtype_composition$neutrophil <- factor(
  sample_subtype_composition$neutrophil,
  levels = neutrophil_levels
)

write.csv(sample_subtype_composition, file.path(fig5_table_dir, "Reviewer15_sample_level_neutrophil_composition.csv"), row.names = FALSE)

p_reviewer15_composition <- ggplot(
  sample_subtype_composition,
  aes(x = group, y = percent, colour = group)
) +
  geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.4, fill = NA) +
  geom_jitter(width = 0.12, size = 2.5) +
  facet_wrap(~neutrophil, ncol = 3, scales = "free_y") +
  scale_colour_manual(values = group_colors) +
  theme_global +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(x = NULL, y = "Neutrophils (%)")

save_plot_pdf(p_reviewer15_composition, file.path(fig5_figure_dir, "Reviewer15_neutrophil_subtypes_sample_level.pdf"), width = 8, height = 6)



# 15.2 IL1R2 and ISG co-expression by sample -----------------------------
sample_coexpression <- cell_metadata %>%
  group_by(sample_id, group, virus, severity) %>%
  summarise(
    total_neutrophils = n(),
    IL1R2_positive_cells = sum(IL1R2_positive, na.rm = TRUE),
    ISG_positive_cells = sum(ISG_positive, na.rm = TRUE),
    double_positive_cells = sum(IL1R2_positive & ISG_positive, na.rm = TRUE),
    IL1R2_positive_percent = 100 * mean(IL1R2_positive, na.rm = TRUE),
    ISG_positive_percent = 100 * mean(ISG_positive, na.rm = TRUE),
    double_positive_percent = 100 * mean(IL1R2_positive & ISG_positive, na.rm = TRUE),
    mean_IL1R2 = mean(IL1R2_expression, na.rm = TRUE),
    mean_IFN_score = mean(IFN_combined_score, na.rm = TRUE),
    mean_inflammatory_score = mean(Inflammatory_score, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(sample_coexpression, file.path(fig5_table_dir, "Reviewer15_IL1R2_ISG_coexpression_by_sample.csv"), row.names = FALSE)

coexpression_colors <- c(
  "IL1R2- / ISG-" = "#D9D9D9",
  "IL1R2+ / ISG-" = "#F5B47F",
  "IL1R2- / ISG+" = "#4581B2",
  "IL1R2+ / ISG+" = "#CB5C5B"
)

umap_coordinates <- Embeddings(neut, reduction = "umap.neut")[, 1:2, drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("cell")

colnames(umap_coordinates)[2:3] <- c("UMAP_1", "UMAP_2")

umap_plot_data <- umap_coordinates %>%
  left_join(
    cell_metadata %>%
      select(cell, IL1R2_expression, IFN_combined_score, IL1R2_ISG_status),
    by = "cell"
  ) %>%
  arrange(IL1R2_ISG_status)

p_il1r2_umap <- ggplot(umap_plot_data, aes(UMAP_1, UMAP_2, colour = IL1R2_expression)) +
  geom_point(size = 0.6) +
  scale_colour_gradient(low = "#F2F2F2", high = "#D73027") +
  theme_void() + labs(colour = "IL1R2")

p_ifn_umap <- ggplot(umap_plot_data, aes(UMAP_1, UMAP_2, colour = IFN_combined_score)) +
  geom_point(size = 0.6) +
  scale_colour_gradient(low = "#F2F2F2", high = "#2166AC") +
  theme_void() + labs(colour = "IFN score")

p_coexpression_umap <- ggplot(umap_plot_data, aes(UMAP_1, UMAP_2, colour = IL1R2_ISG_status)) +
  geom_point(size = 0.6) +
  scale_colour_manual(values = coexpression_colors, drop = FALSE) +
  theme_void() + labs(colour = NULL)

p_reviewer15_umap <- p_il1r2_umap + p_ifn_umap + p_coexpression_umap +
  plot_annotation(tag_levels = "A")

save_plot_pdf(p_reviewer15_umap, file.path(fig5_figure_dir, "Reviewer15_IL1R2_ISG_coexpression_UMAP.pdf"), width = 13, height = 4)

p_double_positive <- ggplot(
  sample_coexpression,
  aes(x = group, y = double_positive_percent, colour = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, fill = NA, linewidth = 0.5) +
  geom_jitter(width = 0.12, size = 3) +
  scale_colour_manual(values = group_colors) +
  theme_global +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(x = NULL, y = "IL1R2+ / ISG+ neutrophils (%)")

save_plot_pdf(p_double_positive, file.path(fig5_figure_dir, "Reviewer15_IL1R2_ISG_double_positive_by_group.pdf"), width = 5, height = 5)



# 15.3 Sample-level pseudobulk programme scores --------------------------
counts_data <- tryCatch(
  GetAssayData(neut, assay = "RNA", layer = "counts"),
  error = function(e) NULL
)

if (is.null(counts_data)) {
  neut[["RNA"]] <- JoinLayers(neut[["RNA"]])
  counts_data <- GetAssayData(neut, assay = "RNA", layer = "counts")
}

sample_factor <- factor(as.character(neut$sample_id))
sample_design <- Matrix::sparse.model.matrix(~0 + sample_factor)
colnames(sample_design) <- levels(sample_factor)

pseudobulk_counts <- counts_data %*% sample_design
pseudobulk_counts_dense <- as.matrix(pseudobulk_counts)

library_size <- colSums(pseudobulk_counts_dense)

pseudobulk_logCPM <- log2(sweep(pseudobulk_counts_dense, 2, library_size, "/") * 1e6 + 1)

saveRDS(pseudobulk_counts, file.path(fig5_object_dir, "Reviewer15_neutrophil_pseudobulk_counts.rds"))

saveRDS(pseudobulk_logCPM, file.path(fig5_object_dir, "Reviewer15_neutrophil_pseudobulk_logCPM.rds"))

sample_program_scores <- data.frame(
  sample_id = colnames(pseudobulk_logCPM),
  IFNA_score = pseudobulk_program_score(pseudobulk_logCPM, IFNA_genes),
  IFNG_score = pseudobulk_program_score(pseudobulk_logCPM, IFNG_genes),
  IFN_combined_score = pseudobulk_program_score(pseudobulk_logCPM, IFN_combined_genes),
  Inflammatory_score = pseudobulk_program_score(pseudobulk_logCPM, Inflammatory_genes),
  IL1R2_logCPM = as.numeric(pseudobulk_logCPM["IL1R2", ]),
  stringsAsFactors = FALSE
) %>%
  left_join(sample_info, by = "sample_id") %>%
  left_join(
    sample_coexpression %>%
      select(sample_id, double_positive_percent),
    by = "sample_id"
  )

write.csv(sample_program_scores, file.path(fig5_table_dir, "Reviewer15_neutrophil_sample_pseudobulk_program_scores.csv"), row.names = FALSE)



# 15.4 Programme scores by clinical group -------------------------------
programme_long <- sample_program_scores %>%
  select(
    sample_id,
    group,
    virus,
    severity,
    IFNA_score,
    IFNG_score,
    IFN_combined_score,
    Inflammatory_score,
    IL1R2_logCPM
  ) %>%
  pivot_longer(
    cols = c(IFNA_score, IFNG_score, IFN_combined_score, Inflammatory_score, IL1R2_logCPM),
    names_to = "programme",
    values_to = "score"
  )

programme_labels <- c(
  IFNA_score = "IFN-alpha programme",
  IFNG_score = "IFN-gamma programme",
  IFN_combined_score = "Combined IFN programme",
  Inflammatory_score = "Inflammatory programme",
  IL1R2_logCPM = "IL1R2 expression"
)

programme_long$programme <- factor(
  programme_labels[programme_long$programme],
  levels = unname(programme_labels)
)

p_programme_group <- ggplot(
  programme_long,
  aes(x = group, y = score, colour = group)) +
  geom_boxplot(outlier.shape = NA, fill = NA, width = 0.55, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 2.5) +
  facet_wrap(~programme, ncol = 3, scales = "free_y") +
  scale_colour_manual(values = group_colors) +
  theme_global +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(x = NULL, y = "Sample-level score")

save_plot_pdf(p_programme_group, file.path(fig5_figure_dir, "Reviewer15_IFN_inflammatory_IL1R2_scores_by_group.pdf"), width = 8, height = 6)



# 15.5 Save reviewer-requested results -----------------------------------
reviewer15_metadata <- neut@meta.data %>%
  as.data.frame() %>%
  select(
    sample_id,
    group,
    virus,
    severity,
    neutrophil,
    IL1R2_expression,
    ISG_detected_number,
    IL1R2_positive,
    ISG_positive,
    IL1R2_ISG_status,
    IFNA_score,
    IFNG_score,
    IFN_combined_score,
    Inflammatory_score
  )

saveRDS(neut, file.path(fig5_object_dir, "Fig5_neutrophil_Reviewer15_scores.rds"))
saveRDS(reviewer15_metadata, file.path(fig5_object_dir, "Reviewer15_cell_level_metadata.rds"))

writeLines(capture.output(sessionInfo()), file.path(fig5_log_dir, "Reviewer15_sessionInfo.txt"))

message("Figure 5 analysis and reviewer-requested analysis completed successfully.")

