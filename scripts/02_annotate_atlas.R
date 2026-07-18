#!/usr/bin/env Rscript
source("config/config.R")
source("R/helpers.R")
assert_packages(c("Seurat", "dplyr", "readr", "ggplot2"))

pnd <- readRDS(path("results", "objects", "01_pnd_clustered.rds"))
pnd <- tryCatch(Seurat::JoinLayers(pnd), error = function(e) pnd)

markers <- Seurat::FindAllMarkers(
  pnd,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
write_table(markers, "cluster_markers.csv")

brain_markers <- c(
  "P2ry12", "Tmem119", "Cx3cr1", "Mog", "Mbp", "Opalin", "Pdgfra", "Cspg4",
  "Aldh1l1", "Aqp4", "Slc1a2", "Cldn5", "Pecam1", "Rgs5", "Pdgfrb",
  "Ttr", "Kcne2", "Foxj1", "Dnah12", "Lyve1", "Cd163", "Ly6g", "S100a8",
  "Cd3d", "Cd3e", "Cd79a", "Cd19", "Gad1", "Gad2", "Slc17a7", "Snap25"
)
brain_markers <- safe_genes(brain_markers, pnd)

p_dot <- Seurat::DotPlot(pnd, features = brain_markers, group.by = "seurat_clusters") +
  Seurat::RotatedAxis() +
  ggplot2::theme_bw(base_size = 10)
save_plot(p_dot, "Figure2A_marker_dotplot.pdf", 16, 8)

# Manual annotation must be checked against markers before the repository is frozen.
cluster_ids <- levels(Seurat::Idents(pnd))
if (!all(cluster_ids %in% names(CLUSTER_LABELS))) {
  stop("CLUSTER_LABELS does not cover all Seurat clusters: ", paste(setdiff(cluster_ids, names(CLUSTER_LABELS)), collapse = ", "))
}
pnd <- Seurat::RenameIdents(pnd, CLUSTER_LABELS[cluster_ids])
pnd$celltype <- as.character(Seurat::Idents(pnd))

feature_genes <- safe_genes(c("Cx3cr1", "Gfap", "Mog", "Pdgfra"), pnd)
p_feature <- Seurat::FeaturePlot(pnd, features = feature_genes, reduction = "umap", ncol = 4, order = TRUE)
save_plot(p_feature, "Figure2B_feature_plots.pdf", 14, 4)

p_umap_celltype <- Seurat::DimPlot(
  pnd,
  reduction = "umap",
  group.by = "celltype",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
)
save_plot(p_umap_celltype, "Figure2C_annotated_UMAP.pdf", 10, 8)

pooled_prop <- pnd@meta.data |>
  dplyr::count(group, celltype, name = "n_cells") |>
  dplyr::group_by(group) |>
  dplyr::mutate(percent = 100 * n_cells / sum(n_cells)) |>
  dplyr::ungroup()
write_table(pooled_prop, "Figure2D_pooled_cell_proportions.csv")

p_pooled <- ggplot2::ggplot(pooled_prop, ggplot2::aes(celltype, percent, fill = group)) +
  ggplot2::geom_col(position = "dodge") +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values = GROUP_COLORS) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::labs(x = NULL, y = "Pooled cell proportion (%)", fill = NULL)
save_plot(p_pooled, "Figure2D_pooled_cell_proportions.pdf", 8, 7)

saveRDS(pnd, path("results", "objects", "02_pnd_annotated.rds"))
