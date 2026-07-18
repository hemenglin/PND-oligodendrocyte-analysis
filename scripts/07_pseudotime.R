#!/usr/bin/env Rscript
source("config/config.R")
source("R/helpers.R")
assert_packages(c("Seurat", "monocle3", "SingleCellExperiment", "dplyr", "tibble", "readr", "ggplot2"))

oligo <- readRDS(path("results", "objects", "03_oligodendrocyte_lineage.rds"))
focus <- subset(oligo, subset = celltype %in% c("Newly formed Oligodendrocytes", "Mature Oligodendrocytes"))
counts <- read_seurat_counts(focus)
cell_metadata <- focus@meta.data
gene_metadata <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))
cds <- monocle3::new_cell_data_set(counts, cell_metadata = cell_metadata, gene_metadata = gene_metadata)

umap <- Seurat::Embeddings(focus, reduction = "umap_oligo")
SingleCellExperiment::reducedDims(cds)$UMAP <- umap[colnames(cds), , drop = FALSE]
cds <- monocle3::cluster_cells(cds, reduction_method = "UMAP", random_seed = SEED)
cds <- monocle3::learn_graph(cds, use_partition = TRUE)

# Select the principal graph node most frequently nearest to NFOL cells.
closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex <- as.matrix(closest_vertex[colnames(cds), , drop = FALSE])[, 1]
nfol_cells <- colnames(cds)[colData(cds)$celltype == "Newly formed Oligodendrocytes"]
root_node <- names(which.max(table(closest_vertex[nfol_cells])))
cds <- monocle3::order_cells(cds, root_pr_nodes = root_node)

pt <- monocle3::pseudotime(cds)
pt_df <- tibble::tibble(
  cell = names(pt),
  pseudotime = as.numeric(pt),
  sample = colData(cds)$sample,
  group = colData(cds)$group,
  celltype = colData(cds)$celltype
)
write_table(pt_df, "Figure5_cell_level_pseudotime.csv")

sample_pt <- pt_df |>
  dplyr::filter(is.finite(pseudotime)) |>
  dplyr::group_by(sample, group) |>
  dplyr::summarise(median_pseudotime = median(pseudotime), n_cells = dplyr::n(), .groups = "drop")
write_table(sample_pt, "Figure5F_sample_median_pseudotime.csv")
wt <- wilcox_sample_test(sample_pt, "median_pseudotime")
write_table(tibble::tibble(p_value = wt$p.value, method = wt$method, root_node = root_node), "Figure5F_pseudotime_statistics.csv")

# Trajectory modules: myelin and NFOL-specific module are distinct analyses.
focus <- add_module_score_reproducible(focus, MYELIN_GENES, "TrajectoryMyelin")
focus <- add_module_score_reproducible(focus, NFOL_MARKERS, "TrajectoryNFOL")
module_df <- focus@meta.data |>
  dplyr::group_by(sample, group) |>
  dplyr::summarise(
    myelin_score = mean(TrajectoryMyelin1),
    nfol_score = mean(TrajectoryNFOL1),
    .groups = "drop"
  )
write_table(module_df, "Figure5I_trajectory_module_scores_per_sample.csv")
saveRDS(cds, path("results", "objects", "07_monocle3_cds.rds"))
