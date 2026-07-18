#!/usr/bin/env Rscript
source("config/config.R")
source("R/helpers.R")
assert_packages(c("Seurat", "dplyr", "readr", "ggplot2", "ggpubr"))

pnd <- readRDS(path("results", "objects", "02_pnd_annotated.rds"))
oligo_states <- c("OPCs", "Newly formed Oligodendrocytes", "Mature Oligodendrocytes")
oligo <- subset(pnd, subset = celltype %in% oligo_states)
Seurat::DefaultAssay(oligo) <- "RNA"
oligo <- Seurat::NormalizeData(oligo, verbose = FALSE)
oligo <- Seurat::FindVariableFeatures(oligo, nfeatures = N_VARIABLE_FEATURES, verbose = FALSE)
oligo <- Seurat::ScaleData(oligo, verbose = FALSE)
oligo <- Seurat::RunPCA(oligo, verbose = FALSE)
oligo <- Seurat::RunUMAP(oligo, dims = seq_len(N_OLIGO_DIMS), seed.use = SEED, reduction.name = "umap_oligo")

p_group <- Seurat::DimPlot(oligo, reduction = "umap_oligo", group.by = "group", cols = GROUP_COLORS)
p_state <- Seurat::DimPlot(oligo, reduction = "umap_oligo", group.by = "celltype", cols = OLIGO_COLORS, label = TRUE, repel = TRUE)
save_plot(p_group, "Figure3A_oligo_UMAP_group.pdf", 6, 5)
save_plot(p_state, "Figure3B_oligo_UMAP_state.pdf", 6, 5)

prop <- oligo@meta.data |>
  dplyr::count(sample, group, celltype, name = "n_cells") |>
  dplyr::group_by(sample) |>
  dplyr::mutate(percent = 100 * n_cells / sum(n_cells)) |>
  dplyr::ungroup()
prop$celltype <- factor(prop$celltype, levels = oligo_states)
write_table(prop, "Figure3C_oligo_state_proportions_per_sample.csv")

prop_stats <- prop |>
  dplyr::group_by(celltype) |>
  dplyr::group_modify(~{
    wt <- wilcox_sample_test(.x, "percent")
    tibble::tibble(p_value = wt$p.value)
  })
write_table(prop_stats, "Figure3C_oligo_state_proportion_statistics.csv")

p_prop <- ggplot2::ggplot(prop, ggplot2::aes(celltype, percent, color = group)) +
  ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = 0.06, dodge.width = 0.65), size = 3) +
  ggplot2::stat_summary(ggplot2::aes(group = group), fun = mean, geom = "crossbar", position = ggplot2::position_dodge(0.65), width = 0.45, color = "black") +
  ggplot2::stat_summary(ggplot2::aes(group = group), fun.data = ggplot2::mean_se, geom = "errorbar", position = ggplot2::position_dodge(0.65), width = 0.15, color = "black") +
  ggplot2::scale_color_manual(values = GROUP_COLORS) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::labs(x = NULL, y = "Cell proportion per sample (%)", color = NULL)
save_plot(p_prop, "Figure3C_oligo_state_proportions.pdf", 7, 5)

saveRDS(oligo, path("results", "objects", "03_oligodendrocyte_lineage.rds"))
