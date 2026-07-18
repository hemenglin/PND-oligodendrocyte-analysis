#!/usr/bin/env Rscript
source("config/config.R")
source("R/helpers.R")
assert_packages(c("Seurat", "dplyr", "tidyr", "readr", "ggplot2", "ggpubr", "rstatix"))

oligo <- readRDS(path("results", "objects", "03_oligodendrocyte_lineage.rds"))
celltypes <- c("Mature Oligodendrocytes", "OPCs", "Newly formed Oligodendrocytes")
labels <- c("Mature oligodendrocytes", "OPCs", "NFOLs")

all_scores <- lapply(seq_along(celltypes), function(i) {
  obj <- subset(oligo, subset = celltype == celltypes[[i]])
  obj <- add_module_score_reproducible(obj, MYELIN_GENES, "MyelinModule")
  out <- summarize_score_by_sample(obj, "MyelinModule1")
  out$cell_population <- labels[[i]]
  out
}) |>
  dplyr::bind_rows()
all_scores$cell_population <- factor(all_scores$cell_population, levels = labels)
write_table(all_scores, "Figure3F-H_myelin_module_scores_per_sample.csv")

stats <- all_scores |>
  dplyr::group_by(cell_population) |>
  dplyr::group_modify(~{
    wt <- wilcox_sample_test(.x, "mean_score")
    tibble::tibble(p_value = wt$p.value, method = wt$method)
  }) |>
  dplyr::ungroup()
write_table(stats, "Figure3F-H_myelin_module_statistics.csv")

p <- ggplot2::ggplot(all_scores, ggplot2::aes(group, mean_score, color = group)) +
  ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.05), size = 3) +
  ggplot2::stat_summary(fun = mean, geom = "crossbar", width = 0.42, color = "black") +
  ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "errorbar", width = 0.15, color = "black") +
  ggplot2::facet_wrap(~cell_population, scales = "free_y", nrow = 1) +
  ggplot2::scale_color_manual(values = GROUP_COLORS) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(legend.position = "none", strip.background = ggplot2::element_blank()) +
  ggplot2::labs(x = NULL, y = "Mean myelin module score per sample")
save_plot(p, "Figure3F-H_myelin_module_scores.pdf", 10, 4)

# Supplementary sample-level marker expression; cells are not statistical replicates.
marker_specs <- list(OPCs = OPC_MARKERS, NFOLs = NFOL_MARKERS)
marker_objects <- list(
  OPCs = subset(oligo, subset = celltype == "OPCs"),
  NFOLs = subset(oligo, subset = celltype == "Newly formed Oligodendrocytes")
)
marker_source <- lapply(names(marker_specs), function(pop) {
  genes <- safe_genes(marker_specs[[pop]], marker_objects[[pop]])
  Seurat::FetchData(marker_objects[[pop]], vars = c("sample", "group", genes)) |>
    tibble::as_tibble() |>
    tidyr::pivot_longer(dplyr::all_of(genes), names_to = "gene", values_to = "expression") |>
    dplyr::group_by(sample, group, gene) |>
    dplyr::summarise(mean_expression = mean(expression), .groups = "drop") |>
    dplyr::mutate(population = pop)
}) |>
  dplyr::bind_rows()
write_table(marker_source, "FigureS1_marker_expression_per_sample.csv")
