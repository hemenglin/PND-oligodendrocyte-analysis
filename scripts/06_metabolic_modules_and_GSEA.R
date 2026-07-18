#!/usr/bin/env Rscript
source("config/config.R")
source("R/helpers.R")
assert_packages(c("Seurat", "DESeq2", "clusterProfiler", "org.Mm.eg.db", "AnnotationDbi", "dplyr", "readr", "ggplot2"))

oligo <- readRDS(path("results", "objects", "03_oligodendrocyte_lineage.rds"))
mature <- subset(oligo, subset = celltype == "Mature Oligodendrocytes")
res <- readr::read_csv(path("results", "tables", "Table_S2_mature_oligodendrocyte_pseudobulk_DE.csv"), show_col_types = FALSE)

gene_set_file <- path("data", "gene_sets", "metabolic_gene_sets.csv")
if (!file.exists(gene_set_file)) {
  stop("Create data/gene_sets/metabolic_gene_sets.csv with columns module and gene. See README.")
}
gene_sets_df <- readr::read_csv(gene_set_file, show_col_types = FALSE)
gene_sets <- split(gene_sets_df$gene, gene_sets_df$module)

module_source <- lapply(names(gene_sets), function(module_name) {
  obj <- add_module_score_reproducible(mature, gene_sets[[module_name]], paste0("Module_", make.names(module_name)))
  score_col <- grep(paste0("^Module_", make.names(module_name)), colnames(obj@meta.data), value = TRUE)[1]
  summarize_score_by_sample(obj, score_col) |>
    dplyr::mutate(module = module_name)
}) |>
  dplyr::bind_rows()
write_table(module_source, "Figure4A_metabolic_module_scores_per_sample.csv")

module_stats <- module_source |>
  dplyr::group_by(module) |>
  dplyr::group_modify(~tibble::tibble(p_value = wilcox_sample_test(.x, "mean_score")$p.value)) |>
  dplyr::ungroup() |>
  dplyr::mutate(p_adjust = p.adjust(p_value, method = "BH"))
write_table(module_stats, "Figure4A_metabolic_module_statistics.csv")

p_modules <- ggplot2::ggplot(module_source, ggplot2::aes(group, mean_score, color = group)) +
  ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.05), size = 2.5) +
  ggplot2::stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "black") +
  ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "errorbar", width = 0.14, color = "black") +
  ggplot2::facet_wrap(~module, scales = "free_y") +
  ggplot2::scale_color_manual(values = GROUP_COLORS) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(legend.position = "none", strip.background = ggplot2::element_blank()) +
  ggplot2::labs(x = NULL, y = "Mean module score per sample")
save_plot(p_modules, "Figure4A_metabolic_module_scores.pdf", 10, 7)

# Rank by DESeq2 Wald statistic, as reported in STAR Methods.
rank_df <- res |>
  dplyr::filter(!is.na(gene), !is.na(stat)) |>
  dplyr::distinct(gene, .keep_all = TRUE)
map <- clusterProfiler::bitr(rank_df$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
rank_mapped <- rank_df |>
  dplyr::inner_join(map, by = c("gene" = "SYMBOL")) |>
  dplyr::group_by(ENTREZID) |>
  dplyr::slice_max(abs(stat), n = 1, with_ties = FALSE) |>
  dplyr::ungroup()
gene_list <- rank_mapped$stat
names(gene_list) <- rank_mapped$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

gsea <- clusterProfiler::gseGO(
  geneList = gene_list,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  keyType = "ENTREZID",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE,
  seed = TRUE
)
gsea_df <- as.data.frame(gsea)
write_table(gsea_df, "Table_S4_GO_BP_GSEA_complete.csv")
saveRDS(gsea, path("results", "objects", "06_GO_BP_GSEA.rds"))
