#!/usr/bin/env Rscript
source("config/config.R")
source("R/helpers.R")
assert_packages(c("Seurat", "DESeq2", "Matrix", "dplyr", "tibble", "tidyr", "readr", "ggplot2"))

oligo <- readRDS(path("results", "objects", "03_oligodendrocyte_lineage.rds"))
mature <- subset(oligo, subset = celltype == "Mature Oligodendrocytes")

pb_counts <- aggregate_counts_by_sample(mature, sample_col = "sample")
pb_meta <- sample_metadata(mature)
pb_meta <- pb_meta[match(colnames(pb_counts), pb_meta$sample), , drop = FALSE]
rownames(pb_meta) <- pb_meta$sample

if (!identical(colnames(pb_counts), rownames(pb_meta))) stop("Pseudobulk count and metadata order mismatch.")

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(pb_counts),
  colData = pb_meta,
  design = ~ group
)
keep <- rowSums(DESeq2::counts(dds) >= 10) >= 2
dds <- dds[keep, ]
dds <- DESeq2::DESeq(dds, quiet = TRUE)
res <- DESeq2::results(dds, contrast = c("group", "Surgery", "Control")) |>
  as.data.frame() |>
  tibble::rownames_to_column("gene") |>
  dplyr::arrange(padj)
write_table(res, "Table_S2_mature_oligodendrocyte_pseudobulk_DE.csv")
saveRDS(dds, path("results", "objects", "04_mature_oligodendrocyte_dds.rds"))

vsd <- DESeq2::vst(dds, blind = FALSE)
vst_mat <- SummarizedExperiment::assay(vsd)
display_genes <- intersect(MYELIN_DISPLAY_GENES, rownames(vst_mat))
heat <- vst_mat[display_genes, , drop = FALSE] |>
  as.data.frame() |>
  tibble::rownames_to_column("gene") |>
  tidyr::pivot_longer(-gene, names_to = "sample", values_to = "vst_expression") |>
  dplyr::left_join(pb_meta, by = "sample") |>
  dplyr::group_by(gene) |>
  dplyr::mutate(row_z_score = as.numeric(scale(vst_expression))) |>
  dplyr::ungroup()
heat$sample <- factor(heat$sample, levels = rev(SAMPLE_ORDER))
heat$gene <- factor(heat$gene, levels = display_genes)
write_table(heat, "Figure3D_myelin_pseudobulk_heatmap_source.csv")

p_heat <- ggplot2::ggplot(heat, ggplot2::aes(gene, sample, fill = row_z_score)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.4) +
  ggplot2::scale_fill_gradient2(low = "#4E79A7", mid = "white", high = "#E15759", midpoint = 0, name = "Row Z-score") +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "italic")) +
  ggplot2::labs(x = NULL, y = NULL)
save_plot(p_heat, "Figure3D_myelin_pseudobulk_heatmap.pdf", 8, 5)

myelin_stats <- res |>
  dplyr::filter(gene %in% display_genes) |>
  dplyr::mutate(
    significance = dplyr::case_when(
      !is.na(padj) & padj < 0.05 ~ "FDR < 0.05",
      !is.na(pvalue) & pvalue < 0.05 ~ "Nominal p < 0.05; FDR >= 0.05",
      TRUE ~ "Not significant"
    )
  ) |>
  dplyr::arrange(log2FoldChange)
write_table(myelin_stats, "Figure3E_myelin_gene_statistics.csv")
myelin_stats$gene <- factor(myelin_stats$gene, levels = myelin_stats$gene)

p_lollipop <- ggplot2::ggplot(myelin_stats, ggplot2::aes(log2FoldChange, gene)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = log2FoldChange, yend = gene), color = "grey60") +
  ggplot2::geom_point(ggplot2::aes(color = significance), size = 3) +
  ggplot2::scale_color_manual(values = c(
    "FDR < 0.05" = "#D55E5E",
    "Nominal p < 0.05; FDR >= 0.05" = "#F28E2B",
    "Not significant" = "#4E79A7"
  )) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic"), legend.title = ggplot2::element_blank()) +
  ggplot2::labs(x = "log2 fold change (Surgery vs Control)", y = NULL)
save_plot(p_lollipop, "Figure3E_myelin_gene_lollipop.pdf", 7, 5)
