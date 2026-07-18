assert_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing packages: ", paste(missing, collapse = ", "),
         "\nInstall them before running the pipeline.", call. = FALSE)
  }
}

read_seurat_counts <- function(object, assay = "RNA") {
  tryCatch(
    SeuratObject::GetAssayData(object, assay = assay, layer = "counts"),
    error = function(e) SeuratObject::GetAssayData(object, assay = assay, slot = "counts")
  )
}

read_seurat_data <- function(object, assay = "RNA") {
  tryCatch(
    SeuratObject::GetAssayData(object, assay = assay, layer = "data"),
    error = function(e) SeuratObject::GetAssayData(object, assay = assay, slot = "data")
  )
}

safe_genes <- function(genes, object) {
  intersect(unique(genes), rownames(object))
}

aggregate_counts_by_sample <- function(object, sample_col = "sample", assay = "RNA") {
  counts <- read_seurat_counts(object, assay)
  sample_id <- object[[sample_col, drop = TRUE]]
  if (anyNA(sample_id)) stop("Missing sample identifiers in metadata column: ", sample_col)
  out <- vapply(unique(sample_id), function(s) {
    Matrix::rowSums(counts[, sample_id == s, drop = FALSE])
  }, FUN.VALUE = numeric(nrow(counts)))
  rownames(out) <- rownames(counts)
  out
}

sample_metadata <- function(object, sample_col = "sample", group_col = "group") {
  md <- object@meta.data[, c(sample_col, group_col), drop = FALSE]
  md <- unique(md)
  names(md) <- c("sample", "group")
  if (anyDuplicated(md$sample)) stop("A sample maps to more than one group.")
  md$group <- factor(md$group, levels = GROUP_LEVELS)
  md
}

add_module_score_reproducible <- function(object, genes, name) {
  genes <- safe_genes(genes, object)
  if (!length(genes)) stop("None of the requested genes are present for module: ", name)
  set.seed(SEED)
  Seurat::AddModuleScore(
    object = object,
    features = list(genes),
    name = name,
    nbin = MODULE_NBIN,
    ctrl = MODULE_CTRL,
    assay = "RNA",
    seed = SEED
  )
}

summarize_score_by_sample <- function(object, score_col) {
  object@meta.data |>
    tibble::rownames_to_column("cell") |>
    dplyr::group_by(sample, group) |>
    dplyr::summarise(
      mean_score = mean(.data[[score_col]], na.rm = TRUE),
      median_score = stats::median(.data[[score_col]], na.rm = TRUE),
      n_cells = dplyr::n(),
      .groups = "drop"
    )
}

wilcox_sample_test <- function(df, value_col) {
  stats::wilcox.test(
    stats::reformulate("group", response = value_col),
    data = df,
    exact = FALSE
  )
}

save_plot <- function(plot, filename, width, height) {
  ggplot2::ggsave(
    filename = path("results", "figures", filename),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf
  )
}

write_table <- function(x, filename) {
  readr::write_csv(x, path("results", "tables", filename), na = "")
}

capture_environment <- function() {
  writeLines(capture.output(sessionInfo()), path("environment", "sessionInfo.txt"))
  if (requireNamespace("renv", quietly = TRUE)) {
    renv::snapshot(project = PROJECT_ROOT, prompt = FALSE)
  }
}
