#!/usr/bin/env Rscript
source("config/config.R")
source("R/helpers.R")
assert_packages(c("Seurat", "harmony", "dplyr", "readr", "ggplot2", "patchwork", "Matrix"))

set.seed(SEED)
raw_dir <- path("data", "raw", "GSE267933")
if (!dir.exists(raw_dir)) {
  stop("Place the six 10x sample directories under data/raw/GSE267933/ before running.")
}

samples <- SAMPLE_ORDER
missing_dirs <- samples[!dir.exists(file.path(raw_dir, samples))]
if (length(missing_dirs)) stop("Missing sample directories: ", paste(missing_dirs, collapse = ", "))

seurat_list <- lapply(samples, function(s) {
  counts <- Seurat::Read10X(file.path(raw_dir, s))
  obj <- Seurat::CreateSeuratObject(
    counts = counts,
    project = s,
    min.cells = 3,
    min.features = 200
  )
  obj$sample <- s
  obj$group <- ifelse(s %in% CONTROL_SAMPLES, "Control", "Surgery")
  obj
})
names(seurat_list) <- samples

pnd <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = samples)
pnd$sample <- factor(pnd$sample, levels = SAMPLE_ORDER)
pnd$group <- factor(pnd$group, levels = GROUP_LEVELS)
pnd[["percent.mt"]] <- Seurat::PercentageFeatureSet(pnd, pattern = "^mt-")

qc_before <- pnd@meta.data |>
  tibble::rownames_to_column("cell")
write_table(qc_before, "Figure1_QC_before_filtering.csv")

p_qc_before <- Seurat::VlnPlot(
  pnd,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample",
  ncol = 3,
  pt.size = 0
)
save_plot(p_qc_before, "Figure1A_QC_before.pdf", 14, 5)

pnd <- subset(
  pnd,
  subset = nFeature_RNA > QC_MIN_FEATURES &
    nFeature_RNA < QC_MAX_FEATURES &
    percent.mt < QC_MAX_MT
)

qc_after <- pnd@meta.data |>
  tibble::rownames_to_column("cell")
write_table(qc_after, "Figure1_QC_after_filtering.csv")

p_qc_after <- Seurat::VlnPlot(
  pnd,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample",
  ncol = 3,
  pt.size = 0
)
save_plot(p_qc_after, "Figure1B_QC_after.pdf", 14, 5)

pnd <- Seurat::NormalizeData(pnd, verbose = FALSE)
pnd <- Seurat::FindVariableFeatures(
  pnd,
  selection.method = "vst",
  nfeatures = N_VARIABLE_FEATURES,
  verbose = FALSE
)
pnd <- Seurat::ScaleData(pnd, vars.to.regress = "percent.mt", verbose = FALSE)
pnd <- Seurat::RunPCA(pnd, verbose = FALSE)

p_elbow <- Seurat::ElbowPlot(pnd, ndims = 50)
save_plot(p_elbow, "Figure1C_ElbowPlot.pdf", 6, 4)

p_pca_group <- Seurat::DimPlot(pnd, reduction = "pca", group.by = "group", cols = GROUP_COLORS)
p_pca_sample <- Seurat::DimPlot(pnd, reduction = "pca", group.by = "sample")
save_plot(p_pca_group, "Figure1D_PCA_group.pdf", 6, 5)
save_plot(p_pca_sample, "Figure1E_PCA_sample.pdf", 7, 5)

pnd <- harmony::RunHarmony(pnd, group.by.vars = "sample", verbose = FALSE)
pnd <- Seurat::RunUMAP(pnd, reduction = "harmony", dims = seq_len(N_HARMONY_DIMS), seed.use = SEED)
pnd <- Seurat::FindNeighbors(pnd, reduction = "harmony", dims = seq_len(N_HARMONY_DIMS), verbose = FALSE)
pnd <- Seurat::FindClusters(pnd, resolution = CLUSTER_RESOLUTION, random.seed = SEED, verbose = FALSE)

saveRDS(pnd, path("results", "objects", "01_pnd_clustered.rds"))
capture_environment()
