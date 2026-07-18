#!/usr/bin/env Rscript
source("config/config.R")
source("R/helpers.R")
assert_packages(c("dplyr", "purrr", "readr", "tibble", "tidyr", "clusterProfiler", "org.Mm.eg.db"))

bulk_dir <- path("data", "raw", "GSE174412")
files <- list.files(bulk_dir, pattern = "GSM.*\\.txt(\\.gz)?$", full.names = TRUE)
if (!length(files)) stop("Place GSE174412 processed sample files under data/raw/GSE174412/.")

read_one <- function(file) {
  con <- if (grepl("\\.gz$", file)) gzfile(file) else file
  x <- read.delim(con, header = TRUE, sep = "\t", check.names = FALSE)
  if (ncol(x) < 2) stop("Expected at least two columns in ", basename(file))
  out <- x[, 1:2]
  names(out) <- c("gene_id", tools::file_path_sans_ext(tools::file_path_sans_ext(basename(file))))
  out
}
expr <- purrr::reduce(lapply(files, read_one), dplyr::full_join, by = "gene_id")

# Remove Ensembl version suffix before mapping.
expr$gene_id_clean <- sub("\\..*$", "", expr$gene_id)
map <- clusterProfiler::bitr(unique(expr$gene_id_clean), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
expr_symbol <- expr |>
  dplyr::left_join(map, by = c("gene_id_clean" = "ENSEMBL")) |>
  dplyr::filter(!is.na(SYMBOL))
write_table(expr_symbol, "FigureS2_GSE174412_mapped_expression.csv")

message("Complete the sample-to-group mapping in data/metadata/GSE174412_samples.csv, then rerun the statistical section.")
