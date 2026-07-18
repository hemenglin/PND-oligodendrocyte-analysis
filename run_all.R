#!/usr/bin/env Rscript
scripts <- sprintf("scripts/%02d_%s.R", 1:10, c(
  "preprocess_atlas", "annotate_atlas", "oligodendrocyte_lineage",
  "pseudobulk_myelin", "module_scores", "metabolic_modules_and_GSEA",
  "pseudotime", "cellchat", "external_bulk_validation", "capture_environment"
))
for (script in scripts) {
  message("\n===== Running ", script, " =====")
  status <- system2(file.path(R.home("bin"), "Rscript"), script)
  if (status != 0) stop("Pipeline stopped at ", script, call. = FALSE)
}
