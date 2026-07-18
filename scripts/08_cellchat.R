#!/usr/bin/env Rscript
source("config/config.R")
source("R/helpers.R")
assert_packages(c("Seurat", "CellChat", "dplyr", "readr", "future"))

pnd <- readRDS(path("results", "objects", "02_pnd_annotated.rds"))
chat_types <- c(
  "Activated Microglia", "Inflammatory Microglia", "Microglia", "Astrocytes",
  "OPCs", "Newly formed Oligodendrocytes", "Mature Oligodendrocytes"
)
pnd_chat <- subset(pnd, subset = celltype %in% chat_types)

run_cellchat <- function(object, group_name) {
  object <- subset(object, subset = group == group_name)
  expr <- read_seurat_data(object)
  meta <- object@meta.data[, "celltype", drop = FALSE]
  chat <- CellChat::createCellChat(object = expr, meta = meta, group.by = "celltype")
  chat@DB <- CellChat::CellChatDB.mouse
  chat <- CellChat::subsetData(chat)
  future::plan("sequential")
  chat <- CellChat::identifyOverExpressedGenes(chat)
  chat <- CellChat::identifyOverExpressedInteractions(chat)
  chat <- CellChat::computeCommunProb(chat, raw.use = FALSE, population.size = TRUE)
  chat <- CellChat::filterCommunication(chat, min.cells = 10)
  chat <- CellChat::computeCommunProbPathway(chat)
  chat <- CellChat::aggregateNet(chat)
  chat
}

chat_control <- run_cellchat(pnd_chat, "Control")
chat_surgery <- run_cellchat(pnd_chat, "Surgery")
merged <- CellChat::mergeCellChat(list(Control = chat_control, Surgery = chat_surgery), add.names = GROUP_LEVELS)

control_lr <- CellChat::subsetCommunication(chat_control)
surgery_lr <- CellChat::subsetCommunication(chat_surgery)
control_lr$group <- "Control"
surgery_lr$group <- "Surgery"
write_table(dplyr::bind_rows(control_lr, surgery_lr), "Table_S5_CellChat_ligand_receptor_complete.csv")

summary_df <- tibble::tibble(
  group = GROUP_LEVELS,
  interaction_count = c(sum(chat_control@net$count), sum(chat_surgery@net$count)),
  interaction_strength = c(sum(chat_control@net$weight), sum(chat_surgery@net$weight))
)
write_table(summary_df, "Figure6A_CellChat_network_summary.csv")

# Sample-level expression is independent of pooled CellChat inference.
candidate_ligands <- c("F11r", "Grn", "Pros1", "Psap")
candidate_receptors <- c("Jam3", "Sort1", "Tyro3", "Gpr37")
expr_source <- lapply(list(
  ligands = list(types = c("Activated Microglia", "Inflammatory Microglia", "Microglia"), genes = candidate_ligands),
  receptors = list(types = "Mature Oligodendrocytes", genes = candidate_receptors)
), function(spec) {
  obj <- subset(pnd_chat, subset = celltype %in% spec$types)
  genes <- safe_genes(spec$genes, obj)
  Seurat::FetchData(obj, vars = c("sample", "group", genes)) |>
    tibble::as_tibble() |>
    tidyr::pivot_longer(dplyr::all_of(genes), names_to = "gene", values_to = "expression") |>
    dplyr::group_by(sample, group, gene) |>
    dplyr::summarise(mean_expression = mean(expression), .groups = "drop")
}) |>
  dplyr::bind_rows(.id = "panel")
write_table(expr_source, "Figure6F-G_candidate_expression_per_sample.csv")

saveRDS(merged, path("results", "objects", "08_CellChat_merged.rds"))
