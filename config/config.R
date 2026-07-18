# Central configuration -------------------------------------------------------
PROJECT_ROOT <- normalizePath(Sys.getenv("PND_PROJECT_ROOT", unset = "."), mustWork = FALSE)

path <- function(...) file.path(PROJECT_ROOT, ...)

dir.create(path("results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(path("results", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(path("results", "objects"), recursive = TRUE, showWarnings = FALSE)
dir.create(path("results", "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(path("environment"), recursive = TRUE, showWarnings = FALSE)

SEED <- 123L
set.seed(SEED)

CONTROL_SAMPLES <- c("GSM8281758", "GSM8281759", "GSM8281760")
SURGERY_SAMPLES <- c("GSM8281761", "GSM8281762", "GSM8281763")
SAMPLE_ORDER <- c(CONTROL_SAMPLES, SURGERY_SAMPLES)
GROUP_LEVELS <- c("Control", "Surgery")

QC_MIN_FEATURES <- 500L
QC_MAX_FEATURES <- 7000L
QC_MAX_MT <- 20
N_VARIABLE_FEATURES <- 3000L
N_HARMONY_DIMS <- 30L
CLUSTER_RESOLUTION <- 0.5
N_OLIGO_DIMS <- 20L

GROUP_COLORS <- c(Control = "#4E79A7", Surgery = "#1BB3B7")
OLIGO_COLORS <- c(
  OPCs = "#59A14F",
  `Newly formed Oligodendrocytes` = "#F28E2B",
  `Mature Oligodendrocytes` = "#8E63CE"
)

MYELIN_GENES <- c("Plp1", "Mbp", "Mobp", "Mog", "Mag", "Opalin")
MYELIN_DISPLAY_GENES <- c("Plp1", "Mbp", "Mobp", "Mog", "Mag", "Opalin", "Cnp", "Mal", "Ermn", "Fa2h", "Ugt8a")
OPC_MARKERS <- c("Pdgfra", "Cspg4", "Sox10", "Olig1", "Olig2")
NFOL_MARKERS <- c("Tcf7l2", "Bcas1", "Enpp6", "Nkx2-2")

# AddModuleScore parameters reported in STAR Methods.
MODULE_NBIN <- 24L
MODULE_CTRL <- 100L

# Manual cluster labels. Confirm against marker plots before publication.
CLUSTER_LABELS <- c(
  `0` = "Activated Microglia",
  `1` = "Inflammatory Microglia",
  `2` = "Mature Oligodendrocytes",
  `3` = "Newly formed Oligodendrocytes",
  `4` = "Astrocytes",
  `5` = "Choroid Plexus Epithelial Cells",
  `6` = "OPCs",
  `7` = "Microglia",
  `8` = "Astrocytes",
  `9` = "Endothelial Cells",
  `10` = "Border-associated Macrophages",
  `11` = "Ependymal Cells",
  `12` = "Neutrophils",
  `13` = "T Cells",
  `14` = "GABAergic Neurons",
  `15` = "Glutamatergic Neurons",
  `16` = "B Cells",
  `17` = "Choroid Plexus-related Epithelial Cells",
  `18` = "Pericytes"
)
