# PND mature oligodendrocyte single-cell reanalysis

This repository contains the reproducible R workflow for the manuscript:

> Single-cell transcriptomics identifies a mitochondrial energy-related transcriptional signature in mature oligodendrocytes in an aged mouse model of postoperative neurocognitive disorder

## Public datasets

- Hippocampal single-cell RNA-seq: GEO **GSE267933**
- Independent prefrontal cortex bulk RNA-seq: GEO **GSE174412**

Raw public data are not redistributed. Download the files from GEO and place them in the directories described below.

## Repository structure

```text
config/                  Central analysis parameters and gene sets
R/                       Reusable helper functions
scripts/                 Numbered analysis scripts
  01_preprocess_atlas.R
  02_annotate_atlas.R
  03_oligodendrocyte_lineage.R
  04_pseudobulk_myelin.R
  05_module_scores.R
  06_metabolic_modules_and_GSEA.R
  07_pseudotime.R
  08_cellchat.R
  09_external_bulk_validation.R
  10_capture_environment.R
data/raw/                 User-downloaded GEO data; ignored by Git
data/metadata/            Sample annotations
data/gene_sets/           Predefined gene sets
results/figures/          Generated figures
results/tables/           Source data and statistical results
results/objects/          Intermediate R objects; ignored by Git
environment/              sessionInfo and renv lockfile
legacy_PND_code.R         Original exploratory script retained for audit history
```

## Required input layout

```text
data/raw/GSE267933/
  GSM8281758/
  GSM8281759/
  GSM8281760/
  GSM8281761/
  GSM8281762/
  GSM8281763/
```

Each sample directory should contain the corresponding 10x matrix files readable by `Seurat::Read10X()`.

Place the processed GSE174412 sample files under:

```text
data/raw/GSE174412/
```

Create `data/metadata/GSE174412_samples.csv` with at least:

```text
sample,group
GSMxxxxxx,Control
GSMxxxxxx,PND
```

## Running the workflow

Run from the repository root. The project deliberately does not use `setwd()` or absolute Windows paths.

```r
Sys.setenv(PND_PROJECT_ROOT = normalizePath("."))
source("scripts/01_preprocess_atlas.R")
```

Run scripts in numerical order, or use:

```bash
Rscript run_all.R
```

Script 02 contains a manual cluster-to-cell-type mapping in `config/config.R`. Confirm this mapping against the generated marker dot plot before freezing a release.

## Statistical unit

Biological samples, not individual cells, are the unit of group-level inference. The workflow therefore uses:

- sample-level lineage proportions;
- sample-aggregated pseudobulk counts for DESeq2;
- sample-averaged module scores;
- sample-level median pseudotime;
- sample-level candidate ligand/receptor expression.

CellChat networks are reconstructed from pooled cells within each group and are interpreted as descriptive and hypothesis-generating, not as replicate-level differential communication tests.

## Myelin module

The same predefined myelin module is scored separately in mature oligodendrocytes, OPCs, and NFOLs:

```text
Plp1, Mbp, Mobp, Mog, Mag, Opalin
```

This replaces the inconsistent exploratory code in which OPC- and NFOL-specific programs were used for Figure 3G/H.

## Reproducibility

- Random seed: `123`
- `AddModuleScore`: `nbin = 24`, `ctrl = 100`
- DESeq2 design: `~ group`
- Pseudobulk gene filter: raw count >= 10 in at least two samples
- GO-BP GSEA ranking statistic: DESeq2 Wald statistic
- Multiple testing: Benjamini-Hochberg where stated

Run `scripts/10_capture_environment.R` after a successful analysis. If `renv` is installed, the script also snapshots the package environment.

## Important checks before public release

1. Verify all cluster labels in `config/config.R` against marker genes.
2. Confirm that the pseudotime root and cell-state distribution support the reported orientation.
3. Fill the complete metabolic gene-set file rather than relying on undocumented objects from an R workspace.
4. Complete GSE174412 sample metadata and the external-validation statistical section.
5. Run the full workflow from a clean R session and compare every generated figure/table with the submitted manuscript.
6. Replace any manuscript placeholder with the final GitHub URL and archived Zenodo DOI.

## License

MIT License. Public GEO data remain subject to their original terms and attribution requirements.
