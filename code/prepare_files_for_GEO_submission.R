# Prepare files for GEO submission.
# Peter Hickey
# 2020-07-24

library(here)
library(DropletUtils)

outdir <- here("GEO")
dir.create(outdir)

# Raw data -------------------------------------------------------------------

# BAM
file.copy(
  from = here(
    "extdata",
    "jamesC_10X_040919",
    "CellRanger_GEX_HTO",
    "HTLV_GEX_HTO",
    "outs",
    "possorted_genome_bam.bam"),
  to = outdir,
  recursive = FALSE,
  overwrite = FALSE)

# Count matrix data
file.copy(
  from = here(
    "extdata",
    "jamesC_10X_040919",
    "CellRanger_GEX_HTO",
    "HTLV_GEX_HTO",
    "outs",
    "raw_feature_bc_matrix"),
  to = outdir,
  recursive = TRUE,
  overwrite = FALSE)

# Demultiplexed data -----------------------------------------------------------

sce <- readRDS(here("data", "SCEs", "C057_Cooney.demultiplexed.SCE.rds"))

# Gene counts
write10xCounts(
  path = file.path(outdir, "demultiplexed"),
  x = counts(sce),
  barcodes = sce$Barcode,
  gene.id = rowData(sce)$ID,
  gene.symbol = rowData(sce)$Symbol,
  gene.type = rowData(sce)$Type,
  version = "3")

# HTO counts
write10xCounts(
  path = file.path(outdir, "hto"),
  x = as(counts(altExp(sce)), "dgCMatrix"),
  barcodes = sce$Barcode,
  gene.id = rowData(altExp(sce))$ID,
  gene.symbol = rowData(altExp(sce))$Symbol,
  gene.type = rowData(altExp(sce))$Type,
  version = "3")

# colData
coldata_df <- as.data.frame(
  colData(sce)[, c("Barcode", "hto_cluster", "HTO", "Sample", "Treatment")])
write.csv(
  x = coldata_df,
  file = file.path(outdir, "demultiplexed", "coldata.csv.gz"),
  row.names = FALSE,
  quote = FALSE)

# TCR files --------------------------------------------------------------------

file.copy(
  from = here("data", "CellRanger", "HTLV.filtered_contig_annotations.csv"),
  to = file.path(outdir, "filtered_contig_annotations.csv"))
R.utils::gzip(file.path(outdir, "filtered_contig_annotations.csv"))
