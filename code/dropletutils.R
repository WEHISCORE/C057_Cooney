# Create SingleCellExperiment object from 10X data and identify empty droplets
# using DropletUtils.
# Peter Hickey
# 2020-06-25

# Setup ------------------------------------------------------------------------

library(DropletUtils)
library(here)

dir.create(here("data", "SCEs"), recursive = TRUE)
dir.create(here("data", "emptyDrops"))

# Construct SingleCellExperiment object ----------------------------------------

# GEX
sce <- read10xCounts(
  samples = here(
    "extdata",
    "jamesC_10X_040919",
    "CellRanger_GEX_HTO",
    "HTLV_GEX_HTO",
    "outs",
    "raw_feature_bc_matrix"))
sce$project <- factor("C057_Cooney")
stopifnot(!anyDuplicated(sce$Barcode))
colnames(sce) <- paste0(sce$project, ".", sce$Barcode)
sce <- splitAltExps(
  sce,
  rowData(sce)$Type,
  "Gene Expression")

# VDJ
tcr <- read.csv(
  here("data", "CellRanger", "HTLV.filtered_contig_annotations.csv"))

tra <- tcr[tcr$chain == "TRA", ]
trb <- tcr[tcr$chain == "TRB", ]
sce$TRA <- split(DataFrame(tra), factor(tra$barcode, sce$Barcode))
sce$TRB <- split(DataFrame(trb), factor(trb$barcode, sce$Barcode))

# Identify empty droplets ------------------------------------------------------

set.seed(666)
empties <- emptyDrops(counts(sce), lower = 100)

# Check if more permutations are needed; see
# https://osca.bioconductor.org/quality-control.html#testing-for-empty-droplets
more_permutations_needed <- table(
  Sig = empties$FDR <= 0.001,
  Limited = empties$Limited)[1, 2] > 0
stopifnot(!more_permutations_needed)

# Save outputs -----------------------------------------------------------------

saveRDS(
  sce,
  file = here("data", "SCEs", "C057_Cooney.CellRanger.SCE.rds"),
  compress = "xz")

saveRDS(
  object = empties,
  file = here("data", "emptyDrops", "C057_Cooney.emptyDrops.rds"),
  compress = "xz")
writeLines(
  text = sce[, which(empties$FDR <= 0.001)][["Barcode"]],
  con = here("data", "emptyDrops", "C057_Cooney.barcodes.txt"))
