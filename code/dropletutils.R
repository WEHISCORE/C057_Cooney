# Create SingleCellExperiment objects from 10X data using DropletUtils.
# Peter Hickey
# 2019-06-04

# Setup ------------------------------------------------------------------------

library(DropletUtils)
library(here)

# Key variables ----------------------------------------------------------------

dir.create(here("data", "SCEs"), recursive = TRUE)

# Create and save SingleCellExperiment object(s) -------------------------------

sce <- read10xCounts(
  samples = here(
    "extdata",
    "jamesC_10X_040919",
    "CellRanger_GEX_HTO",
    "HTLV_GEX_HTO",
    "outs",
    "raw_feature_bc_matrix"))
saveRDS(
  object = sce,
  file = here("data", "SCEs", "HTLV_GEX_HTO.CellRanger.SCE.rds"),
  compress = "xz")

# TODO: VDJ?
