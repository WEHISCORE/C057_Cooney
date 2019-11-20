library(here)
library(DropletUtils)

sce <- readRDS(here("data", "SCEs", "HTLV_GEX_HTO.CellRanger.SCE.rds"))
rna <- sce[grepl("ENS", rownames(sce)), ]

set.seed(666)
empties <- emptyDrops(counts(rna), lower = 100)

# NOTE: Check if more permutations are needed (https://bioconductor.org/packages/devel/workflows/vignettes/simpleSingleCell/inst/doc/tenx.html#32_examining_cell-calling_diagnostics)
more_permutations_needed <- table(
  Sig = empties$FDR <= 0.001,
  Limited = empties$Limited)[1, 2] > 0
stopifnot(!more_permutations_needed)

dir.create(here("data", "emptyDrops"))
saveRDS(empties, here("data", "emptyDrops", "empties.rds"))
