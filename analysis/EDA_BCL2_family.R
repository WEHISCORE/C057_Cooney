# EDA of BCL2 family as requested by James
# 2024-09-20
# Peter Hickey

# TODO: Need to resolve the device issue thing (dev.off(), etc.) so that I can
#       knit this document, but not urgent.


# NOTE: This script does not use the renv analysis environment used by the rest
#       of this project. Instead, it used R 4.4, BioC 3.19, and edgeR 4.2.1.
renv::deactivate()

library(here)
library(scater)
library(edgeR)
library(scater)

# BCL family members provided by James
bcl2_df <- read.csv(here("data/group-1057.csv"))

outdir <- here("output/BCL2_family")
dir.create(outdir)

sce <- readRDS(here("data/SCEs/C057_Cooney.annotated.SCE.rds"))

# Some useful colours
sample_colours <- setNames(
  unique(sce$sample_colours),
  unique(names(sce$sample_colours)))
treatment_colours <- setNames(
  unique(sce$treatment_colours),
  unique(names(sce$treatment_colours)))
cluster_colours <- setNames(
  unique(sce$cluster_colours),
  unique(names(sce$cluster_colours)))

# Single-cell ------------------------------------------------------------------

# NOTE: dev.off()-shenanigans required to get PDFs to properly save within the
#       rmarkdown render (similar to
#       https://github.com/raivokolde/pheatmap/issues/37 but I came up with a
#       different solution)
dev.list()
p <- plotExpression(
  sce,
  bcl2_df$Approved.symbol,
  x = "Sample",
  colour_by = "Treatment",
  ncol = 5) +
  scale_colour_manual(values = treatment_colours, name = "Treatment") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ggtitle("Ignoring cycling")
dev.list()
ggsave(
  file.path(outdir, "ignoring_cycling_subset.single_cell.pdf"),
  p,
  width = 12,
  height = 10)
dev.list()
p
dev.list()
dev.off(3)

p <- plotExpression(
  sce[, sce$cycling_subset == "Not cycling"],
  bcl2_df$Approved.symbol,
  x = "Sample",
  colour_by = "Treatment",
  ncol = 5) +
  scale_colour_manual(values = treatment_colours, name = "Treatment") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ggtitle("Not cycling")
ggsave(
  file.path(outdir, "not_cycling.single_cell.pdf"),
  p,
  width = 12,
  height = 10)
dev.list()
p
dev.list()
dev.off(3)

p <- plotExpression(
  sce[, sce$cycling_subset == "Cycling"],
  bcl2_df$Approved.symbol,
  x = "Sample",
  colour_by = "Treatment",
  ncol = 5) +
  scale_colour_manual(values = treatment_colours, name = "Treatment") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ggtitle("Cycling")
ggsave(
  file.path(outdir, "cycling.single_cell.pdf"),
  p,
  width = 12,
  height = 10)
dev.list()
p
dev.list()
dev.off(3)

# Pseudobulk: Ignoring `cycling_subset` ----------------------------------------

se <- aggregateAcrossCells(
  sce,
  id = colData(sce)[, c("Sample", "Treatment")],
  coldata.merge = FALSE,
  use.dimred = FALSE,
  use.altexps = FALSE)
se <- logNormCounts(se)

plotHeatmap(
  se,
  features = bcl2_df$Approved.symbol,
  center = FALSE,
  color = hcl.colors(101, "Inferno"),
  order_columns_by = c("Treatment", "Sample"),
  column_annotation_colors = list(
    Treatment = treatment_colours,
    Sample = sample_colours),
  cluster_rows = TRUE,
  fontsize = 8,
  filename = file.path(
    outdir,
    "ignoring_cycling_subset.pseudobulk.heatmap.pdf"))
dev.list()
dev.off(3)
dev.list()

plotHeatmap(
  se,
  features = bcl2_df$Approved.symbol,
  center = TRUE,
  scale = TRUE,
  zlim = c(-2, 2),
  symmetric = TRUE,
  color = hcl.colors(101, "Blue-Red 3"),
  order_columns_by = c("Treatment", "Sample"),
  column_annotation_colors = list(
    Treatment = treatment_colours,
    Sample = sample_colours),
  cluster_rows = TRUE,
  fontsize = 8,
  filename = file.path(
    outdir,
    "ignoring_cycling_subset.pseudobulk.row-normalized_heatmap.pdf"))
dev.list()
dev.off(3)
dev.list()

y <- SE2DGEList(se)
colnames(y) <- se$Sample
y$samples$Treatment <- relevel(y$samples$Treatment, "Uninfected")
design <- model.matrix(~Treatment, y$samples)

keep0 <- filterByExpr(y, design = design)
keep <- keep0
keep[bcl2_df$Approved.symbol] <- TRUE
# NOTE: Some versions of this analysis included all BCL family members, even
#       those that did not pass the filterByExpr() filter. However, these genes
#       by definition have so few counts as to be meaningless.
table(keep0, keep)
y$counts[keep & !keep0, ]
y <- y[keep0, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, "TreatmentInfected")

y_index <- ids2indices(
  list(`BCL2-family` = bcl2_df$Approved.symbol),
  rownames(y))
fry(
  y,
  index = y_index,
  design = design,
  coef = "TreatmentInfected")

pdf(
  file.path(
    outdir,
    "ignoring_cycling_subset.pseudobulk.gene_set_test.pdf"),
  width = 8,
  height = 4)
par(mfrow = c(1, 2))
y_status <- rep("Other", nrow(lrt))
y_status[rownames(lrt) %in% bcl2_df$Approved.symbol] <- "In"
plotMD(
  lrt,
  status = y_status,
  values = "In",
  hl.col = "red",
  legend = "bottomright",
  main = "Infected vs. Uninfected")
abline(h = 0, col = "darkgrey")
barcodeplot(
  statistics = lrt$table$logFC,
  index = rownames(lrt) %in% bcl2_df$Approved.symbol)
dev.off()

# Non-cycling subset -----------------------------------------------------------

not_cycling_sce <- readRDS(
  here("data", "SCEs", "C057_Cooney.not_cycling.annotated.SCE.rds"))

not_cycling_se <- aggregateAcrossCells(
  not_cycling_sce,
  id = colData(not_cycling_sce)[, c("Sample", "Treatment")],
  coldata.merge = FALSE,
  use.dimred = FALSE,
  use.altexps = FALSE)
not_cycling_se <- logNormCounts(not_cycling_se)

plotHeatmap(
  not_cycling_se,
  features = bcl2_df$Approved.symbol,
  center = FALSE,
  color = hcl.colors(101, "Inferno"),
  order_columns_by = c("Treatment", "Sample"),
  column_annotation_colors = list(
    Treatment = treatment_colours,
    Sample = sample_colours),
  cluster_rows = TRUE,
  fontsize = 8,
  filename = file.path(
    outdir,
    "non-cycling_subset.pseudobulk.heatmap.pdf"))

plotHeatmap(
  not_cycling_se,
  # NOTE: Have to exclude any genes not expressed in any sample.
  features = bcl2_df$Approved.symbol[
    rowSums(counts(not_cycling_se)[bcl2_df$Approved.symbol, ]) > 0],
  center = TRUE,
  scale = TRUE,
  zlim = c(-2, 2),
  symmetric = TRUE,
  color = hcl.colors(101, "Blue-Red 3"),
  order_columns_by = c("Treatment", "Sample"),
  column_annotation_colors = list(
    Treatment = treatment_colours,
    Sample = sample_colours),
  cluster_rows = TRUE,
  fontsize = 8,
  filename = file.path(
    outdir,
    "non-cycling_subset.pseudobulk.row-normalized_heatmap.pdf"))

y <- SE2DGEList(not_cycling_se)
colnames(y) <- not_cycling_se$Sample
y$samples$Treatment <- relevel(y$samples$Treatment, "Uninfected")
design <- model.matrix(~Treatment, y$samples)

keep0 <- filterByExpr(y, design = design)
keep <- keep0
keep[bcl2_df$Approved.symbol] <- TRUE
# NOTE: Some versions of this analysis included all BCL family members, even
#       those that did not pass the filterByExpr() filter. However, these genes
#       by definition have so few counts as to be meaningless.
table(keep0, keep)
y$counts[keep & !keep0, ]
y <- y[keep0, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, "TreatmentInfected")

y_index <- ids2indices(
  list(`BCL2-family` = bcl2_df$Approved.symbol),
  rownames(y))
fry(
  y,
  index = y_index,
  design = design,
  coef = "TreatmentInfected")

pdf(
  file.path(
    outdir,
    "non-cycling_subset.pseudobulk.gene_set_test.pdf"),
  width = 8,
  height = 4)
par(mfrow = c(1, 2))
y_status <- rep("Other", nrow(lrt))
y_status[rownames(lrt) %in% bcl2_df$Approved.symbol] <- "In"
plotMD(
  lrt,
  status = y_status,
  values = "In",
  hl.col = "red",
  legend = "bottomright")
abline(h = 0, col = "darkgrey")
barcodeplot(
  statistics = lrt$table$logFC,
  index = rownames(lrt) %in% bcl2_df$Approved.symbol)
dev.off()

# Cycling subset ---------------------------------------------------------------

cycling_sce <- readRDS(
  here("data", "SCEs", "C057_Cooney.cycling.annotated.SCE.rds"))

cycling_se <- aggregateAcrossCells(
  cycling_sce,
  id = colData(cycling_sce)[, c("Sample", "Treatment")],
  coldata.merge = FALSE,
  use.dimred = FALSE,
  use.altexps = FALSE)
cycling_se <- logNormCounts(cycling_se)

plotHeatmap(
  cycling_se,
  features = bcl2_df$Approved.symbol,
  center = FALSE,
  color = hcl.colors(101, "Inferno"),
  order_columns_by = c("Treatment", "Sample"),
  column_annotation_colors = list(
    Treatment = treatment_colours,
    Sample = sample_colours),
  cluster_rows = TRUE,
  fontsize = 8,
  filename = file.path(
    outdir,
    "cycling_subset.pseudobulk.heatmap.pdf"))

plotHeatmap(
  cycling_se,
  # NOTE: Have to exclude any genes not expressed in any sample.
  features = bcl2_df$Approved.symbol[
    rowSums(counts(cycling_se)[bcl2_df$Approved.symbol, ]) > 0],
  center = TRUE,
  scale = TRUE,
  zlim = c(-2, 2),
  symmetric = TRUE,
  color = hcl.colors(101, "Blue-Red 3"),
  order_columns_by = c("Treatment", "Sample"),
  column_annotation_colors = list(
    Treatment = treatment_colours,
    Sample = sample_colours),
  cluster_rows = TRUE,
  fontsize = 8,
  filename = file.path(
    outdir,
    "cycling_subset.pseudobulk.row-normalized_heatmap.pdf"))

y <- SE2DGEList(cycling_se)
colnames(y) <- cycling_se$Sample
y$samples$Treatment <- relevel(y$samples$Treatment, "Uninfected")
design <- model.matrix(~Treatment, y$samples)

keep0 <- filterByExpr(y, design = design)
keep <- keep0
keep[bcl2_df$Approved.symbol] <- TRUE
# NOTE: Some versions of this analysis included all BCL family members, even
#       those that did not pass the filterByExpr() filter. However, these genes
#       by definition have so few counts as to be meaningless.
table(keep0, keep)
y$counts[keep & !keep0, ]
y <- y[keep0, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, "TreatmentInfected")

y_index <- ids2indices(
  list(`BCL2-family` = bcl2_df$Approved.symbol),
  rownames(y))
fry(
  y,
  index = y_index,
  design = design,
  coef = "TreatmentInfected")

pdf(
  file.path(
    outdir,
    "cycling_subset.pseudobulk.gene_set_test.pdf"),
  width = 8,
  height = 4)
par(mfrow = c(1, 2))
y_status <- rep("Other", nrow(lrt))
y_status[rownames(lrt) %in% bcl2_df$Approved.symbol] <- "In"
plotMD(
  lrt,
  status = y_status,
  values = "In",
  hl.col = "red",
  legend = "bottomright")
abline(h = 0, col = "darkgrey")
barcodeplot(
  statistics = lrt$table$logFC,
  index = rownames(lrt) %in% bcl2_df$Approved.symbol)
dev.off()

# Session info -----------------------------------------------------------------

sessionInfo()
