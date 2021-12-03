# EDA of BCL2 family as requested by James
# 2021-12-03
# Peter Hickey

library(here)
library(scater)
library(edgeR)

# BCL family members provided by James
bcl2_df <- read.csv(here("data/group-1057.csv"))

outdir <- here("output/BCL2_family")
dir.create(outdir)

# Ignoring `cycling_subset` ----------------------------------------------------

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

plotHeatmap(
  se,
  features = bcl2_df$Approved.symbol,
  center = TRUE,
  # zlim = c(-3, 3),
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

y <- SE2DGEList(se)
colnames(y) <- se$Sample
y$samples$Treatment <- relevel(y$samples$Treatment, "Uninfected")
design <- model.matrix(~Treatment, y$samples)

keep0 <- filterByExpr(y, design = design)
keep <- keep0
keep[bcl2_df$Approved.symbol] <- TRUE
y <- y[keep, ,keep.lib.sizes = FALSE]
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

write.csv(
  cbind(
    as.data.frame(y$counts[bcl2_df$Approved.symbol, ]),
    keep = keep0[bcl2_df$Approved.symbol]),
  file.path(
    outdir, "ignoring_cycling_subset.pseudobulk.raw_counts.csv"))

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
  features = bcl2_df$Approved.symbol,
  center = TRUE,
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
y <- y[keep, ,keep.lib.sizes = FALSE]
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

write.csv(
  cbind(
    as.data.frame(y$counts[bcl2_df$Approved.symbol, ]),
    keep = keep0[bcl2_df$Approved.symbol]),
  file.path(
    outdir, "non-cycling_subset.pseudobulk.raw_counts.csv"))

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
  features = bcl2_df$Approved.symbol,
  center = TRUE,
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
y <- y[keep, ,keep.lib.sizes = FALSE]
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

write.csv(
  cbind(
    as.data.frame(y$counts[bcl2_df$Approved.symbol, ]),
    keep = keep0[bcl2_df$Approved.symbol]),
  file.path(
    outdir, "cycling_subset.pseudobulk.raw_counts.csv"))
