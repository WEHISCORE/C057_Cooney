library(tricycle)
library(here)
sce <- readRDS(here("data/SCEs/C057_Cooney.annotated.SCE.rds"))
# sce <- readRDS(here("data/SCEs/C057_Cooney.cycling.annotated.SCE.rds"))
# sce <- readRDS(here("data/SCEs/C057_Cooney.not_cycling.annotated.SCE.rds"))

sce <- project_cycle_space(
  x = sce,
  gname.type = "SYMBOL",
  species = "human")

library(ggplot2)
library(scattermore)
library(scater)
scater::plotReducedDim(sce, dimred = "tricycleEmbedding") +
  labs(x = "Projected PC1", y = "Projected PC2") +
  ggtitle(sprintf("Projected cell cycle space (n=%d)", ncol(sce))) +
  theme_bw(base_size = 14)

sce <- estimate_cycle_position(sce)

top2a.idx <- which(rownames(sce) == "TOP2A")
fit.l <- fit_periodic_loess(
  sce$tricyclePosition,
  assay(sce, 'logcounts')[top2a.idx, ],
  plot = TRUE,
  x_lab = "Cell cycle position \u03b8",
  y_lab = "log2(Top2a)",
  fig.title = paste0(
    "Expression of Top2a along \u03b8 (n=",
    ncol(sce),
    ")"))
fit.l

sce <- estimate_cycle_stage(
  sce,
  gname.type = "SYMBOL",
  species = "human")
table(sce$CCStage, useNA = "ifany", sce$cluster)

scater::plotReducedDim(
  sce,
  dimred = "tricycleEmbedding",
  colour_by = "CCStage") +
  labs(
    x = "Projected PC1",
    y = "Projected PC2",
    title = paste0("Projected cell cycle space (n=", ncol(sce), ")")) +
  theme_bw(base_size = 14)

scater::plotReducedDim(
  sce,
  dimred = "UMAP",
  colour_by = "CCStage",
  point_size = 0.5,
  point_alpha = 1) +
  labs(
    title = paste0("Projected cell cycle space (n=", ncol(sce), ")")) +
  theme_cowplot()

plot_ccposition_den(
  sce$tricyclePosition,
  sce$Sample,
  "Sample",
  bw = 10,
  fig.title = "Kernel density of \u03b8") +
  theme_bw(base_size = 14)

plot_ccposition_den(
  sce$tricyclePosition,
  sce$Sample,
  "Sample",
  type = "circular",
  bw = 10,
  fig.title = "Kernel density of \u03b8") +
  theme_bw(base_size = 14)

library(cowplot)
p <- plot_emb_circle_scale(
  sce,
  dimred = "UMAP",
  point.size = 1,
  point.alpha = 1) +
  theme_cowplot()
legend <- circle_scale_legend(text.size = 5, alpha = 0.9)
plot_grid(p, legend, ncol = 2, rel_widths = c(1, 0.4))

plotColData(sce, "tricyclePosition", x = "cluster")
plotColData(sce, "CCStage", x = "cluster")
sce$tmp <- sce$CCStage
sce$tmp[is.na(sce$tmp)] <- "NA"
plotColData(sce, "tmp", x = "cluster")
