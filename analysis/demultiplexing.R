library(scater)
library(here)
library(cowplot)
library(dplyr)
library(janitor)

source(here("analysis", "helper_functions.R"))

sce <- readRDS(here("data", "SCEs", "HTLV_GEX_HTO.CellRanger.SCE.rds"))
colnames(sce) <- sce$Barcode

empties <- readRDS(here("data", "emptyDrops", "empties.rds"))
sce <- sce[, which(empties$FDR <= 0.01)]

# Split the data up into RNA and HTOs
rna <- sce[grepl("ENS", rownames(sce)), ]
hto <- sce[!grepl("ENS", rownames(sce)), ]
counts(hto) <- as.matrix(counts(hto))

# Select cell barcodes detected by both RNA and HTOs
# NOTE: In this case, all droplets are retained because emptyDrops() has
#       already been run on the RNA and because the HTOs were sequenced quite
#       deeply.
rna_bcs <- colSums(assay(rna)) > 0
hto_bcs <- colSums(assay(hto)) > 0
joint_bcs <- rna_bcs & hto_bcs
rna <- rna[, joint_bcs]
hto <- hto[, joint_bcs]

# We know that only the human HTOs (human_1, ..., human_5) were used.
hto <- hto[paste0("human_", 1:5), ]

# Compute size factors
# (1) Library size factors
sf_lib <- librarySizeFactors(hto)
summary(sf_lib)
# (2) Geometric mean of counts
sf_geo <- exp(colMeans(log1p(counts(hto))))
sf_geo <- sf_geo / mean(sf_geo)
summary(sf_geo)

library(scran)
tagdata <- normalizeCounts(counts(hto), sf_lib, log = TRUE)
# NOTE: `d = NA` means no dimensionality reduction is performed.
g <- buildSNNGraph(tagdata, k = 100, d = NA)
clusters <- igraph::cluster_louvain(g)$membership
plot(
  x = sf_lib,
  y = sf_geo,
  log="xy",
  col = clusters,
  xlab = "Library size factors (tag)",
  ylab = "Geometric mean size factors (tag)",
  pch = 16,
  cex = 0.5)
abline(0, 1, col = "grey", lty = 2, lwd = 2)

sizeFactors(hto) <- sf_lib
hto <- normalizeSCE(hto)

set.seed(1010010)
hto <- runTSNE(hto)
set.seed(4657)
hto <- runUMAP(hto)

g <- buildSNNGraph(hto, k = 100, d = NA)
clusters <- igraph::cluster_louvain(g)$membership
hto$cluster <- factor(clusters)
cluster_colours <- setNames(
  colorblindr::palette_OkabeIto,
  levels(hto$cluster))
cluster_colours <- cluster_colours[!is.na(names(cluster_colours))]

averaged <- sweep(
  sumCountsAcrossCells(hto, clusters, exprs_values = "logcounts"),
  2,
  table(hto$cluster),
  "/")
# Sanity check
stopifnot(all(apply(averaged, 2, which.max) == c(3, 1, 2, 1, 4, 5, 1)))

library(pheatmap)
mat <- averaged - rowMeans(averaged)
rownames(mat) <- rownames(hto)
pheatmap(
  mat = mat,
  breaks = seq(-max(abs(mat)), max(abs(mat)), length.out = 101),
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  annotation_col = data.frame(
    cluster = colnames(averaged),
    row.names = colnames(averaged)),
  annotation_colors = list(cluster = cluster_colours),
  color = colorspace::divergex_hcl(n = 100, palette = "RdBu", rev = TRUE),
  main = "Row-normalized average log-expression of HTOs")

.plotHeatmap(
  object = hto,
  features = rownames(hto),
  columns = order(hto$cluster),
  colour_columns_by = "cluster",
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  center = TRUE,
  symmetric = TRUE,
  show_colnames = FALSE,
  fontsize = 8,
  annotation_colors = list(cluster = cluster_colours),
  color = colorspace::divergex_hcl(n = 100, palette = "RdBu", rev = TRUE),
  main = "Row-normalized log-expression of HTOs")

plotTSNE(hto, colour_by = "cluster", point_size = 0) +
  geom_point(aes(colour = colour_by), size = 0.5) +
  scale_colour_manual(values = cluster_colours, name = "cluster") +
  ggtitle("t-SNE: HTOs") +
  theme_cowplot()

plotUMAP(hto, colour_by = "cluster", point_size = 0) +
  geom_point(aes(colour = colour_by), size = 0.5) +
  scale_colour_manual(values = cluster_colours, name = "cluster") +
  ggtitle("UMAP: HTOs") +
  theme_cowplot()

# Map clusters to HTOs to samples
cluster_to_hto_df <- data.frame(
  cluster = levels(hto$cluster),
  HTO = c(
    paste0("human_", c(3, 1, 2, paste0(1:5, collapse = "_"), 4, 5)),
    NA_character_),
  stringsAsFactors = FALSE)
hto_to_sample_df <- data.frame(
  HTO = paste0("human_", c(1:5, paste0(1:5, collapse = "_"))),
  Sample = c(paste0("infected_", 1:3), paste0("uninfected_", 1:3)),
  Treatment = c(rep("Infected", 3), rep("Uninfected", 3)),
  Replicate = rep(1:3, 2),
  stringsAsFactors = FALSE)
hto_df <- left_join(cluster_to_hto_df, hto_to_sample_df)
knitr::kable(hto_df, caption = "Table linking HTO clusters to HTO to samples.")

stopifnot(identical(colnames(sce), colnames(hto)))
colData(sce) <- colData(sce) %>%
  as.data.frame() %>%
  select(-Sample) %>%
  mutate(cluster = hto$cluster) %>%
  left_join(hto_df) %>%
  select(-cluster) %>%
  DataFrame()

# NOTE: James was targetting ~2000 cells/sample
tabyl(
  as.data.frame(colData(sce)[, c("Treatment", "Replicate")]),
  Treatment,
  Replicate) %>%
  adorn_title("combined") %>%
  knitr::kable(
    caption = "Tabulation of samples by treatment and replicate (mouse). The samples are unpaired. 'NA' values are those samples who HTOs could not be demultiplexed.")

saveRDS(
  sce,
  here("data", "SCEs", "C057_Cooney.demultiplexed.SCE.rds"),
  compress = "xz")
