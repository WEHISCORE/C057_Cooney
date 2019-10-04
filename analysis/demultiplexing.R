library(scater)
library(here)
source(here("analysis", "helper_functions.R"))

sce <- readRDS(here("data", "SCEs", "HTLV_GEX_HTO.CellRanger.SCE.rds"))
colnames(sce) <- sce$Barcode

empties <- readRDS(here("data", "emptyDrops", "empties.rds"))
sce <- sce[, which(empties$FDR <= 0.01)]

# Split the data up into RNA and HTOs
rna <- sce[grepl("ENS", rownames(sce)), ]
hto <- sce[!grepl("ENS", rownames(sce)), ]
# Select cell barcodes detected by both RNA and HTOs
rna_bcs <- colSums(assay(rna)) > 0
hto_bcs <- colSums(assay(hto)) > 0
joint_bcs <- rna_bcs & hto_bcs
rna <- rna[, joint_bcs]
hto <- hto[, joint_bcs]

hto <- hto[paste0("human_", 1:5), ]

sizeFactors(hto) <- librarySizeFactors(hto)
hto <- normalizeSCE(hto)

hto <- runPCA(hto, ncomponents = 5)
set.seed(11278)
hto <- runTSNE(hto)

library(scran)
g <- buildSNNGraph(hto, k = 60, use.dimred = "PCA")
clust <- igraph::cluster_louvain(g)$membership
sort(table(clust), decreasing = TRUE)
hto$cluster <- factor(clust)
cluster_colours <- setNames(
  Polychrome::dark.colors(nlevels(hto$cluster)),
  levels(hto$cluster))

.plotHeatmap(
  object = hto,
  features = rownames(hto),
  columns = order(hto$cluster),
  colour_columns_by = "cluster",
  cluster_cols = FALSE,
  center = FALSE,
  symmetric = FALSE,
  show_colnames = FALSE,
  fontsize = 8,
  annotation_colors = list(cluster = cluster_colours),
  color = viridis::inferno(100))

.plotHeatmap(
  object = hto,
  features = rownames(hto),
  columns = order(hto$cluster),
  colour_columns_by = "cluster",
  cluster_cols = FALSE,
  center = TRUE,
  symmetric = TRUE,
  show_colnames = FALSE,
  fontsize = 8,
  annotation_colors = list(cluster = cluster_colours),
  color = colorspace::divergex_hcl(n = 100, palette = "RdBu", rev = TRUE))

plotTSNE(hto, colour_by = "cluster", point_size = 0) +
  geom_point(aes(colour = colour_by), size = 0.5) +
  scale_colour_manual(values = cluster_colours, name = "cluster")

stopifnot(identical(colnames(sce), colnames(hto)))
sce$hto_cluster <- hto$cluster
saveRDS(
  sce,
  here("data", "SCEs", "C057_Cooney.demultiplexed.SCE.rds"),
  compress = "xz")
