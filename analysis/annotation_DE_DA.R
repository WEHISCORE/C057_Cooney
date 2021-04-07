# Setup ------------------------------------------------------------------------

# Takes off from cell selection

library(SingleCellExperiment)
library(here)
library(cowplot)

sce <- readRDS(here("data/SCEs/C057_Cooney.cells_selected.SCE.rds"))

# data frames containing co-ordinates and factors for creating reduced
# dimensionality plots.
umap_df <- cbind(
  data.frame(
    x = reducedDim(sce, "UMAP")[, 1],
    y = reducedDim(sce, "UMAP")[, 2]),
  as.data.frame(colData(sce)[, !colnames(colData(sce)) %in% c("TRA", "TRB")]))

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

# Some useful gene sets
mito_set <- rownames(sce)[which(rowData(sce)$CHR == "MT")]
ribo_set <- grep("^RP(S|L)", rownames(sce), value = TRUE)
# NOTE: A more curated approach for identifying ribosomal protein genes
#       (https://github.com/Bioconductor/OrchestratingSingleCellAnalysis-base/blob/ae201bf26e3e4fa82d9165d8abf4f4dc4b8e5a68/feature-selection.Rmd#L376-L380)
library(msigdbr)
c2_sets <- msigdbr(species = "Homo sapiens", category = "C2")
ribo_set <- union(
  ribo_set,
  c2_sets[c2_sets$gs_name == "KEGG_RIBOSOME", ]$human_gene_symbol)

library(scran)
set.seed(1000)
var_fit <- modelGeneVarByPoisson(sce, block = sce$batch)
hvg <- getTopHVGs(var_fit, var.threshold = 0)
is_mito <- hvg %in% mito_set
is_ribo <- hvg %in% ribo_set
hvg <- hvg[!(is_mito | is_ribo)]

set.seed(1010)
var_fit.sample <- modelGeneVarByPoisson(sce, block = sce$Sample)
hvg.sample <- getTopHVGs(var_fit.sample, var.threshold = 0)
is_mito <- hvg.sample %in% mito_set
is_ribo <- hvg.sample %in% ribo_set
hvg.sample <- hvg.sample[!(is_mito | is_ribo)]

# Re-process data --------------------------------------------------------------

library(scater)
# NOTE: This is a no-op if the same hvg set is used from cell selection.
set.seed(11235)
sce <- denoisePCA(sce, var_fit, subset.row = hvg)
set.seed(8875)
sce <- runUMAP(sce, dimred = "PCA")
umap_df <- cbind(
  data.frame(
    x = reducedDim(sce, "UMAP")[, 1],
    y = reducedDim(sce, "UMAP")[, 2]),
  as.data.frame(colData(sce)[, !colnames(colData(sce)) %in% c("TRA", "TRB")]))
set.seed(8111)
snn_gr <- buildSNNGraph(sce, use.dimred = "PCA")
clusters <- igraph::cluster_louvain(snn_gr)
sce$cluster <- factor(clusters$membership)
umap_df$cluster <- sce$cluster
cluster_colours <- setNames(
  Polychrome::glasbey.colors(nlevels(sce$cluster) + 1)[-1],
  levels(sce$cluster))
sce$cluster_colours <- cluster_colours[sce$cluster]

# Diagnosing batch effects -----------------------------------------------------

# Samples were collected in a single capture, but there are some strong
# differences between Venetoclax and Control samples.
library(patchwork)
# Some clusters are highly treatment specific (just another way of showing
# what is shown in the plots at the end of 'cell selection').
table(cluster = sce$cluster, treatment = sce$Treatment)
uncorrected.p1 <- plotUMAP(sce, colour_by = "Treatment", point_size = 0.5) +
  scale_colour_manual(values = treatment_colours, name = "Treatment") +
  ggtitle("Uncorrected")
uncorrected.p2 <- plotUMAP(sce, colour_by = "cluster", point_size = 0.5) +
  scale_colour_manual(values = cluster_colours, name = "cluster") +
  ggtitle("Uncorrected")
uncorrected.p1 + uncorrected.p2

plot_grid(
  ggplot(as.data.frame(colData(sce)[, c("cluster", "Sample")])) +
    geom_bar(
      aes(x = cluster, fill = Sample),
      position = position_fill(reverse = TRUE)) +
    coord_flip() +
    ylab("Frequency") +
    scale_fill_manual(values = sample_colours) +
    theme_cowplot(font_size = 8),
  ggplot(as.data.frame(colData(sce)[, c("cluster", "Treatment")])) +
    geom_bar(
      aes(x = cluster, fill = Treatment),
      position = position_fill(reverse = TRUE)) +
    coord_flip() +
    ylab("Frequency") +
    scale_fill_manual(values = treatment_colours) +
    theme_cowplot(font_size = 8),
  ggplot(as.data.frame(colData(sce)[, "cluster", drop = FALSE])) +
    geom_bar(aes(x = cluster, fill = cluster)) +
    coord_flip() +
    ylab("Number of cells") +
    scale_fill_manual(values = cluster_colours) +
    theme_cowplot(font_size = 8) +
    guides(fill = FALSE),
  align = "h",
  ncol = 3)

# Linear regression with rescaleBatches() --------------------------------------

library(batchelor)
rescaled <- rescaleBatches(sce, batch = sce$Treatment)
library(scran)
set.seed(666)
rescaled <- denoisePCA(
  rescaled,
  var_fit,
  subset.row = hvg,
  assay.type = "corrected")
snn.gr <- buildSNNGraph(rescaled, use.dimred="PCA")
clusters.resc <- igraph::cluster_louvain(snn.gr)$membership
table(cluster=clusters.resc, treatment=rescaled$batch)
rescaled <- runUMAP(rescaled, dimred="PCA")
rescaled$batch <- factor(rescaled$batch)
rescaled$cluster <- factor(clusters.resc)
rescaled.p1 <- plotUMAP(rescaled, colour_by = "batch", point_size = 0.5) +
  scale_colour_manual(values = treatment_colours, name = "Treatment") +
  ggtitle("Rescaled")
rescaled.p2 <- plotUMAP(rescaled, colour_by = "cluster", point_size = 0.5) +
  ggtitle("Rescaled")
rescaled.p1 + rescaled.p2

plot_grid(
  ggplot(as.data.frame(
    cbind(
      colData(sce)[, "Sample", drop = FALSE],
      colData(rescaled)[, "cluster", drop = FALSE]))) +
    geom_bar(
      aes(x = cluster, fill = Sample),
      position = position_fill(reverse = TRUE)) +
    coord_flip() +
    ylab("Frequency") +
    scale_fill_manual(values = sample_colours) +
    theme_cowplot(font_size = 8),
  ggplot(as.data.frame(
    cbind(
      colData(sce)[, "Treatment", drop = FALSE],
      colData(rescaled)[, "cluster", drop = FALSE]))) +
    geom_bar(
      aes(x = cluster, fill = Treatment),
      position = position_fill(reverse = TRUE)) +
    coord_flip() +
    ylab("Frequency") +
    scale_fill_manual(values = treatment_colours) +
    theme_cowplot(font_size = 8),
  ggplot(as.data.frame(colData(rescaled)[, "cluster", drop = FALSE])) +
    geom_bar(aes(x = cluster, fill = cluster)) +
    coord_flip() +
    ylab("Number of cells") +
    theme_cowplot(font_size = 8) +
    guides(fill = FALSE),
  align = "h",
  ncol = 3)

# Linear regression with regressBatches() --------------------------------------

set.seed(10001)
residuals <- regressBatches(sce, batch = sce$Treatment, d = 50,
                            subset.row=hvg, correct.all=TRUE)

snn.gr <- buildSNNGraph(residuals, use.dimred="corrected")
clusters.resid <- igraph::cluster_louvain(snn.gr)$membership
table(Cluster=clusters.resid, treatment=residuals$batch)
residuals <- runUMAP(residuals, dimred="corrected")
residuals$batch <- factor(residuals$batch)
residuals$cluster <- factor(clusters.resid)
residuals.p1 <- plotUMAP(residuals, colour_by = "batch", point_size = 0.5) +
  scale_colour_manual(values = treatment_colours, name = "Treatment") +
  ggtitle("Residuals")
residuals.p2 <- plotUMAP(residuals, colour_by = "cluster", point_size = 0.5) +
  ggtitle("Residuals")
residuals.p1 + residuals.p2

plot_grid(
  ggplot(as.data.frame(
    cbind(
      colData(sce)[, "Sample", drop = FALSE],
      colData(residuals)[, "cluster", drop = FALSE]))) +
    geom_bar(
      aes(x = cluster, fill = Sample),
      position = position_fill(reverse = TRUE)) +
    coord_flip() +
    ylab("Frequency") +
    scale_fill_manual(values = sample_colours) +
    theme_cowplot(font_size = 8),
  ggplot(as.data.frame(
    cbind(
      colData(sce)[, "Treatment", drop = FALSE],
      colData(residuals)[, "cluster", drop = FALSE]))) +
    geom_bar(
      aes(x = cluster, fill = Treatment),
      position = position_fill(reverse = TRUE)) +
    coord_flip() +
    ylab("Frequency") +
    scale_fill_manual(values = treatment_colours) +
    theme_cowplot(font_size = 8),
  ggplot(as.data.frame(colData(residuals)[, "cluster", drop = FALSE])) +
    geom_bar(aes(x = cluster, fill = cluster)) +
    coord_flip() +
    ylab("Number of cells") +
    theme_cowplot(font_size = 8) +
    guides(fill = FALSE),
  align = "h",
  ncol = 3)

# MNN by treatment -------------------------------------------------------------

set.seed(1000101001)
mnn_treatment.out <- fastMNN(
  sce,
  batch = sce$Treatment,
  d = 50,
  k = 20,
  subset.row = hvg)
snn.gr <- buildSNNGraph(mnn_treatment.out, use.dimred="corrected")
clusters.mnn_treatment <- igraph::cluster_louvain(snn.gr)$membership
table(Cluster=clusters.mnn_treatment, treatment=mnn_treatment.out$batch)
mnn_treatment.out <- runUMAP(mnn_treatment.out, dimred="corrected")
mnn_treatment.out$batch <- factor(mnn_treatment.out$batch)
mnn_treatment.out$cluster <- factor(clusters.mnn_treatment)
mnn_treatment.p1 <- plotUMAP(
  mnn_treatment.out,
  colour_by = "batch",
  point_size = 0.5) +
  scale_colour_manual(values = treatment_colours, name = "Treatment") +
  ggtitle("MNN (treatment)")
mnn_treatment.p2 <- plotUMAP(
  mnn_treatment.out,
  colour_by = "cluster",
  point_size = 0.5) +
  ggtitle("MNN (treatment)")
mnn_treatment.p1 + mnn_treatment.p2

metadata(mnn_treatment.out)$merge.info$lost.var

plot_grid(
  ggplot(as.data.frame(
    cbind(
      colData(sce)[, "Sample", drop = FALSE],
      colData(mnn_treatment.out)[, "cluster", drop = FALSE]))) +
    geom_bar(
      aes(x = cluster, fill = Sample),
      position = position_fill(reverse = TRUE)) +
    coord_flip() +
    ylab("Frequency") +
    scale_fill_manual(values = sample_colours) +
    theme_cowplot(font_size = 8),
  ggplot(as.data.frame(
    cbind(
      colData(sce)[, "Treatment", drop = FALSE],
      colData(mnn_treatment.out)[, "cluster", drop = FALSE]))) +
    geom_bar(
      aes(x = cluster, fill = Treatment),
      position = position_fill(reverse = TRUE)) +
    coord_flip() +
    ylab("Frequency") +
    scale_fill_manual(values = treatment_colours) +
    theme_cowplot(font_size = 8),
  ggplot(as.data.frame(colData(mnn_treatment.out)[, "cluster", drop = FALSE])) +
    geom_bar(aes(x = cluster, fill = cluster)) +
    coord_flip() +
    ylab("Number of cells") +
    theme_cowplot(font_size = 8) +
    guides(fill = FALSE),
  align = "h",
  ncol = 3)

# MNN by sample -----------------------------------------------------------------

set.seed(1000101001)
mnn_sample.out <- fastMNN(
  sce,
  batch = sce$Sample,
  d = 50,
  k = 20,
  subset.row = hvg.sample)
snn.gr <- buildSNNGraph(mnn_sample.out, use.dimred="corrected")
clusters.mnn_sample <- igraph::cluster_louvain(snn.gr)$membership
table(Cluster = clusters.mnn_sample, treatment = mnn_sample.out$batch)
mnn_sample.out <- runUMAP(mnn_sample.out, dimred="corrected")
mnn_sample.out$batch <- factor(mnn_sample.out$batch)
mnn_sample.out$cluster <- factor(clusters.mnn_sample)
mnn_sample.p1 <- plotUMAP(
  mnn_sample.out,
  colour_by = I(sce$Treatment),
  point_size = 0.5) +
  scale_colour_manual(values = treatment_colours, name = "Treatment") +
  ggtitle("MNN (sample)")
mnn_sample.p2 <- plotUMAP(
  mnn_sample.out,
  colour_by = "cluster",
  point_size = 0.5) +
  ggtitle("MNN (sample)")
mnn_sample.p3 <- plotUMAP(
  mnn_sample.out,
  colour_by = "batch",
  point_size = 0.5) +
  scale_colour_manual(values = sample_colours, name = "sample") +
  ggtitle("MNN (sample)")
mnn_sample.p1 + mnn_sample.p2 + mnn_sample.p3 + plot_layout(ncol = 2)

metadata(mnn_sample.out)$merge.info$lost.var

plot_grid(
  ggplot(as.data.frame(
    cbind(
      colData(sce)[, "Sample", drop = FALSE],
      colData(mnn_sample.out)[, "cluster", drop = FALSE]))) +
    geom_bar(
      aes(x = cluster, fill = Sample),
      position = position_fill(reverse = TRUE)) +
    coord_flip() +
    ylab("Frequency") +
    scale_fill_manual(values = sample_colours) +
    theme_cowplot(font_size = 8),
  ggplot(as.data.frame(
    cbind(
      colData(sce)[, "Treatment", drop = FALSE],
      colData(mnn_sample.out)[, "cluster", drop = FALSE]))) +
    geom_bar(
      aes(x = cluster, fill = Treatment),
      position = position_fill(reverse = TRUE)) +
    coord_flip() +
    ylab("Frequency") +
    scale_fill_manual(values = treatment_colours) +
    theme_cowplot(font_size = 8),
  ggplot(as.data.frame(colData(mnn_sample.out)[, "cluster", drop = FALSE])) +
    geom_bar(aes(x = cluster, fill = cluster)) +
    coord_flip() +
    ylab("Number of cells") +
    theme_cowplot(font_size = 8) +
    guides(fill = FALSE),
  align = "h",
  ncol = 3)

# Summary of data integration results ------------------------------------------

uncorrected.p1 + rescaled.p1 + residuals.p1 + mnn_treatment.p1 +
  mnn_sample.p1 + plot_layout(guides = "collect")
uncorrected.p2 + rescaled.p2 + residuals.p2 + mnn_treatment.p2 +
  mnn_sample.p2

# NOTE: Uncorrected data has strong treatment-specific differences (and
#       presumably donor-specific differences).
#       Correcting for treatment using any of the 4 approaches removes/dampens
#       those treatment-specific differences.
#       However, it's a little hard to interpret if this is a good or bad thing.
#       Therefore, will just use uncorrected data in what follows while
#       remembering that this means there will be very treatment-specific
#       (and perhaps donor-specific) clusters.

tmp <- sce[, sce$cluster %in% c("1", "4", "5", "7", "2", "10")]
m <- findMarkers(tmp, tmp$Treatment, direction = "up")

p1 <- plotUMAP(
  sce,
  colour_by = "Treatment",
  text_by = "cluster",
  point_size = 0.5,
  point_alpha = 1) +
  scale_colour_manual(values = treatment_colours, name = "Treatment")
p2 <- plotUMAP(
  sce,
  colour_by = rownames(m$Infected)[1],
  text_by = "cluster",
  point_size = 0.5,
  point_alpha = 1) +
  scale_color_viridis_c(option = "magma", name = rownames(m$Infected)[1])
p3 <- plotUMAP(
  sce,
  colour_by = rownames(m$Infected)[2],
  text_by = "cluster",
  point_size = 0.5,
  point_alpha = 1) +
  scale_color_viridis_c(option = "magma", name = rownames(m$Infected)[2])
p4 <- plotUMAP(
  sce,
  colour_by = rownames(m$Infected)[3],
  text_by = "cluster",
  point_size = 0.5,
  point_alpha = 1) +
  scale_color_viridis_c(option = "magma", name = rownames(m$Infected)[3])
p5 <- plotUMAP(
  sce,
  colour_by = rownames(m$Infected)[4],
  text_by = "cluster",
  point_size = 0.5,
  point_alpha = 1) +
  scale_color_viridis_c(option = "magma", name = rownames(m$Infected)[4])
p6 <- plotUMAP(
  sce,
  colour_by = rownames(m$Uninfected)[1],
  text_by = "cluster",
  point_size = 0.5,
  point_alpha = 1) +
  scale_color_viridis_c(option = "magma", name = rownames(m$Uninfected)[1])
p7 <- plotUMAP(
  sce,
  colour_by = rownames(m$Uninfected)[2],
  text_by = "cluster",
  point_size = 0.5,
  point_alpha = 1) +
  scale_color_viridis_c(option = "magma", name = rownames(m$Uninfected)[2])
p8 <- plotUMAP(
  sce,
  colour_by = rownames(m$Uninfected)[3],
  text_by = "cluster",
  point_size = 0.5,
  point_alpha = 1) +
  scale_color_viridis_c(option = "magma", name = rownames(m$Uninfected)[3])
p9 <- plotUMAP(
  sce,
  colour_by = rownames(m$Uninfected)[4],
  text_by = "cluster",
  point_size = 0.5,
  point_alpha = 1) +
  scale_color_viridis_c(option = "magma", name = rownames(m$Uninfected)[4])
p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + plot_layout(ncol = 3)

plotUMAP(
  sce,
  colour_by = "GNLY",
  point_size = 0.5,
  point_alpha = 1,
  other_fields = "Sample") +
  scale_color_viridis_c(option = "magma", name = "GNLY") +
  facet_wrap(~Sample, ncol = 3)

# Clustering at different resolutions ------------------------------------------

set.seed(126)
snn_gr <- buildSNNGraph(sce, use.dimred = "PCA", k = 20)
clusters <- igraph::cluster_louvain(snn_gr)
sce$cluster_k20 <- factor(clusters$membership)
snn_gr <- buildSNNGraph(sce, use.dimred = "PCA", k = 50)
clusters <- igraph::cluster_louvain(snn_gr)
sce$cluster_k50 <- factor(clusters$membership)
snn_gr <- buildSNNGraph(sce, use.dimred = "PCA", k = 100)
clusters <- igraph::cluster_louvain(snn_gr)
sce$cluster_k100 <- factor(clusters$membership)

p1 <- plotUMAP(
  sce,
  colour_by = "cluster",
  point_size = 0.5,
  point_alpha = 1,
  text_by = "cluster") +
  scale_colour_manual(values = cluster_colours, name = "cluster") +
  ggtitle("Uncorrected: k = 10") +
  guides(colour = FALSE)
p2 <- plotUMAP(
  sce,
  colour_by = "cluster_k20",
  point_size = 0.5,
  point_alpha = 1,
  text_by = "cluster_k20") +
  ggtitle("Uncorrected: k = 20") +
  guides(colour = FALSE)
p3 <- plotUMAP(
  sce,
  colour_by = "cluster_k50",
  point_size = 0.5,
  point_alpha = 1,
  text_by = "cluster_k50") +
  ggtitle("Uncorrected: k = 50") +
  guides(colour = FALSE)
p4 <- plotUMAP(
  sce,
  colour_by = "cluster_k100",
  point_size = 0.5,
  point_alpha = 1,
  text_by = "cluster_k100") +
  ggtitle("Uncorrected: k = 100") +
  guides(colour = FALSE)
p1 + p2 + p3 + p4 + plot_layout(ncol = 2)

# NOTE: In all clusterings we observe an 'out group' on the right hand side.

# Some clusters are readily annotated by uniquely upregulated markers ----------

plotUMAP(
  sce,
  colour_by = "cluster",
  point_size = 0.5,
  point_alpha = 1,
  text_by = "cluster") +
  scale_colour_manual(values = cluster_colours, name = "cluster") +
  ggtitle("Uncorrected") +
  guides(colour = guide_legend(override.aes = list(size = 5)))

out <- pairwiseTTests(sce, sce$cluster, direction = "up", lfc = 0.5)
# TODO: Block?
top_markers <- getTopMarkers(
  out$statistics,
  out$pairs,
  pairwise = FALSE,
  pval.type = "all")
features <- unlist(top_markers)
sce$libsize <- log10(sce$sum)
plotHeatmap(
  sce,
  features,
  order_columns_by = c("cluster", "libsize"),
  cluster_rows = FALSE,
  color = hcl.colors(101, "Blue-Red 3"),
  center = TRUE,
  zlim = c(-3, 3),
  column_annotation_colors = list(
    cluster = cluster_colours),
  annotation_row = data.frame(
    cluster = names(features),
    row.names = features),
  main = "Uniquely upregulated (logFC > 0.5) cluster marker genes")

# Look at cycling cluster(s) ---------------------------------------------------

# NOTE: Cell cycle could really be proliferation of T-cells upon activation

# Appears to be a cycling and non-cycling clustering driving the results.
# Cycling cluster(s) typified by upregulation of a large number of genes, larger
# library sizes, and cyclin expression.
# Non-cycling clusters typified by opposite and no upregulation of specific
# markers.

markers <- findMarkers(
  sce,
  sce$cluster,
  direction = "up",
  block = sce$Sample,
  lfc = 0.5,
  pval.type = "all")
cluster_8_markers <- markers[["8"]][markers[["8"]]$FDR < 0.05, ]
features <- rownames(cluster_8_markers)
rowData(sce)[features, ]

library(limma)
go <- goana(
  de = unlist(rowData(sce)[features, "NCBI.ENTREZID"]),
  species = "Hs")
topGO(go)

cycling_subset <- "8"
sce$cycling_subset <- ifelse(
  sce$cluster %in% cycling_subset,
  "Cycling",
  "Not cycling")
cycling_subset_markers <- findMarkers(
  sce,
  sce$cycling_subset,
  direction = "up",
  block = sce$Sample,
  lfc = 0.5)


p1 <- plotUMAP(
  sce,
  colour_by = "cluster",
  text_by = "cluster",
  point_size = 0.5,
  point_alpha = 1) +
  scale_colour_manual(values = cluster_colours) +
  guides(colour = FALSE)
p2 <- plotUMAP(
  sce,
  colour_by = "cycling_subset",
  text_by = "cluster",
  point_size = 0.5,
  point_alpha = 1)
p3 <- plotColData(
  sce,
  "sum",
  x = "cluster",
  colour_by = "cycling_subset",
  point_size = 0.5,
  point_alpha = 1)
(p1 + p2) / p3 + plot_layout(guides = "collect")

cyclin_genes <- grep("^CCN[ABDE][0-9]$", rowData(sce)$Symbol)
cyclin_genes <- sort(rownames(sce)[cyclin_genes])
cyclin_genes

features <- unlist(
  List(
    cycling_markers = head(
      rownames(cycling_subset_markers[["Not cycling"]]),
      25),
    not_cycling_markers = head(
      rownames(cycling_subset_markers[["Cycling"]]),
      25),
    cyclins = cyclin_genes))
plotHeatmap(
  sce,
  features,
  order_columns_by = c("cycling_subset", "cluster", "libsize"),
  cluster_rows = FALSE,
  color = hcl.colors(101, "Blue-Red 3"),
  center = TRUE,
  zlim = c(-3, 3),
  column_annotation_colors = list(cluster = cluster_colours),
  annotation_row = data.frame(
    GeneSet = names(features),
    row.names = features))

plot_grid(
  ggplot(as.data.frame(colData(sce)[, c("cycling_subset", "Sample")])) +
    geom_bar(
      aes(x = cycling_subset, fill = Sample),
      position = position_fill(reverse = TRUE)) +
    coord_flip() +
    ylab("Frequency") +
    scale_fill_manual(values = sample_colours) +
    theme_cowplot(font_size = 8),
  ggplot(as.data.frame(colData(sce)[, c("cycling_subset", "Treatment")])) +
    geom_bar(
      aes(x = cycling_subset, fill = Treatment),
      position = position_fill(reverse = TRUE)) +
    coord_flip() +
    ylab("Frequency") +
    scale_fill_manual(values = treatment_colours) +
    theme_cowplot(font_size = 8),
  ggplot(as.data.frame(colData(sce)[, "cycling_subset", drop = FALSE])) +
    geom_bar(aes(x = cycling_subset, fill = cycling_subset)) +
    coord_flip() +
    ylab("Number of cells") +
    theme_cowplot(font_size = 8) +
    guides(fill = FALSE),
  align = "h",
  ncol = 3)

# Using the cyclins ------------------------------------------------------------

# - Cyclin A is expressed across S and G2
# - Cyclin B is expressed highest in late G2 and mitosis (Morgan 2007).
# - Cyclin D is expressed throughout but peaks at G1
# - Cyclin E is expressed highest in the G1/S transition

plotHeatmap(
  sce,
  order_columns_by = c("cycling_subset", "cluster"),
  cluster_rows = FALSE,
  features = sort(cyclin_genes),
  color = viridisLite::magma(101),
  column_annotation_colors = list(cluster = cluster_colours))

markers <- findMarkers(
  sce,
  sce$cluster,
  subset.row = cyclin_genes,
  test.type = "wilcox",
  direction = "up")
markers[[cycling_subset]]

# Using the cyclone() classifier -----------------------------------------------

# hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds",
#                                 package="scran"))
#
# # Using Ensembl IDs to match up with the annotation in 'hs.pairs'.
# NOTE: This was run and the results saved.
# set.seed(100)
# assignments <- cyclone(
#   sce,
#   hs.pairs,
#   gene.names = rowData(sce)$ENSEMBL.GENEID)
assignments <- readRDS(here("tmp/cyclone_assignments.rds"))
colData(sce) <- cbind(colData(sce), DataFrame(assignments$score))
sce$phases <- assignments$phases

plot(
  x = sce$G1,
  y = sce$G2M,
  xlab = "G1 score",
  ylab = "G2/M score",
  pch = 16,
  col = sce$cluster_colours)

par(mfrow = c(2, 2))
boxplot(
  sce$G1 ~ sce$cluster,
  ylab = "score",
  col = cluster_colours,
  xlab = "cluster",
  main = "Cyclone: G1")
boxplot(
  sce$G2M ~ sce$cluster,
  ylab = "score",
  col = cluster_colours,
  xlab = "cluster",
  main = "Cyclone: G2M")
boxplot(
  sce$S ~ sce$cluster,
  ylab = "score",
  col = cluster_colours,
  xlab = "cluster",
  main = "Cyclone: S")

proportions(table(sce$cluster, assignments$phases), 1)
proportions(table(sce$cycling_subset, sce$phases), 1)

plotUMAP(
  sce,
  colour_by = "G1",
  text_by = "cluster",
  point_size = 0.5,
  point_alpha = 1) +
  plotUMAP(
    sce,
    colour_by = "G2M",
    text_by = "cluster",
    point_size = 0.5,
    point_alpha = 1) +
  plotUMAP(
    sce,
    colour_by = "S",
    text_by = "cluster",
    point_size = 0.5,
    point_alpha = 1) +
  plotUMAP(
    sce,
    colour_by = "cluster",
    text_by = "cluster",
    point_size = 0.5,
    point_alpha = 1) +
  scale_colour_manual(values = cluster_colours, name = "cluster") +
  plot_layout(ncol = 2)

# NOTE: Not sure how much I trust these results given the classifier was trained
#       on quite a different dataset.

# Subset analysis based on cycling cluster(s) ----------------------------------

set.seed(9391)
list_of_sce <- quickSubCluster(
  sce,
  groups = sce$cycling_subset,
  prepFUN = function(x) {
    var_fit <- modelGeneVarByPoisson(x)
    hvg <- getTopHVGs(var_fit, var.threshold = 0)
    is_mito <- hvg %in% mito_set
    is_ribo <- hvg %in% ribo_set
    hvg <- hvg[!(is_mito | is_ribo)]
    # NOTE: Keep the original dimensionality reduction around for downstream
    #       plotting.
    reducedDimNames(x) <- paste0("original_", reducedDimNames(x))
    denoisePCA(x, var_fit, subset.row = hvg)
  },
  clusterFUN = function(x) {
    snn_gr <- buildSNNGraph(x, use.dimred = "PCA", k = 50)
    factor(igraph::cluster_louvain(snn_gr)$membership)
  })
# It's also useful to have per-sample UMAP representations.
set.seed(17127)
list_of_sce <- lapply(list_of_sce, runUMAP, dimred = "PCA")
# It's also useful to have SingleR DICE with fine labels for subclusters
library(SingleR)
library(celldex)
ref <- DatabaseImmuneCellExpressionData()
labels_fine <- ref$label.fine
# NOTE: This code doesn't necessarily generalise beyond the DICE main labels.
label_fine_collapsed_colours <- setNames(
  c(
    Polychrome::glasbey.colors(nlevels(factor(labels_fine)) + 1)[-1],
    "orange"),
  c(levels(factor(labels_fine)), "other"))
list_of_sce <- lapply(list_of_sce, function(x) {
  pred_subcluster_fine <- SingleR(
    test = x,
    ref = ref,
    labels = labels_fine,
    clusters = x$subcluster)
  x$label_subcluster_fine <- factor(
    pred_subcluster_fine[x$subcluster, "pruned.labels"])
  x
})
list_of_sce <- lapply(list_of_sce, function(x) {
  pred_cell_fine <- SingleR(
    test = x,
    ref = ref,
    labels = labels_fine)
  x$label_cell_fine <- factor(pred_cell_fine$pruned.labels)
  x
})

plot_grid(
  plotlist = lapply(names(list_of_sce), function(n) {
    x <- list_of_sce[[n]]
    plotUMAP(x, colour_by = "subcluster", point_size = 0.5, point_alpha = 1) +
      ggtitle(n)
  }))

plot_grid(
  plotlist = lapply(names(list_of_sce), function(n) {
    x <- list_of_sce[[n]]
    plotUMAP(
      x,
      colour_by = "label_subcluster_fine",
      text_by = "subcluster",
      point_size = 0.5,
      point_alpha = 1) +
      ggtitle(n)
  }))

plot_grid(
  plotlist = lapply(names(list_of_sce), function(n) {
    x <- list_of_sce[[n]]
    plotUMAP(
      x,
      colour_by = "label_cell_fine",
      text_by = "subcluster",
      point_size = 0.5,
      point_alpha = 1) +
      ggtitle(n)
  }))

# Cycling subset analysis-------------------------------------------------------

cycling_sce <- list_of_sce[["Cycling"]]

plot_grid(
  plotlist = c(
    lapply(c("G1", "G2M", "S"), function(phase) {
      plotUMAP(
        cycling_sce,
        colour_by = phase,
        text_by = "subcluster",
        point_size = 0.5,
        point_alpha = 1)
    }),
    list(
      ggplot(
        as.data.frame(
          colData(cycling_sce)[, c("subcluster", "phases")])) +
        geom_bar(
          aes(x = subcluster, fill = phases),
          position = position_fill(reverse = TRUE)) +
        coord_flip() +
        ylab("Frequency") +
        theme_cowplot(font_size = 8)))) +
  plot_annotation(title = "Cycling cluster")

out <- pairwiseTTests(
  cycling_sce,
  cycling_sce$subcluster,
  direction = "up",
  block = cycling_sce$Sample)
top_markers <- getTopMarkers(
  out$statistics,
  out$pairs,
  pairwise = FALSE,
  pval.type = "all",
  n = 10)
features <- unlist(top_markers)
markers <- findMarkers(
  cycling_sce,
  cycling_sce$subcluster,
  direction = "up",
  pval.type = "all",
  block = cycling_sce$Sample)
plotHeatmap(
  cycling_sce,
  features,
  order_columns_by = c("subcluster", "Treatment", "Sample"),
  cluster_rows = FALSE,
  color = hcl.colors(101, "Blue-Red 3"),
  center = TRUE,
  zlim = c(-3, 3),
  column_annotation_colors = list(
    Treatment = treatment_colours,
    Sample = sample_colours),
  annotation_row = data.frame(
    subcluster = names(features),
    row.names = features),
  main = "Cycling: Uniquely upregulated cluster marker genes")

plot_grid(
  ggplot(as.data.frame(colData(cycling_sce)[, c("subcluster", "Sample")])) +
    geom_bar(
      aes(x = subcluster, fill = Sample),
      position = position_fill(reverse = TRUE)) +
    coord_flip() +
    ylab("Frequency") +
    scale_fill_manual(values = sample_colours) +
    theme_cowplot(),
  ggplot(as.data.frame(colData(cycling_sce)[, c("subcluster", "Treatment")])) +
    geom_bar(
      aes(x = subcluster, fill = Treatment),
      position = position_fill(reverse = TRUE)) +
    coord_flip() +
    ylab("Frequency") +
    scale_fill_manual(values = treatment_colours) +
    theme_cowplot(),
  ggplot(as.data.frame(colData(cycling_sce)[, "subcluster", drop = FALSE])) +
    geom_bar(aes(x = subcluster, fill = subcluster)) +
    coord_flip() +
    ylab("Number of cells") +
    theme_cowplot() +
    guides(fill = FALSE),
  align = "h",
  ncol = 3) +
  plot_annotation(title = "Cycling")

plotUMAP(
  cycling_sce,
  colour_by = "subcluster",
  text_by = "subcluster",
  point_size = 0.5,
  point_alpha = 1) +
  plotUMAP(
    cycling_sce,
    colour_by = "Treatment",
    text_by = "subcluster",
    point_size = 0.5,
    point_alpha = 1) +
  scale_colour_manual(values = treatment_colours) +
  plot_annotation("Cycling subset")

# Not cycling subset analysis --------------------------------------------------

not_cycling_sce <- list_of_sce[["Not cycling"]]

plot_grid(
  plotlist = c(
    lapply(c("G1", "G2M", "S"), function(phase) {
      plotUMAP(
        not_cycling_sce,
        colour_by = phase,
        text_by = "subcluster",
        point_size = 0.5,
        point_alpha = 1)
    }),
    list(
      ggplot(
        as.data.frame(
          colData(not_cycling_sce)[, c("subcluster", "phases")])) +
        geom_bar(
          aes(x = subcluster, fill = phases),
          position = position_fill(reverse = TRUE)) +
        coord_flip() +
        ylab("Frequency") +
        theme_cowplot(font_size = 8)))) +
  plot_annotation(title = "Not cycling subset")

out <- pairwiseTTests(
  not_cycling_sce,
  not_cycling_sce$subcluster,
  direction = "up",
  block = not_cycling_sce$Sample)
top_markers <- getTopMarkers(
  out$statistics,
  out$pairs,
  pairwise = FALSE,
  pval.type = "all",
  n = 10)
features <- unlist(top_markers)
markers <- findMarkers(
  not_cycling_sce,
  not_cycling_sce$subcluster,
  direction = "up",
  pval.type = "all",
  block = not_cycling_sce$Sample)
plotHeatmap(
  not_cycling_sce,
  features,
  order_columns_by = c("subcluster", "Treatment", "Sample"),
  cluster_rows = FALSE,
  color = hcl.colors(101, "Blue-Red 3"),
  center = TRUE,
  zlim = c(-3, 3),
  column_annotation_colors = list(
    Treatment = treatment_colours,
    Sample = sample_colours),
  annotation_row = data.frame(
    subcluster = names(features),
    row.names = features),
  main = "Not cycling: Uniquely upregulated cluster marker genes")

plot_grid(
  ggplot(as.data.frame(colData(not_cycling_sce)[, c("subcluster", "Sample")])) +
    geom_bar(
      aes(x = subcluster, fill = Sample),
      position = position_fill(reverse = TRUE)) +
    coord_flip() +
    ylab("Frequency") +
    scale_fill_manual(values = sample_colours) +
    theme_cowplot(),
  ggplot(as.data.frame(colData(not_cycling_sce)[, c("subcluster", "Treatment")])) +
    geom_bar(
      aes(x = subcluster, fill = Treatment),
      position = position_fill(reverse = TRUE)) +
    coord_flip() +
    ylab("Frequency") +
    scale_fill_manual(values = treatment_colours) +
    theme_cowplot(),
  ggplot(as.data.frame(colData(not_cycling_sce)[, "subcluster", drop = FALSE])) +
    geom_bar(aes(x = subcluster, fill = subcluster)) +
    coord_flip() +
    ylab("Number of cells") +
    theme_cowplot() +
    guides(fill = FALSE),
  align = "h",
  ncol = 3) +
  plot_annotation(title = "Not cycling")

plotUMAP(
  not_cycling_sce,
  colour_by = "subcluster",
  text_by = "subcluster",
  point_size = 0.5,
  point_alpha = 1) +
  plotUMAP(
    not_cycling_sce,
    colour_by = "Treatment",
    text_by = "subcluster",
    point_size = 0.5,
    point_alpha = 1) +
  scale_colour_manual(values = treatment_colours) +
  plot_annotation("Not cycling subset")

# DE analysis using cycling vs. non-cycling labels -----------------------------

library(edgeR)
plotMDS <- scater::plotMDS

# NOTE: Have to drop the complicated TRA and TRB DFrameList objects.
tmp <- sce
tmp$TRA <- NULL
tmp$TRB <- NULL
summed <- aggregateAcrossCells(
  tmp,
  id = colData(tmp)[, c("cycling_subset", "Sample")])
sizeFactors(summed) <- NULL
summed <- logNormCounts(summed)

summed <- runMDS(summed)
p1 <- plotMDS(
  summed,
  colour_by = "Treatment",
  shape_by = "cycling_subset",
  point_size = 3,
  point_alpha = 2) +
  scale_colour_manual(values = treatment_colours, name = "Treatment")
p2 <- plotMDS(
  summed,
  colour_by = "Sample",
  shape_by = "cycling_subset",
  point_size = 3,
  point_alpha = 2) +
  scale_colour_manual(values = sample_colours, name = "Treatment")
p1 + p2 + plot_layout(guides = "collect")

summed.filt <- summed[, summed$ncells >= 10]

de.results <- pseudoBulkDGE(
  summed.filt,
  label = summed.filt$cycling_subset,
  design = ~Treatment,
  coef = "TreatmentUninfected",
  condition = summed.filt$Treatment,
  include.intermediates = TRUE)

cur.results <- de.results[["Not cycling"]]
cur.results[order(cur.results$PValue),]

sizeFactors(summed.filt) <- NULL
plotExpression(
  logNormCounts(summed.filt),
  features = "IL4R",
  x = "Treatment",
  colour_by = "Treatment",
  other_fields = "cycling_subset") +
  facet_wrap(~cycling_subset) +
  scale_colour_manual(values = treatment_colours, name = "Treatment")

y.cur.results <- metadata(cur.results)$y
plotBCV(y.cur.results)

metadata(de.results)$failed

is.de <- decideTestsPerLabel(de.results, threshold=0.05)
summarizeTestsPerLabel(is.de)

# Upregulated across most labels
up.de <- is.de > 0 & !is.na(is.de)
head(sort(rowMeans(up.de), decreasing=TRUE), 10)

# Downregulated across most labels
down.de <- is.de < 0 & !is.na(is.de)
head(sort(rowMeans(down.de), decreasing=TRUE), 10)

de.specific <- pseudoBulkSpecific(
  summed.filt,
  label = summed.filt$cycling_subset,
  design = ~Treatment,
  coef = "TreatmentUninfected",
  condition = summed.filt$Treatment)

cur.specific <- de.specific[["Not cycling"]]
cur.specific <- cur.specific[order(cur.specific$PValue),]
cur.specific

# DE analysis using cluster labels ---------------------------------------------

library(edgeR)
plotMDS <- scater::plotMDS

# NOTE: Have to drop the complicated TRA and TRB DFrameList objects.
tmp <- sce
tmp$TRA <- NULL
tmp$TRB <- NULL
summed <- aggregateAcrossCells(
  tmp,
  id = colData(tmp)[, c("cluster", "Sample")])
sizeFactors(summed) <- NULL
summed <- logNormCounts(summed)

summed <- runMDS(summed)
p1 <- plotMDS(
  summed,
  colour_by = "cluster",
  shape_by = "Treatment",
  point_size = 3,
  point_alpha = 2) +
  scale_colour_manual(values = cluster_colours, name = "cluster")
p2 <- plotMDS(
  summed,
  colour_by = "cluster",
  shape_by = "Sample",
  point_size = 3,
  point_alpha = 2) +
  scale_colour_manual(values = cluster_colours, name = "cluster")
p1 + p2 + plot_layout(guides = "collect")

summed.filt <- summed[,summed$ncells >= 10]

de.results <- pseudoBulkDGE(
  summed.filt,
  label = summed.filt$cluster,
  design = ~Treatment,
  coef = "TreatmentUninfected",
  condition = summed.filt$Treatment,
  include.intermediates = TRUE)

cur.results <- de.results[["8"]]
cur.results[order(cur.results$PValue),]

sizeFactors(summed.filt) <- NULL
plotExpression(
  logNormCounts(summed.filt),
  features = "CADM1",
  x = "Treatment",
  colour_by = "Treatment",
  other_fields = "cluster") +
  facet_wrap(~cluster) +
  scale_colour_manual(values = treatment_colours, name = "Treatment")

y.cur.results <- metadata(cur.results)$y
plotBCV(y.cur.results)

metadata(de.results)$failed

is.de <- decideTestsPerLabel(de.results, threshold=0.05)
summarizeTestsPerLabel(is.de)

# Upregulated across most labels
up.de <- is.de > 0 & !is.na(is.de)
head(sort(rowMeans(up.de), decreasing=TRUE), 10)

# Downregulated across most labels
down.de <- is.de < 0 & !is.na(is.de)
head(sort(rowMeans(down.de), decreasing=TRUE), 10)

de.specific <- pseudoBulkSpecific(
  summed.filt,
  label = summed.filt$cluster,
  design = ~Treatment,
  coef = "TreatmentUninfected",
  condition = summed.filt$Treatment)

cur.specific <- de.specific[["8"]]
cur.specific <- cur.specific[order(cur.specific$PValue),]
cur.specific

# DE analysis using subcluster labels in non-cycling clusters ------------------

library(edgeR)
plotMDS <- scater::plotMDS

# NOTE: Have to drop the complicated TRA and TRB DFrameList objects.
tmp <- not_cycling_sce
tmp$TRA <- NULL
tmp$TRB <- NULL
summed <- aggregateAcrossCells(
  tmp,
  id = colData(tmp)[, c("subcluster", "Sample")])
sizeFactors(summed) <- NULL
summed <- logNormCounts(summed)

summed <- runMDS(summed)
p1 <- plotMDS(
  summed,
  colour_by = "subcluster",
  shape_by = "Treatment",
  point_size = 3,
  point_alpha = 2)
p2 <- plotMDS(
  summed,
  colour_by = "subcluster",
  shape_by = "Sample",
  point_size = 3,
  point_alpha = 2)
p1 + p2 + plot_layout(guides = "collect")

summed.filt <- summed[,summed$ncells >= 10]

de.results <- pseudoBulkDGE(
  summed.filt,
  label = summed.filt$subcluster,
  design = ~Treatment,
  coef = "TreatmentUninfected",
  condition = summed.filt$Treatment,
  include.intermediates = TRUE)

cur.results <- de.results[["Not cycling.4"]]
cur.results[order(cur.results$PValue),]

sizeFactors(summed.filt) <- NULL
plotExpression(
  logNormCounts(summed.filt),
  features = "CD74",
  x = "Treatment",
  colour_by = "Treatment",
  other_fields = "subcluster") +
  facet_wrap(~subcluster) +
  scale_colour_manual(values = treatment_colours, name = "Treatment")

y.cur.results <- metadata(cur.results)$y
plotBCV(y.cur.results)

metadata(de.results)$failed

is.de <- decideTestsPerLabel(de.results, threshold=0.05)
summarizeTestsPerLabel(is.de)

# Upregulated across most labels
up.de <- is.de > 0 & !is.na(is.de)
head(sort(rowMeans(up.de), decreasing=TRUE), 10)

# Downregulated across most labels
down.de <- is.de < 0 & !is.na(is.de)
head(sort(rowMeans(down.de), decreasing=TRUE), 10)

de.specific <- pseudoBulkSpecific(
  summed.filt,
  label = summed.filt$subcluster,
  design = ~Treatment,
  coef = "TreatmentUninfected",
  condition = summed.filt$Treatment)

cur.specific <- de.specific[["Not cycling.4"]]
cur.specific <- cur.specific[order(cur.specific$PValue),]
cur.specific

# DE analysis using subcluster labels in cycling clusters ----------------------

library(edgeR)
plotMDS <- scater::plotMDS

# NOTE: Have to drop the complicated TRA and TRB DFrameList objects.
tmp <- cycling_sce
tmp$TRA <- NULL
tmp$TRB <- NULL
summed <- aggregateAcrossCells(
  tmp,
  id = colData(tmp)[, c("subcluster", "Sample")])
sizeFactors(summed) <- NULL
summed <- logNormCounts(summed)

summed <- runMDS(summed)
p1 <- plotMDS(
  summed,
  colour_by = "subcluster",
  shape_by = "Treatment",
  point_size = 3,
  point_alpha = 2)
p2 <- plotMDS(
  summed,
  colour_by = "subcluster",
  shape_by = "Sample",
  point_size = 3,
  point_alpha = 2)
p1 + p2 + plot_layout(guides = "collect")

summed.filt <- summed[,summed$ncells >= 10]

de.results <- pseudoBulkDGE(
  summed.filt,
  label = summed.filt$subcluster,
  design = ~Treatment,
  coef = "TreatmentUninfected",
  condition = summed.filt$Treatment,
  include.intermediates = TRUE)

cur.results <- de.results[["Cycling.1"]]
cur.results[order(cur.results$PValue),]

sizeFactors(summed.filt) <- NULL
plotExpression(
  logNormCounts(summed.filt),
  features = "GLUL",
  x = "Treatment",
  colour_by = "Treatment",
  other_fields = "subcluster") +
  facet_wrap(~subcluster) +
  scale_colour_manual(values = treatment_colours, name = "Treatment")

y.cur.results <- metadata(cur.results)$y
plotBCV(y.cur.results)

metadata(de.results)$failed

is.de <- decideTestsPerLabel(de.results, threshold=0.05)
summarizeTestsPerLabel(is.de)

# Upregulated across most labels
up.de <- is.de > 0 & !is.na(is.de)
head(sort(rowMeans(up.de), decreasing=TRUE), 10)

# Downregulated across most labels
down.de <- is.de < 0 & !is.na(is.de)
head(sort(rowMeans(down.de), decreasing=TRUE), 10)

de.specific <- pseudoBulkSpecific(
  summed.filt,
  label = summed.filt$subcluster,
  design = ~Treatment,
  coef = "TreatmentUninfected",
  condition = summed.filt$Treatment)

cur.specific <- de.specific[["Cycling.1"]]
cur.specific <- cur.specific[order(cur.specific$PValue),]
cur.specific

# DA analysis using cycling vs. non-cycling labels -----------------------------

# NOTE: Can't use quasi-likelihood framework with only two labels.

abundances <- table(sce$cycling_subset, sce$Sample)
abundances <- unclass(abundances)
abundances

extra.info <- colData(sce)[match(colnames(abundances), sce$Sample),]
y.ab <- DGEList(abundances, samples=extra.info)

keep <- filterByExpr(y.ab, group=y.ab$samples$Treatment)
summary(keep)
y.ab <- y.ab[keep,]

design <- model.matrix(~Treatment, y.ab$samples)

y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion)
plotBCV(y.ab, cex=1)

fit.ab <- glmFit(y.ab)
res <- glmLRT(fit.ab, coef = "TreatmentUninfected")
summary(decideTests(res))

topTags(res)

p1 <- plotUMAP(sce, colour_by = "cycling_subset", text_by = "cycling_subset")
p2 <- plotUMAP(sce, colour_by = "Treatment") +
  scale_colour_manual(values = treatment_colours)
p1 + p2
p3 <- ggplot(as.data.frame(colData(sce)[, c("cycling_subset", "Sample")])) +
    geom_bar(
      aes(x = cycling_subset, fill = Sample),
      position = position_fill(reverse = TRUE)) +
    coord_flip() +
    ylab("Frequency") +
    scale_fill_manual(values = sample_colours) +
    theme_cowplot()
p4 <- ggplot(as.data.frame(colData(sce)[, c("cycling_subset", "Treatment")])) +
    geom_bar(
      aes(x = cycling_subset, fill = Treatment),
      position = position_fill(reverse = TRUE)) +
    coord_flip() +
    ylab("Frequency") +
    scale_fill_manual(values = treatment_colours) +
    theme_cowplot()
p5 <- ggplot(as.data.frame(colData(sce)[, "cycling_subset", drop = FALSE])) +
    geom_bar(aes(x = cycling_subset, fill = cycling_subset)) +
    coord_flip() +
    ylab("Number of cells") +
    theme_cowplot() +
    guides(fill = FALSE)
p3 + p4 + p5 + plot_layout(ncol = 3)

# DA analysis using cluster labels ---------------------------------------------

abundances <- table(sce$cluster, sce$Sample)
abundances <- unclass(abundances)
abundances

extra.info <- colData(sce)[match(colnames(abundances), sce$Sample),]
y.ab <- DGEList(abundances, samples=extra.info)

keep <- filterByExpr(y.ab, group=y.ab$samples$Treatment)
summary(keep)
y.ab <- y.ab[keep,]

design <- model.matrix(~Treatment, y.ab$samples)

y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion)
plotBCV(y.ab, cex=1)

fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)

# NOTE: With only 2 labels this doesn't work / make sense
summary(fit.ab$df.prior)
plotQLDisp(fit.ab, cex=1)

res <- glmQLFTest(fit.ab, coef = "TreatmentUninfected")
summary(decideTests(res))

topTags(res)

p1 <- plotUMAP(sce, colour_by = "cluster", text_by = "cluster") +
  scale_colour_manual(values = cluster_colours)
p2 <- plotUMAP(sce, colour_by = "Treatment") +
  scale_colour_manual(values = treatment_colours)
p1 + p2
p3 <- ggplot(as.data.frame(colData(sce)[, c("cluster", "Sample")])) +
  geom_bar(
    aes(x = cluster, fill = Sample),
    position = position_fill(reverse = TRUE)) +
  coord_flip() +
  ylab("Frequency") +
  scale_fill_manual(values = sample_colours) +
  theme_cowplot()
p4 <- ggplot(as.data.frame(colData(sce)[, c("cluster", "Treatment")])) +
  geom_bar(
    aes(x = cluster, fill = Treatment),
    position = position_fill(reverse = TRUE)) +
  coord_flip() +
  ylab("Frequency") +
  scale_fill_manual(values = treatment_colours) +
  theme_cowplot()
p5 <- ggplot(as.data.frame(colData(sce)[, "cluster", drop = FALSE])) +
  geom_bar(aes(x = cluster, fill = cluster)) +
  coord_flip() +
  ylab("Number of cells") +
  theme_cowplot() +
  guides(fill = FALSE) +
  scale_fill_manual(values = cluster_colours, name = "cluster")
p3 + p4 + p5 + plot_layout(ncol = 3)

# DA analysis using subcluster labels in non-cycling clusters ------------------

tmp <- not_cycling_sce
abundances <- table(tmp$subcluster, tmp$Sample)
abundances <- unclass(abundances)
abundances

extra.info <- colData(tmp)[match(colnames(abundances), tmp$Sample),]
y.ab <- DGEList(abundances, samples=extra.info)

keep <- filterByExpr(y.ab, group=y.ab$samples$Treatment)
summary(keep)
y.ab <- y.ab[keep,]

design <- model.matrix(~Treatment, y.ab$samples)

y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion)
plotBCV(y.ab, cex=1)

fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)

# NOTE: With only 2 labels this doesn't work / make sense
summary(fit.ab$df.prior)
plotQLDisp(fit.ab, cex=1)

res <- glmQLFTest(fit.ab, coef = "TreatmentUninfected")
summary(decideTests(res))

topTags(res)

p1 <- plotUMAP(
  tmp,
  colour_by = "subcluster",
  text_by = "subcluster",
  point_size = 0.5,
  point_alpha = 1)
p2 <- plotUMAP(
  tmp,
  colour_by = "Treatment",
  point_size = 0.5,
  point_alpha = 1) +
  scale_colour_manual(values = treatment_colours)
p1 + p2
p3 <- ggplot(as.data.frame(colData(tmp)[, c("subcluster", "Sample")])) +
  geom_bar(
    aes(x = subcluster, fill = Sample),
    position = position_fill(reverse = TRUE)) +
  coord_flip() +
  ylab("Frequency") +
  scale_fill_manual(values = sample_colours) +
  theme_cowplot()
p4 <- ggplot(as.data.frame(colData(tmp)[, c("subcluster", "Treatment")])) +
  geom_bar(
    aes(x = subcluster, fill = Treatment),
    position = position_fill(reverse = TRUE)) +
  coord_flip() +
  ylab("Frequency") +
  scale_fill_manual(values = treatment_colours) +
  theme_cowplot()
p5 <- ggplot(as.data.frame(colData(tmp)[, "subcluster", drop = FALSE])) +
  geom_bar(aes(x = subcluster, fill = subcluster)) +
  coord_flip() +
  ylab("Number of cells") +
  theme_cowplot() +
  guides(fill = FALSE)
p3 + p4 + p5 + plot_layout(ncol = 3)

# DA analysis using subcluster labels in cycling clusters ----------------------

tmp <- cycling_sce
abundances <- table(tmp$subcluster, tmp$Sample)
abundances <- unclass(abundances)
abundances

extra.info <- colData(tmp)[match(colnames(abundances), tmp$Sample),]
y.ab <- DGEList(abundances, samples=extra.info)

keep <- filterByExpr(y.ab, group=y.ab$samples$Treatment)
summary(keep)
y.ab <- y.ab[keep,]

design <- model.matrix(~Treatment, y.ab$samples)

y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion)
plotBCV(y.ab, cex=1)

fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)

# NOTE: With only 2 labels this doesn't work / make sense
summary(fit.ab$df.prior)
plotQLDisp(fit.ab, cex=1)

res <- glmQLFTest(fit.ab, coef = "TreatmentUninfected")
summary(decideTests(res))

topTags(res)

p1 <- plotUMAP(
  tmp,
  colour_by = "subcluster",
  text_by = "subcluster",
  point_size = 0.5,
  point_alpha = 1)
p2 <- plotUMAP(
  tmp,
  colour_by = "Treatment",
  point_size = 0.5,
  point_alpha = 1) +
  scale_colour_manual(values = treatment_colours)
p1 + p2
p3 <- ggplot(as.data.frame(colData(tmp)[, c("subcluster", "Sample")])) +
  geom_bar(
    aes(x = subcluster, fill = Sample),
    position = position_fill(reverse = TRUE)) +
  coord_flip() +
  ylab("Frequency") +
  scale_fill_manual(values = sample_colours) +
  theme_cowplot()
p4 <- ggplot(as.data.frame(colData(tmp)[, c("subcluster", "Treatment")])) +
  geom_bar(
    aes(x = subcluster, fill = Treatment),
    position = position_fill(reverse = TRUE)) +
  coord_flip() +
  ylab("Frequency") +
  scale_fill_manual(values = treatment_colours) +
  theme_cowplot()
p5 <- ggplot(as.data.frame(colData(tmp)[, "subcluster", drop = FALSE])) +
  geom_bar(aes(x = subcluster, fill = subcluster)) +
  coord_flip() +
  ylab("Number of cells") +
  theme_cowplot() +
  guides(fill = FALSE)
p3 + p4 + p5 + plot_layout(ncol = 3)

# TODOs ------------------------------------------------------------------------

# - [x] Cell cycle could really be proliferation of T-cells upon activation.
# - [x] Are samples paired?
#   - No
# - [ ] Heatmaps and CSVs of DEGs (UP TO HERE)
# - [x] Block on treatment/donor during marker gene detection?
# - [ ] Convert to Rmd (and split into annotation and DE/DA)
