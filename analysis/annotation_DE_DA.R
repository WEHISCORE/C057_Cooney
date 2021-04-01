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

# Re-process original data -----------------------------------------------------

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
    ylab("Number of droplets") +
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
# TODO: Check denolisePCA works
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
    ylab("Number of droplets") +
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
    ylab("Number of droplets") +
    theme_cowplot(font_size = 8) +
    guides(fill = FALSE),
  align = "h",
  ncol = 3)

# MNN by treatment -------------------------------------------------------------

set.seed(1000101001)
mnn_treatment.out <- fastMNN(sce, batch = sce$Treatment, d=50, k=20, subset.row=hvg)
snn.gr <- buildSNNGraph(mnn_treatment.out, use.dimred="corrected")
clusters.mnn_treatment <- igraph::cluster_louvain(snn.gr)$membership
table(Cluster=clusters.mnn_treatment, treatment=mnn_treatment.out$batch)
mnn_treatment.out <- runUMAP(mnn_treatment.out, dimred="corrected")
mnn_treatment.out$batch <- factor(mnn_treatment.out$batch)
mnn_treatment.out$cluster <- factor(clusters.mnn_treatment)
mnn_treatment.p1 <- plotUMAP(mnn_treatment.out, colour_by = "batch", point_size = 0.5) +
  scale_colour_manual(values = treatment_colours, name = "Treatment") +
  ggtitle("MNN (treatment)")
mnn_treatment.p2 <- plotUMAP(mnn_treatment.out, colour_by = "cluster", point_size = 0.5) +
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
    ylab("Number of droplets") +
    theme_cowplot(font_size = 8) +
    guides(fill = FALSE),
  align = "h",
  ncol = 3)

# MNN by sample -----------------------------------------------------------------

set.seed(1000101001)
mnn_sample.out <- fastMNN(sce, batch = sce$Sample, d=50, k=20, subset.row=hvg.sample)
snn.gr <- buildSNNGraph(mnn_sample.out, use.dimred="corrected")
clusters.mnn_sample <- igraph::cluster_louvain(snn.gr)$membership
table(Cluster=clusters.mnn_sample, treatment=mnn_sample.out$batch)
mnn_sample.out <- runUMAP(mnn_sample.out, dimred="corrected")
mnn_sample.out$batch <- factor(mnn_sample.out$batch)
mnn_sample.out$cluster <- factor(clusters.mnn_sample)
mnn_sample.p1 <- plotUMAP(mnn_sample.out, colour_by = I(sce$Treatment), point_size = 0.5) +
  scale_colour_manual(values = treatment_colours, name = "Treatment") +
  ggtitle("MNN (sample)")
mnn_sample.p2 <- plotUMAP(mnn_sample.out, colour_by = "cluster", point_size = 0.5) +
  ggtitle("MNN (sample)")
mnn_sample.p3 <- plotUMAP(mnn_sample.out, colour_by = "batch", point_size = 0.5) +
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
    ylab("Number of droplets") +
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

p1 <- plotUMAP(sce, colour_by = "cluster") +
  scale_colour_manual(values = cluster_colours, name = "cluster") +
  ggtitle("Uncorrected: k = 10")
p2 <- plotUMAP(sce, colour_by = "cluster_k20") +
  ggtitle("Uncorrected: k = 20")
p3 <- plotUMAP(sce, colour_by = "cluster_k50") +
  ggtitle("Uncorrected: k = 50")
p4 <- plotUMAP(sce, colour_by = "cluster_k100") +
  ggtitle("Uncorrected: k = 100")
p1 + p2 + p3 + p4 + plot_layout(ncol = 2)

# NOTE: In all clusterings we observe an 'out group' on the right hand side.

# Look at cycling cluster(s) ---------------------------------------------------

# NOTE: Cell cycle could really be proliferation of T-cells upon activation

# Appears to be a cycling and non-cycling clustering driving the results.
# Cycling cluster(s) typified by upregulation of a large number of genes, larger
# library sizes, and cyclin expression.
# Non-cycling clusters typified by opposite and no upregulation of specific
# markers.

cycling_cluster <- "8"
sce$cycling_cluster <- ifelse(
  sce$cluster %in% cycling_cluster,
  "Cycling",
  "Not cycling")
plotUMAP(sce, colour_by = "cycling_cluster")
z <- findMarkers(sce, sce$cycling_cluster, direction = "up")

cyclin_genes <- grep("^CCN[ABDE][0-9]$", rowData(sce)$Symbol)
cyclin_genes <- sort(rownames(sce)[cyclin_genes])
cyclin_genes

sce$ls <- log10(sce$sum)
plotHeatmap(
  sce,
  c(
    head(rownames(z[["Not cycling"]]), 25),
    head(rownames(z[["Cycling"]]), 25),
    cyclin_genes),
  order_columns_by = c("cycling_cluster", "cluster"),
  cluster_rows = FALSE,
  color = hcl.colors(101, "Blue-Red 3"),
  center = TRUE,
  zlim = c(-3, 3),
  colour_columns_by = "ls",
  colour_annotation_colours = list(
    cluster = cluster_colours))

plot_grid(
  ggplot(as.data.frame(colData(sce)[, c("cycling_cluster", "Sample")])) +
    geom_bar(
      aes(x = cycling_cluster, fill = Sample),
      position = position_fill(reverse = TRUE)) +
    coord_flip() +
    ylab("Frequency") +
    scale_fill_manual(values = sample_colours) +
    theme_cowplot(font_size = 8),
  ggplot(as.data.frame(colData(sce)[, c("cycling_cluster", "Treatment")])) +
    geom_bar(
      aes(x = cycling_cluster, fill = Treatment),
      position = position_fill(reverse = TRUE)) +
    coord_flip() +
    ylab("Frequency") +
    scale_fill_manual(values = treatment_colours) +
    theme_cowplot(font_size = 8),
  ggplot(as.data.frame(colData(sce)[, "cycling_cluster", drop = FALSE])) +
    geom_bar(aes(x = cycling_cluster, fill = cycling_cluster)) +
    coord_flip() +
    ylab("Number of droplets") +
    theme_cowplot(font_size = 8) +
    guides(fill = FALSE),
  align = "h",
  ncol = 3)

# Subset analysis based on cycling cluster(s) ----------------------------------

set.seed(9391)
list_of_sce <- quickSubCluster(
  sce,
  groups = sce$cycling_cluster,
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
    # TODO: Experiment with k
    snn_gr <- buildSNNGraph(x, use.dimred = "PCA", k = 50)
    factor(igraph::cluster_louvain(snn_gr)$membership)
  })
# It's also useful to have per-sample UMAP representations.
set.seed(17127)
list_of_sce <- lapply(list_of_sce, runUMAP, dimred = "PCA")

plot_grid(
  plotlist = lapply(names(list_of_sce), function(n) {
    x <- list_of_sce[[n]]
    plotUMAP(x, colour_by = "subcluster") +
      ggtitle(n)
  }))

x <- list_of_sce[["Not cycling"]]
m <- findMarkers(x, x$subcluster, direction = "up", pval.type = "all")
sapply(m, function(x) sum(x$FDR < 0.05))
features <- unique(unlist(lapply(m, function(x) head(rownames(x), 10))))
plotHeatmap(
  x,
  features = features,
  order_columns_by = c("subcluster", "Treatment", "Sample"),
  color = hcl.colors(101, "Blue-Red 3"),
  center = TRUE,
  zlim = c(-3, 3),
  column_annotation_colors = list(
    Treatment = treatment_colours,
    Sample = sample_colours),
  cluster_rows = FALSE,
  annotation_row = data.frame(
    subcluster = rep(names(m), each = length(features) / length(m)),
    row.names = features))

x <- list_of_sce[["Cycling"]]
m <- findMarkers(x, x$subcluster, direction = "up", pval.type = "all")
features <- unique(unlist(lapply(m, function(x) head(rownames(x), 10))))
plotHeatmap(
  x,
  features = features,
  order_columns_by = c("subcluster", "Treatment", "Sample"),
  color = hcl.colors(101, "Blue-Red 3"),
  center = TRUE,
  zlim = c(-3, 3),
  column_annotation_colors = list(
    Treatment = treatment_colours,
    Sample = sample_colours),
  cluster_rows = FALSE,
  annotation_row = data.frame(
    subcluster = rep(names(m), each = length(features) / length(m)),
    row.names = features))

# Using the cyclins ------------------------------------------------------------

# - Cyclin A is expressed across S and G2
# - Cyclin B is expressed highest in late G2 and mitosis (Morgan 2007).
# - Cyclin D is expressed throughout but peaks at G1
# - Cyclin E is expressed highest in the G1/S transition

plotHeatmap(
  sce,
  order_columns_by = c("cycling_cluster", "cluster"),
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
markers[[cycling_cluster]]

# Using the cyclone() classifier -----------------------------------------------

# TODO: Save output from run on vc7-shared and load here.
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

plot(
  x = assignments$score$G1,
  y = assignments$score$G2M,
  xlab = "G1 score",
  ylab = "G2/M score",
  pch = 16,
  col = sce$cluster_colours)

boxplot(assignments$scores$G1 ~ sce$cluster, ylab = "G1")
boxplot(assignments$scores$G2M ~ sce$cluster, ylab = "G2M")

proportions(table(sce$cluster, assignments$phases), 1)
proportions(table(sce$cycling_cluster, assignments$phases), 1)

plotUMAP(
  sce,
  colour_by = data.frame(G1 = assignments$score$G1),
  text_by = "cluster") +
  plotUMAP(
    sce,
    colour_by = data.frame(G2M = assignments$score$G2M),
    text_by = "cluster") +
  plotUMAP(sce, colour_by = "cluster", text_by = "cluster") +
  scale_colour_manual(values = cluster_colours, name = "cluster") +
  plot_layout(ncol = 2)

# DE analysis using cycling vs. non-cycling labels -----------------------------

library(edgeR)
plotMDS <- scater::plotMDS

# NOTE: Have to drop the complicated TRA and TRB DFrameList objects.
tmp <- sce
tmp$TRA <- NULL
tmp$TRB <- NULL
summed <- aggregateAcrossCells(
  tmp,
  id = colData(tmp)[, c("cycling_cluster", "Sample")])
sizeFactors(summed) <- NULL
summed <- logNormCounts(summed)

summed <- runMDS(summed)
p1 <- plotMDS(summed, colour_by = "Treatment", shape_by = "cycling_cluster") +
  scale_colour_manual(values = treatment_colours, name = "Treatment")
p2 <- plotMDS(summed, colour_by = "Sample", shape_by = "cycling_cluster") +
  scale_colour_manual(values = sample_colours, name = "Sample")
p1 + p2 + plot_layout(guides = "collect")

summed.filt <- summed[,summed$ncells >= 10]

de.results <- pseudoBulkDGE(
  summed.filt,
  label = summed.filt$cycling_cluster,
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
  other_fields = "cycling_cluster") +
  facet_wrap(~cycling_cluster) +
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
  label = summed.filt$cycling_cluster,
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
p1 <- plotMDS(summed, colour_by = "cluster", shape_by = "Treatment") +
  scale_colour_manual(values = cluster_colours, name = "cluster")
p2 <- plotMDS(summed, colour_by = "cluster", shape_by = "Sample") +
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
tmp <- list_of_sce$`Not cycling`
tmp$TRA <- NULL
tmp$TRB <- NULL
summed <- aggregateAcrossCells(
  tmp,
  id = colData(tmp)[, c("subcluster", "Sample")])
sizeFactors(summed) <- NULL
summed <- logNormCounts(summed)

summed <- runMDS(summed)
p1 <- plotMDS(summed, colour_by = "subcluster", shape_by = "Treatment")
p2 <- plotMDS(summed, colour_by = "subcluster", shape_by = "Sample")
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
tmp <- list_of_sce$Cycling
tmp$TRA <- NULL
tmp$TRB <- NULL
summed <- aggregateAcrossCells(
  tmp,
  id = colData(tmp)[, c("subcluster", "Sample")])
sizeFactors(summed) <- NULL
summed <- logNormCounts(summed)

summed <- runMDS(summed)
p1 <- plotMDS(summed, colour_by = "subcluster", shape_by = "Treatment")
p2 <- plotMDS(summed, colour_by = "subcluster", shape_by = "Sample")
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

# TODO: Need to think about how to do this with only 2 labels.

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
plotUMAP(sce, colour_by = "cluster", text_by = "cluster") +
  scale_colour_manual(values = cluster_colours) +
  plotUMAP(sce, colour_by = "Treatment") +
  scale_colour_manual(values = treatment_colours)

# DA analysis using subcluster labels in non-cycling clusters ------------------

tmp <- list_of_sce$`Not cycling`
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

plotUMAP(tmp, colour_by = "subcluster", text_by = "subcluster") +
  plotUMAP(tmp, colour_by = "Treatment") +
  scale_colour_manual(values = treatment_colours)

# DA analysis using subcluster labels in cycling clusters ----------------------

tmp <- list_of_sce$Cycling
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

plotUMAP(tmp, colour_by = "subcluster", text_by = "subcluster") +
  plotUMAP(tmp, colour_by = "Treatment") +
  scale_colour_manual(values = treatment_colours)

# TODOs ------------------------------------------------------------------------

# - [ ] Cell cycle could really be proliferation of T-cells upon activation.
# - [ ] Are samples paired?
# - [ ] Ambient RNA
# - [ ] Heatmaps and CSVs of DEGs
