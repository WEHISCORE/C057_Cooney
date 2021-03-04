# Takes off from cell selection

# MNN --------------------------------------------------------------------------

# Exploratory, probably over-corrects because we have such a homogeneous
# population?

library(batchelor)
sce$batch <- sce$Sample
set.seed(666)
mnn_out <- fastMNN(
  multiBatchNorm(sce, batch = sce$batch),
  batch = sce$batch,
  cos.norm = FALSE,
  d = ncol(reducedDim(sce, "PCA")),
  auto.merge = TRUE,
  subset.row = hvg,
  correct.all = TRUE)
reducedDim(sce, "corrected") <- reducedDim(mnn_out, "corrected")
assay(sce, "reconstructed") <- assay(mnn_out, "reconstructed")

var_loss <- metadata(mnn_out)$merge.info$lost.var
rownames(var_loss) <- sprintf("merge_%d", seq_len(nrow(var_loss)))
knitr::kable(
  100 * var_loss,
  digits = 1,
  caption = "Percentage of estimated biological variation lost within each donor at each step of the auto merge. Ideally, all these values should be small (e.g., < 10%).")

set.seed(853)
sce <- runUMAP(sce, dimred = "corrected")
set.seed(4759)
snn_gr <- buildSNNGraph(sce, use.dimred = "corrected", k = 5)
clusters <- igraph::cluster_louvain(snn_gr)
sce$cluster <- factor(clusters$membership)
cluster_colours <- setNames(
  scater:::.get_palette("tableau20"),
  levels(sce$cluster))
sce$cluster_colours <- cluster_colours[sce$cluster]

umap_df <- cbind(
  data.frame(
    x = reducedDim(sce, "UMAP")[, 1],
    y = reducedDim(sce, "UMAP")[, 2]),
  as.data.frame(colData(sce)[, !colnames(colData(sce)) %in% c("TRA", "TRB")]))

plot_grid(
  ggplot(aes(x = x, y = y), data = umap_df) +
    geom_point(aes(colour = cluster), size = 0.25) +
    scale_colour_manual(values = cluster_colours) +
    theme_cowplot(font_size = 8) +
    xlab("Dimension 1") +
    ylab("Dimension 2"),
  ggplot(aes(x = x, y = y), data = umap_df) +
    geom_point(aes(colour = Treatment), size = 0.25) +
    theme_cowplot(font_size = 8) +
    xlab("Dimension 1") +
    ylab("Dimension 2") +
    scale_colour_manual(values = treatment_colours),
  ggplot(aes(x = x, y = y), data = umap_df) +
    geom_point(aes(colour = Sample), size = 0.25) +
    theme_cowplot(font_size = 8) +
    xlab("Dimension 1") +
    ylab("Dimension 2") +
    scale_colour_manual(values = sample_colours),
  ncol = 2,
  align = "v")

plot_grid(
  ggplot(as.data.frame(colData(sce)[, c("cluster", "Sample")])) +
    geom_bar(
      aes(x = cluster, fill = Sample),
      position = position_fill(reverse = TRUE)) +
    coord_flip() +
    ylab("Frequency") +
    scale_fill_manual(values = sample_colours) +
    theme_cowplot(font_size = 5),
  ggplot(as.data.frame(colData(sce)[, c("cluster", "Treatment")])) +
    geom_bar(
      aes(x = cluster, fill = Treatment),
      position = position_fill(reverse = TRUE)) +
    coord_flip() +
    ylab("Frequency") +
    scale_fill_manual(values = treatment_colours) +
    theme_cowplot(font_size = 5),
  ggplot(as.data.frame(colData(sce)[, "cluster", drop = FALSE])) +
    geom_bar(aes(x = cluster, fill = cluster)) +
    coord_flip() +
    ylab("Number of droplets") +
    scale_fill_manual(values = cluster_colours) +
    theme_cowplot(font_size = 5) +
    guides(fill = FALSE),
  align = "h",
  ncol = 3)

# Look at cycling cluster(s) ---------------------------------------------------

# Appears to be a cycling and non-cycling clustering driving the results.
# Cycling cluster(s) typified by upregulation of a large number of genes, larger
# library sizes, and cyclin expression.
# Non-cycling clusters typified by opposite and no upregulation of specific
# markers.

# NOTE: cycling clusters depend MNN.
sce$cycling_cluster <- sce$cluster %in% c("2", "8")
# sce$cycling_cluster <- sce$cluster %in% c("7")
plotUMAP(sce, colour_by = "cycling_cluster")
z <- findMarkers(sce, sce$cycling_cluster, direction = "up")

cyclin_genes <- grep("^CCN[ABDE][0-9]$", rowData(sce)$Symbol)
cyclin_genes <- sort(rownames(sce)[cyclin_genes])
cyclin_genes

sce$ls <- log10(sce$sum)
plotHeatmap(
  sce,
  c(head(rownames(z[["FALSE"]]), 25), head(rownames(z[["TRUE"]]), 25), cyclin_genes),
  order_columns_by = c("cycling_cluster", "cluster"),
  cluster_rows = FALSE,
  color = hcl.colors(101, "Blue-Red 3"),
  center = TRUE,
  zlim = c(-3, 3),
  colour_columns_by = "ls")

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
    hvg <- hvg[!is_ribo]
    # NOTE: Keep the original dimensionality reduction around for downstream
    #       plotting.
    reducedDimNames(x) <- paste0("original_", reducedDimNames(x))
    sce <- denoisePCA(
      x,
      var_fit,
      subset.row = hvg)
    # TODO: MNN or not-MNN
    mnn_out <- fastMNN(
      multiBatchNorm(sce, batch = sce$batch),
      batch = sce$batch,
      cos.norm = FALSE,
      d = ncol(reducedDim(sce, "PCA")),
      auto.merge = TRUE,
      subset.row = hvg,
      correct.all = TRUE)
    reducedDim(sce, "corrected") <- reducedDim(mnn_out, "corrected")
    assay(sce, "reconstructed") <- assay(mnn_out, "reconstructed")
    sce
  },
  clusterFUN = function(x) {
    # TODO: MNN or not-MNN
    snn_gr <- buildSNNGraph(x, use.dimred = "corrected")
    factor(igraph::cluster_louvain(snn_gr)$membership)
  })
# It's also useful to have per-sample UMAP representations.
set.seed(17127)
# TODO: MNN or not-MNN
list_of_sce <- lapply(list_of_sce, runUMAP, dimred = "corrected")

plot_grid(
  plotlist = lapply(names(list_of_sce), function(n) {
    x <- list_of_sce[[n]]
    plotUMAP(x, colour_by = "subcluster") +
      ggtitle(ifelse(n, "Cycling", "Not cycling"))
  }))

x <- list_of_sce[["FALSE"]]
m <- findMarkers(x, x$subcluster, direction = "up", pval.type = "some")
plotHeatmap(
  x,
  unique(unlist(lapply(m, function(x) head(rownames(x), 10)))),
  order_columns_by = c("subcluster", "Treatment"),
  color = hcl.colors(101, "Blue-Red 3"),
  center = TRUE,
  zlim = c(-3, 3))

x <- list_of_sce[["TRUE"]]
m <- findMarkers(x, x$subcluster, direction = "up", pval.type = "some")
plotHeatmap(
  x,
  unique(unlist(lapply(m, function(x) head(rownames(x), 10)))),
  order_columns_by = c("subcluster", "Treatment"),
  color = hcl.colors(101, "Blue-Red 3"),
  center = TRUE,
  zlim = c(-3, 3))

# TODO: How to handle TRA/TRB-upregulated cluster in non-cycling?
