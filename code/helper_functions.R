# Helper function to coerce a DataFrame to a data.frame while preserving column
# names as is.
.adf <- function(x) {
  setNames(as.data.frame(x), colnames(x))
}

# Helper function to Combine data from 2 SCEs using gene names.
# NOTE: This assumes more than I'd like about the rowData and doesn't do much
#       checking of these assumptions.
.combine <- function(x, y, rowData_by = c("ENSEMBL", "SYMBOL", "CHR")) {
  if (is.null(rowData_by)) {
    rowData <- dplyr::full_join(
      .adf(rowData(x)) %>%
        tibble::rownames_to_column(var = "gene"),
      .adf(rowData(y)) %>%
        tibble::rownames_to_column(var = "gene")) %>%
      tibble::column_to_rownames("gene") %>%
      DataFrame(., row.names = rownames(.))
  } else {
    rowData <- dplyr::full_join(
      .adf(rowData(x)[, rowData_by, drop = FALSE]),
      .adf(rowData(y)[, rowData_by, drop = FALSE]),
      by = rowData_by) %>%
      DataFrame(row.names = scater::uniquifyFeatureNames(
        .$ENSEMBL,
        .$SYMBOL))
    rownames(x) <- rownames(rowData)[match(rowData(x)$ENSEMBL, rowData$ENSEMBL)]
    rownames(y) <- rownames(rowData)[match(rowData(y)$ENSEMBL, rowData$ENSEMBL)]
  }

  colData <- rbind(colData(x), colData(y))

  counts <- matrix(
    data = 0L,
    nrow = nrow(rowData), ncol = nrow(colData),
    dimnames = list(rownames(rowData), rownames(colData)))
  counts[rownames(x), colnames(x)] <- counts(
    x,
    withDimnames = FALSE)
  counts[rownames(y), colnames(y)] <- counts(
    y,
    withDimnames = FALSE)

  stopifnot(
    identical(
      metadata(x)$scPipe$version,
      metadata(y)$scPipe$version))
  stopifnot(
    identical(
      metadata(x)$scPipe$QC_cols,
      metadata(y)$scPipe$QC_cols))
  stopifnot(
    identical(
      metadata(x)$scPipe$demultiplex_info$status,
      metadata(y)$scPipe$demultiplex_info$status))
  stopifnot(
    identical(
      metadata(x)$scPipe$UMI_dup_info$duplication.number,
      metadata(y)$scPipe$UMI_dup_info$duplication.number))
  stopifnot(identical(metadata(x)$Biomart, metadata(y)$Biomart))
  metadata <- list(
    scPipe = list(
      version = metadata(x)$scPipe$version,
      QC_cols = metadata(x)$scPipe$QC_cols,
      demultiplex_info = data.frame(
        status = metadata(x)$scPipe$demultiplex_info$status,
        count = metadata(x)$scPipe$demultiplex_info$count +
          metadata(y)$scPipe$demultiplex_info$count),
      UMI_dup_info = data.frame(
        duplication.number = metadata(
          x)$scPipe$UMI_dup_info$duplication.number,
        count = metadata(x)$scPipe$UMI_dup_info$count +
          metadata(y)$scPipe$UMI_dup_info$count)),
    Biomart = metadata(x)$Biomart)

  sce <- SingleCellExperiment(
    rowData = rowData,
    colData = colData,
    assays = list(counts = counts),
    metadata = metadata)

  stopifnot(identical(int_metadata(x), int_metadata(y)))
  int_metadata(sce) <- int_metadata(x)
  int_elementMetadata <- dplyr::full_join(
    x = .adf(int_elementMetadata(x)) %>%
      tibble::add_column(gene = rownames(x)),
    y = .adf(int_elementMetadata(y)) %>%
      tibble::add_column(gene = rownames(y))) %>%
    tibble::column_to_rownames("gene") %>%
    DataFrame()
  int_elementMetadata(sce) <- int_elementMetadata

  stopifnot(validObject(sce))
  sce
}

.cbindSCEs <- function(list_of_sce, rowData_by = 1:6) {
  do.call(
    cbind,
    lapply(list_of_sce, function(sce) {
      # NOTE: Some fudging to combine only the necessary bits of each SCE
      #       (basically, don't include any QC metrics).
      rowData(sce) <- rowData(sce)[, rowData_by]
      sce
    }))
}

# NOTE: My best guess of what this function does https://github.com/MarioniLab/compareSingleCell/blob/543aa28e3ae25fad4ffb1d47c27c8a364966095c/vignettes/embryo_expression.Rmd#L76
.sumCountsAcrossCells <- function(sce, cluster_sample) {
  counts <- counts(sce, withDimnames = FALSE)
  edgeR::sumTechReps(counts, cluster_sample)
}

# NOTE: A modified version of scater::plotHeatmap() that allows me to pass
#       `annotation_colors` down to pheatmap::pheatmap().
.plotHeatmap <- function (object, features, columns = NULL, exprs_values = "logcounts",
                          center = FALSE, zlim = NULL, symmetric = FALSE, color = NULL,
                          colour_columns_by = NULL, by_exprs_values = exprs_values,
                          by_show_single = FALSE, show_colnames = TRUE, ...) {
  features_to_use <- scater:::.subset2index(features, object, byrow = TRUE)
  heat.vals <- assay(
    object,
    exprs_values,
    withDimnames = FALSE)[features_to_use, , drop = FALSE]
  rownames(heat.vals) <- rownames(object)[features_to_use]
  if (is.null(colnames(object))) {
    colnames(heat.vals) <- seq_len(ncol(object))
    show_colnames <- FALSE
  } else {
    colnames(heat.vals) <- colnames(object)
  }
  if (!is.null(columns)) {
    heat.vals <- heat.vals[, columns, drop = FALSE]
  }
  if (center) {
    heat.vals <- heat.vals - DelayedMatrixStats::rowMeans2(DelayedArray(heat.vals))
  }
  if (is.null(zlim)) {
    zlim <- range(heat.vals)
  }
  if (symmetric) {
    extreme <- max(abs(zlim))
    zlim <- c(-extreme, extreme)
  }
  heat.vals[heat.vals < zlim[1]] <- zlim[1]
  heat.vals[heat.vals > zlim[2]] <- zlim[2]
  if (is.null(color)) {
    color <- eval(formals(pheatmap::pheatmap)$color, envir = environment(pheatmap::pheatmap))
  }
  color.breaks <- seq(zlim[1], zlim[2], length.out = length(color) +
                        1L)
  if (length(colour_columns_by)) {
    column_variables <- column_colorings <- list()
    for (field in colour_columns_by) {
      colour_by_out <- scater:::.choose_vis_values(object, field,
                                                   mode = "column", search = "any", exprs_values = by_exprs_values,
                                                   discard_solo = !by_show_single)
      if (is.null(colour_by_out$val)) {
        next
      }
      else if (is.numeric(colour_by_out$val)) {
        colour_fac <- cut(colour_by_out$val, 25)
      } else {
        colour_fac <- as.factor(colour_by_out$val)
      }
      nlevs_colour_by <- nlevels(colour_fac)
      # if (nlevs_colour_by <= 10) {
      #   col_scale <- scater:::.get_palette("tableau10medium")
      # } else if (nlevs_colour_by > 10 && nlevs_colour_by <=
      #          20) {
      #   col_scale <- scater:::.get_palette("tableau20")
      # } else {
      #   col_scale <- viridis::viridis(nlevs_colour_by)
      # }
      # col_scale <- col_scale[seq_len(nlevs_colour_by)]
      # names(col_scale) <- levels(colour_fac)
      column_variables[[colour_by_out$name]] <- colour_fac
      # column_colorings[[colour_by_out$name]] <- col_scale
    }
    column_variables <- do.call(data.frame, c(column_variables,
                                              list(row.names = colnames(object))))
  } else {
    # column_variables <- column_colorings <- NULL
  }
  pheatmap::pheatmap(heat.vals, color = color, breaks = color.breaks,
                     annotation_col = column_variables,
                     # annotation_colors = column_colorings,
                     show_colnames = show_colnames, ...)
}

# Plot SingleR scores on a reduced dimension plot from a SCE.
plotScoreReducedDim <- function(results, sce, use_dimred = "TSNE",
                                max.labels = 20, normalize = TRUE, ncol = 5,
                                ...) {
  scores <- results$scores
  rownames(scores) <- rownames(results)
  m <- rowMaxs(scale(t(scores)))
  to.keep <- head(order(m, decreasing = TRUE), max.labels)
  if (normalize) {
    mmax <- rowMaxs(scores)
    mmin <- rowMins(scores)
    scores <- (scores - mmin) / (mmax - mmin)
    scores <- scores ^ 3
  }
  scores <- scores[, to.keep, drop = FALSE]
  cns <- colnames(scores)
  p <- lapply(cns, function(cn) {
    scater::plotReducedDim(
      sce,
      use_dimred = use_dimred,
      colour_by = data.frame(Score = scores[, cn]),
      ...) +
      ggtitle(cn) +
      scale_fill_viridis_c(limits = force(if(normalize) c(0, 1) else NULL)) +
      guides(fill = guide_colourbar(title = "Score"))
  })
    cowplot::plot_grid(plotlist = p, ncol = ncol)
}

# NOTE: A modified version of SingleR::plotScoreHeatmap() that allows me to
#       pass annotations down to pheatmap::pheatmap().
.plotScoreHeatmap <- function (results, cells.use = NULL, labels.use = NULL,
                               clusters = NULL, max.labels = 40,
                               normalize = TRUE, cells.order = NULL,
                               order.by.clusters = FALSE, ...) {
  if (is.null(rownames(results))) {
    rownames(results) <- seq_len(nrow(results))
  }
  scores <- results$scores
  rownames(scores) <- rownames(results)
  if (!is.null(cells.use)) {
    scores <- scores[cells.use, ]
  }
  if (!is.null(labels.use)) {
    scores <- scores[, labels.use]
  }
  m <- rowMaxs(scale(t(scores)))
  to.keep <- head(order(m, decreasing = TRUE), max.labels)
  if (normalize) {
    mmax <- rowMaxs(scores)
    mmin <- rowMins(scores)
    scores <- (scores - mmin)/(mmax - mmin)
    scores <- scores^3
  }
  scores <- scores[, seq_len(ncol(scores)) %in% to.keep, drop = FALSE]
  scores <- t(scores)
  if (!is.null(clusters)) {
    names(clusters) <- rownames(results)
    clusters <- data.frame(Clusters = clusters[colnames(scores)],
                           row.names = colnames(scores))
  }
  cluster_cols <- FALSE
  if (order.by.clusters && !is.null(clusters)) {
    order <- order(clusters$Clusters)
  }
  else if (!is.null(cells.order)) {
    order <- cells.order
  }
  else {
    order <- seq_len(ncol(scores))
    cluster_cols <- TRUE
  }
  args <- list(mat = scores[, order, drop = FALSE], border_color = NA,
               show_colnames = FALSE, clustering_method = "ward.D2",
               cluster_cols = cluster_cols, ...)
  if (!is.null(clusters)) {
    args$annotation_col <- clusters[order, , drop = FALSE]
  }
  do.call(pheatmap::pheatmap, args)
}

# Collapse labels only found in fewer than `cutoff` proportion of cells in all
# patients as 'other'
.collapseLabel <- function(labels, patient, cutoff = 0.01) {
  tmp <- table(labels, patient)
  tmp2 <- apply(tmp, 2, function(x) (x / sum(x)) > cutoff)
  tmp3 <- rownames(tmp2)[rowAnys(tmp2)]
  tmp4 <- ifelse(
    labels %in% tmp3,
    as.character(labels),
    "other")
  factor(tmp4, names(sort(table(tmp4), decreasing = TRUE)))
}
