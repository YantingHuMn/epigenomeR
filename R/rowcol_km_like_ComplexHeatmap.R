# Generate row and col cluster groups like ComplexHeatmap
# Post: Performs k-means clustering on matrix rows and columns, then hierarchically clusters the resulting groups to create ordered cluster assignments similar to ComplexHeatmap package output.
# Parameter: mat: Input numeric matrix for clustering
#            row_k: Number of row clusters
#            col_k: Number of column clusters
#            seed: Random seed for reproducible clustering (default: 42)
#            row_repeats: Number of k-means repetitions for row consensus clustering (default: 1)
#            col_repeats: Number of k-means repetitions for column consensus clustering (default: 1)
#            cluster_row_slices: Whether to hierarchically cluster row groups (default: TRUE)
#            cluster_col_slices: Whether to hierarchically cluster column groups (default: TRUE)
#            distance_method_for_slices: Distance method for hierarchical clustering ("euclidean", "manhattan", etc.)
#            hclust_method_for_slices: Linkage method for hierarchical clustering ("complete", "average", "single", etc.)
#            do_reorder_rows: Whether to reorder rows within clusters based on weights (default: TRUE)
#            reorder_rows_weight: Custom weights for row reordering (default: NULL, uses negative row means)
# Output: List containing row_letter (LETTERS A,B,C... for row clusters) and col_num (numbers 1,2,3... for column clusters)
rowcol_km_like_ComplexHeatmap <- function(mat, row_k, col_k, seed = 42, row_repeats = 1, col_repeats = 1, cluster_row_slices = TRUE, cluster_col_slices = TRUE, distance_method_for_slices = "euclidean", hclust_method_for_slices = "complete", do_reorder_rows = TRUE, reorder_rows_weight = NULL) {

  stopifnot(is.matrix(mat), is.numeric(mat))
  if (is.null(rownames(mat))) rownames(mat) <- sprintf("row_%d", seq_len(nrow(mat)))
  if (is.null(colnames(mat))) colnames(mat) <- sprintf("col_%d", seq_len(ncol(mat)))

  set.seed(seed)

  consensus_kmeans <- function(m, centers, reps) {
    parts <- lapply(seq_len(reps), function(i) {
      clue::as.cl_hard_partition(stats::kmeans(m, centers, iter.max = 50))
    })
    cons <- clue::cl_consensus(clue::cl_ensemble(list = parts))
    as.vector(clue::cl_class_ids(cons))
  }

  # row cluster: output letter
  row_cl <- consensus_kmeans(mat, row_k, row_repeats)

  row_meanmat <- sapply(sort(unique(row_cl)), function(i) {
    colMeans(mat[row_cl == i, , drop = FALSE], na.rm = TRUE)
  })
  if (!is.matrix(row_meanmat)) row_meanmat <- matrix(row_meanmat, nrow = 1)

  if (cluster_row_slices) {
    hc_row <- stats::hclust(stats::dist(t(row_meanmat), method = distance_method_for_slices), method = hclust_method_for_slices)
    w_row  <- colMeans(row_meanmat)
    ord_row <- stats::order.dendrogram(stats::reorder(stats::as.dendrogram(hc_row), w_row, mean))
  } else {
    ord_row <- order(colMeans(row_meanmat))
  }

  row_cl2 <- match(row_cl, ord_row)
  row_cl2 <- factor(row_cl2, levels = seq_along(ord_row))
  row_num <- as.integer(as.character(row_cl2))
  names(row_num) <- rownames(mat)
  row_letter <- LETTERS[row_num]
  names(row_letter) <- rownames(mat)

  # reorder row cluster
  if (isTRUE(do_reorder_rows) && is.null(reorder_rows_weight)) {
    reorder_rows_weight <- -rowMeans(mat, na.rm = TRUE)
  }

  row_order <- integer(0)
  for (sid in seq_along(ord_row)) {
    idx <- which(row_num == sid)
    if (length(idx) <= 1) {
      row_order <- c(row_order, idx)
    } else {
      submat <- mat[idx, , drop = FALSE]
      hc <- stats::hclust(stats::dist(submat, method = distance_method_for_slices), method = hclust_method_for_slices)
      dg <- stats::as.dendrogram(hc)
      if (isTRUE(do_reorder_rows)) dg <- stats::reorder(dg, reorder_rows_weight[idx], mean)
      od <- stats::order.dendrogram(dg)
      row_order <- c(row_order, idx[od])
    }
  }

  row_letter <- row_letter[row_order]

  # col cluster: output num
  col_cl <- consensus_kmeans(t(mat), col_k, col_repeats)

  col_meanmat <- sapply(sort(unique(col_cl)), function(i) {
    rowMeans(mat[, col_cl == i, drop = FALSE], na.rm = TRUE)
  })
  if (!is.matrix(col_meanmat)) col_meanmat <- matrix(col_meanmat, nrow = 1)

  if (cluster_col_slices) {
    hc_col <- stats::hclust(stats::dist(t(col_meanmat), method = distance_method_for_slices), method = hclust_method_for_slices)
    w_col  <- colMeans(col_meanmat)
    ord_col <- stats::order.dendrogram(stats::reorder(stats::as.dendrogram(hc_col), w_col, mean))
  } else {
    ord_col <- order(colMeans(col_meanmat))
  }

  col_cl2 <- match(col_cl, ord_col)
  col_cl2 <- factor(col_cl2, levels = seq_along(ord_col))
  col_num <- as.integer(as.character(col_cl2))
  names(col_num) <- colnames(mat)

  # reorder columns within each column-slice, analogous to rows
  reorder_cols_weight <- -colMeans(mat, na.rm = TRUE)
  col_order <- integer(0)
  for (sid in seq_along(ord_col)) {
    idx <- which(col_num == sid)
    if (length(idx) <= 1) {
      col_order <- c(col_order, idx)
    } else {
      submat <- mat[, idx, drop = FALSE]
      hc <- stats::hclust(stats::dist(t(submat),
                                      method = distance_method_for_slices),
                          method = hclust_method_for_slices)
      dg <- stats::as.dendrogram(hc)
      dg <- stats::reorder(dg, reorder_cols_weight[idx], mean)
      od <- stats::order.dendrogram(dg)
      col_order <- c(col_order, idx[od])
    }
  }
  col_num <- col_num[col_order]

  list(row_letter = row_letter, col_num = col_num)
}
