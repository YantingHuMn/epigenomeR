# Automatically performs k-means clustering, generates cluster files internally, and draw heatmap
# Post: create_cluster_heatmap can be replaced by this function
# Post: Advanced wrapper function that combines automated k-means clustering with heatmap generation. Performs clustering analysis, saves cluster assignments, and generates publication-ready heatmaps in a single workflow.
# Post: include function rowcol_km_like_ComplexHeatmap; apply_cluster_heatmap; combine_count_matrix
# Parameters: count_matrix_file_path: Path to count matrix .feather file (rows=genomic positions, cols=samples)
#             row_km: Number of k-means clusters for rows (genomic features)
#             col_km: Number of k-means clusters for columns (samples)
#             out_dir: Directory to save cluster files and heatmap output
#             seed: Random seed for reproducible clustering (default: 123)
#             plot: whether to generate heatmap plot
#             show_column_names: whether to show column names at the bottom of the heatmap (default: FALSE)
#             lower_range: Lower bound for heatmap color scale (default: NULL, auto-determined)
#             upper_range: Upper bound for heatmap color scale (default: NULL, auto-determined)
#             row_title_fontsize: Font size for row cluster titles (default: NULL)
#             col_title_fontsize: Font size for column cluster titles (default: NULL)
#             legend_title_fontsize: Font size for legend title (default: NULL)
#             legend_label_fontsize: Font size for legend labels (default: NULL)
# Output: None (saves cluster assignment tables and heatmap plots to output directory)
biclustering <- function(cm_path, row_km, col_km, out_dir, seed = 42, plot = TRUE, show_column_names = FALSE, lower_range = NULL, upper_range = NULL, row_title_fontsize = NULL, col_title_fontsize = NULL, legend_title_fontsize = NULL, legend_label_fontsize = NULL) {
  library(arrow)
  library(tibble)
  library(dplyr)

  if (is.null(cm_path) || length(cm_path) == 0) {
    stop("`count_matrix_file_path` is required", call. = FALSE)
  }
  mat <- as.matrix(column_to_rownames(read_feather(cm_path), var = "pos"))

  message("Performing bidirectional k-means clustering...")
  result <- bidirectional_kmeans_clustering(mat = mat, row_k = row_km, col_k = col_km, seed = seed)
  row_letter <- result$row_letter
  col_num <- result$col_num

  df_row <- data.frame(
    feature = names(row_letter),
    label   = unname(row_letter),
    stringsAsFactors = FALSE
  )

  df_col <- data.frame(
    feature = names(col_num),
    label   = unname(col_num),
    stringsAsFactors = FALSE
  )

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  path1 <- file.path(out_dir, "row_table.tsv")
  path2 <- file.path(out_dir, "col_table.tsv")
  write.table(df_row, path1, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(df_col, path2, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Saved cluster assignments to:")
  message("  - row clusters: ", path1)
  message("  - col clusters: ", path1)
  
  if (plot) {
    message("Generating heatmap...")
    mat_sorted <- mat[names(result$row_letter), names(result$col_num)]
    biclustering_heatmap(mat = mat_sorted, row_cluster_file_path = path1, col_cluster_file_path = path2, out_dir = out_dir, show_column_names = show_column_names, lower_range = lower_range, upper_range = upper_range, row_title_fontsize = row_title_fontsize, col_title_fontsize = col_title_fontsize, legend_title_fontsize = legend_title_fontsize, legend_label_fontsize = legend_label_fontsize)
  }

  return(list("row_table" = path1, "col_table" = path2))
}
