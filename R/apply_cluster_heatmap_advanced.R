# Automatically performs k-means clustering, generates cluster files internally, and draw heatmap
# Post: create_cluster_heatmap can be replaced by this function
# Post: Advanced wrapper function that combines automated k-means clustering with heatmap generation. Performs clustering analysis, saves cluster assignments, and generates publication-ready heatmaps in a single workflow.
# Post: include function rowcol_km_like_ComplexHeatmap; apply_cluster_heatmap; combine_count_matrix
# Parameters: count_matrix_file_path: Path to count matrix .feather file (rows=genomic positions, cols=samples)
#             row_km: Number of k-means clusters for rows (genomic features)
#             col_km: Number of k-means clusters for columns (samples)
#             output_dir_path: Directory to save cluster files and heatmap output
#             seed: Random seed for reproducible clustering (default: 123)
#             show_dend_boolean: Whether to show dendrograms in heatmap (default: FALSE)
#             count_matrix_overlap: Whether to merge shared positions across matrices (default: FALSE)
#             lower_range: Lower bound for heatmap color scale (default: NULL, auto-determined)
#             upper_range: Upper bound for heatmap color scale (default: NULL, auto-determined)
#             row_title_fontsize: Font size for row cluster titles (default: NULL)
#             col_title_fontsize: Font size for column cluster titles (default: NULL)
#             legend_title_fontsize: Font size for legend title (default: NULL)
#             legend_label_fontsize: Font size for legend labels (default: NULL)
# Output: None (saves cluster assignment tables and heatmap plots to output directory)
apply_cluster_heatmap_advanced <- function(count_matrix_file_path, row_km, col_km, output_dir_path, seed = 123, show_dend_boolean = FALSE, count_matrix_overlap= FALSE, lower_range = NULL, upper_range = NULL, row_title_fontsize = NULL, col_title_fontsize = NULL, legend_title_fontsize = NULL, legend_label_fontsize = NULL) {
  library(arrow)
  library(tibble)
  mat <- as.matrix(column_to_rownames(read_feather(count_matrix_file_path), var = "pos"))
  result <- rowcol_km_like_ComplexHeatmap(mat = mat, row_k = row_km, col_k = col_km)
  row_letter <- result$row_letter
  col_num <- result$col_num

  library(dplyr)
  df_row_out <- data.frame(
    feature = names(row_letter),
    label   = unname(row_letter),
    stringsAsFactors = FALSE
  )
  df_row_sorted <- df_row_out %>% arrange(label)

  df_col_out <- data.frame(
    feature = names(col_num),
    label   = unname(col_num),
    stringsAsFactors = FALSE
  )
  df_col_sorted <- df_col_out %>% arrange(label)

  if (is.null(count_matrix_file_path) || length(count_matrix_file_path) == 0 || any(is.na(count_matrix_file_path)) || any(!nzchar(count_matrix_file_path))) {
    stop("`count_matrix_file_path` is required", call. = FALSE)
  }

  if (!dir.exists(output_dir_path)) dir.create(output_dir_path, recursive = TRUE, showWarnings = FALSE)

  path1 <- file.path(output_dir_path, "row_table.tsv")
  path2 <- file.path(output_dir_path, "col_table.tsv")

  write.table(df_row_sorted, path1, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(df_col_sorted, path2, sep = "\t", quote = FALSE, row.names = FALSE)

  apply_cluster_heatmap(count_matrix_file_path = count_matrix_file_path, row_cluster_file_path = path1, col_cluster_file_path = path2, output_dir_path = output_dir_path)

}
