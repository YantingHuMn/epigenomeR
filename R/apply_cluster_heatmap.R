# Apply cluster and draw heatmap
# Post: Apply combine_count_matrix to one or multiple count matrices
# Post: Wrapper function that processes single or multiple count matrices to generate clustered heatmaps. Handles data loading, validation, and optional matrix overlap analysis before calling the core heatmap generation function.
# Parameters: count_matrix_file_path: Vector of file paths to count matrix .feather files (rows=genomic positions, cols=samples)
#             row_cluster_file_path: Path to row cluster assignment .tsv file containing position-to-cluster mappings
#             col_cluster_file_path: Path to column cluster assignment .tsv file containing sample-to-cluster mappings
#             output_dir_path: Directory to save plot output (default: NULL, uses input file directory)
#             seed: Random seed for reproducible results (default: 123)
#             show_dend_boolean: Whether to show dendrograms in heatmap (default: FALSE)
#             count_matrix_overlap: Whether to merge shared positions across matrices for comparative analysis (default: FALSE)
#             lower_range: Lower bound for heatmap color scale (default: NULL, auto-determined)
#             upper_range: Upper bound for heatmap color scale (default: NULL, auto-determined)
#             row_title_fontsize: Font size for row cluster titles (default: NULL)
#             col_title_fontsize: Font size for column cluster titles (default: NULL)
#             legend_title_fontsize: Font size for legend title (default: NULL)
#             legend_label_fontsize: Font size for legend labels (default: NULL)
# Output: None (saves heatmap plots to specified directory)
apply_cluster_heatmap <- function(count_matrix_file_path, row_cluster_file_path, col_cluster_file_path, output_dir_path = NULL, seed = 123, show_dend_boolean = FALSE, count_matrix_overlap= FALSE, lower_range = NULL, upper_range = NULL, row_title_fontsize = NULL, col_title_fontsize = NULL, legend_title_fontsize = NULL, legend_label_fontsize = NULL) {
  # load library
  if (!requireNamespace("latex2exp", quietly = TRUE)) {
    stop("The latex2exp package is required but not installed.")
  }

  suppressPackageStartupMessages({
    library(ggplot2)
    library(ComplexHeatmap)
    library(circlize)
    library(tibble)
    library(arrow)
    library(latex2exp)
    library(svglite)
    library(viridis)
    library(glue)
    library(dplyr)
    library(tools)
  })

  for (i in seq_along(count_matrix_file_path)) {
    # load heatmap data -- count_matrix regions by CRF
    if (!file.exists(count_matrix_file_path[i])) {
      stop("The specified input file does not exist.")
    }
  }

  if (length(count_matrix_file_path)  == 1) {
    base_name <- basename(tools::file_path_sans_ext(count_matrix_file_path))
    if (is.null(output_dir_path) || output_dir_path == "") {
      output_dir_path <- dirname(count_matrix_file_path)
    } else {
      output_dir_path <- output_dir_path
    }
    V1V2_wgc = read_feather(count_matrix_file_path)
    if (!"pos" %in% colnames(V1V2_wgc)) {
      stop("The column 'pos' is missing from the input data.")
    }
    V1V2_wgc = column_to_rownames(V1V2_wgc, var="pos")
    combine_count_matrix(count_matrix = V1V2_wgc, row_cluster_file_path = row_cluster_file_path, col_cluster_file_path = col_cluster_file_path, basename = base_name, output_dir_path = output_dir_path, seed = seed, show_dend_boolean = show_dend_boolean, lower_range = lower_range, upper_range = upper_range, row_title_fontsize = row_title_fontsize, col_title_fontsize = col_title_fontsize, legend_title_fontsize = legend_title_fontsize, legend_label_fontsize = legend_label_fontsize)
  } else {
    if (count_matrix_overlap == FALSE) {
      for (i in seq_along(count_matrix_file_path)) {
        base_name <- basename(tools::file_path_sans_ext(count_matrix_file_path[i]))
        if (is.null(output_dir_path) || output_dir_path == "") {
          output_dir_path <- dirname(count_matrix_file_path[i])
        } else {
          output_dir_path <- output_dir_path
        }
        V1V2_wgc = read_feather(count_matrix_file_path[i])
        if (!"pos" %in% colnames(V1V2_wgc)) {
          stop("The column 'pos' is missing from the input data.")
        }
        V1V2_wgc = column_to_rownames(V1V2_wgc, var="pos")
        combine_count_matrix(count_matrix = V1V2_wgc, row_cluster_file_path = row_cluster_file_path, col_cluster_file_path = col_cluster_file_path, basename = base_name, output_dir_path = output_dir_path, seed = seed, show_dend_boolean = show_dend_boolean, lower_range = lower_range, upper_range = upper_range, row_title_fontsize = row_title_fontsize, col_title_fontsize = col_title_fontsize, legend_title_fontsize = legend_title_fontsize, legend_label_fontsize = legend_label_fontsize)
      }
    } else {
      count_matrix_list <- list()
      for (i in seq_along(count_matrix_file_path)) {
        mat <- read_feather(count_matrix_file_path[i])
        if (!"pos" %in% colnames(mat)) {
          stop(glue("The column 'pos' is missing from input {i}."))
        }
        mat <- column_to_rownames(mat, var = "pos")
        count_matrix_list[[i]] <- mat
      }
      shared_pos <- Reduce(intersect, lapply(count_matrix_list, rownames)) # no need to merge, only need overlapped row and col
      shared_cols <- Reduce(intersect, lapply(count_matrix_list, colnames))
      #count_matrix_merged <- do.call(cbind, lapply(count_matrix_list, function(x) x[shared_pos, , drop=FALSE]))
      for (i in seq_along(count_matrix_file_path)) {
        base_name <- paste0(basename(tools::file_path_sans_ext(count_matrix_file_path[i])), "all_combined")
        if (is.null(output_dir_path) || output_dir_path == "") {
          output_dir_path <- dirname(count_matrix_file_path[1])
        } else {
          output_dir_path <- output_dir_path
        }
        V1V2_wgc = read_feather(count_matrix_file_path[i])
        if (!"pos" %in% colnames(V1V2_wgc)) {
          stop("The column 'pos' is missing from the input data.")
        }
        V1V2_wgc = column_to_rownames(V1V2_wgc, var="pos")
        V1V2_wgc <- V1V2_wgc[shared_pos, shared_cols, drop = FALSE]
        combine_count_matrix(count_matrix = V1V2_wgc, row_cluster_file_path = row_cluster_file_path, col_cluster_file_path = col_cluster_file_path, basename = base_name, output_dir_path = output_dir_path, seed = seed, show_dend_boolean = show_dend_boolean, lower_range = lower_range, upper_range = upper_range, row_title_fontsize = row_title_fontsize, col_title_fontsize = col_title_fontsize, legend_title_fontsize = legend_title_fontsize, legend_label_fontsize = legend_label_fontsize)
      }
    }
  }
}
