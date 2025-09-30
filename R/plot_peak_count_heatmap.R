# Post: Generate a heatmap (PDF) of tag-tag pair peak counts from a full matrix (filtered pairs in the lower triangle)
# Parameters: all_df: Data frame with all tag pairs and associated values (2 columns: "pair", "count").
#             filtered_df: Data frame of filtered tag pairs (1st column: "pair") to highlight in the heatmap.
#             group_csv: Optional data frame specifying tag grouping info. Must include columns `tag_names` and `category_names`.
#             by: Delimiter used to split tag pairs (default = "-").
#             tag_names: Column name in `group_csv` that holds tag names (default = "tag").
#             category_names: Column name in `group_csv` that holds grouping categories (default = "category").
#             output_path_dir: Directory to save the output heatmap PDF.
#             target_pair_mapping_df: (Optional) A mapping data frame to rename row/column names (not used if NULL).
# Output: PDF file saved to: <output_path_dir>/AAA_test_qc_heatmap.pdf
plot_peak_count_heatmap <- function(all_df, filtered_df, group_csv = NULL, by = "-", tag_names = "tag", category_names = "category", output_path_dir, target_pair_mapping_df = NULL) {
  # load library
  suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(latex2exp)
    library(circlize)
    library(grid)
  })

  # read in
  if (!is.null(output_path_dir)) {
    if (!dir.exists(output_path_dir)) {
      dir.create(output_path_dir, recursive = TRUE)
    }
    output_pdf_name <- file.path(output_path_dir, "AAA_test_qc_heatmap.pdf")
  }

  if (!is.null(group_csv)) {
    results <- generate_df(all_df = all_df, filtered_df = filtered_df, group_csv = group_csv,tag_names = tag_names, category_names = category_names, by = by)
    peak_count_mat <- results$tags_peak_num_mat
    category_vector <- results$category_vector
    splits <- results$split

    # if (!is.null(target_pair_mapping_df)) {
    #     rownames(peak_count_mat) <- map_target_names(rownames(peak_count_mat), target_pair_mapping_df)
    #     colnames(peak_count_mat) <- map_target_names(colnames(peak_count_mat), target_pair_mapping_df)
    # }

    na_legend <- Legend(labels = "Filtered Targets", legend_gp = gpar(fill = "grey90"), title = "")
    group_col <- c("#8ECFC9", "#FFBE7A", "#FA7F6F", "#82B0D2", "#BEB8DC")

    row_annotation = rowAnnotation(
      foo = anno_block(gp = gpar(fill = group_col),
                       labels = category_vector,
                       labels_gp = gpar(col = "white", fontsize = 11))
    )
    col_annotation = HeatmapAnnotation(
      foo = anno_block(gp = gpar(fill = group_col),
                       labels = category_vector,
                       labels_gp = gpar(col = "white", fontsize = 11))
    )

    ht_opt$ROW_ANNO_PADDING = unit(0.4, "cm")
    ht_opt$COLUMN_ANNO_PADDING = unit(0.4, "cm")

    log_mat <- log2(peak_count_mat + 1)
    log_min <- min(log_mat, na.rm = TRUE)
    log_max <- max(log_mat, na.rm = TRUE)
    log_avg <- (log_min + log_max)/2

    pdf(output_pdf_name, width = 12, height = 12)
    ht <- Heatmap(log_mat, rect_gp = gpar(type = "none"), na_col = "grey90",
                  cluster_rows = FALSE, cluster_columns = FALSE,
                  cell_fun = function(j, i, x, y, w, h, fill) {
                    if (i >= j) {
                      grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                    }
                  },
                  row_labels = TeX(rownames(peak_count_mat)),
                  column_labels = TeX(colnames(peak_count_mat)),
                  row_names_side = "left",
                  column_title = "Target Pair Peak Count",
                  #column_names_rot = 45,
                  col = colorRamp2(c(log_min, log_avg, log_max), c("blue", "white", "red")),
                  heatmap_legend_param = list(
                    title = "log2(peak counts)",
                    color_bar = "continuous",
                    title_gp = gpar(fontsize = 15)
                  ),
                  column_split = splits,
                  row_split = splits,
                  row_title = rep("", length(unique(splits))),
                  left_annotation = row_annotation,
                  bottom_annotation = col_annotation
    )
    draw(ht, annotation_legend_list = list(na_legend), merge_legend = TRUE)
    dev.off()
  } else {
    results <- generate_df(all_df_path = all_df_path, filtered_df_path = filtered_df_path, group_csv = group_csv,tag_names = tag_names, category_names = category_names, by = by)
    peak_count_mat <- results$tags_peak_num_mat

    # if (!is.null(target_pair_mapping_df)) {
    #     rownames(peak_count_mat) <- map_target_names(rownames(peak_count_mat), target_pair_mapping_df)
    #     colnames(peak_count_mat) <- map_target_names(colnames(peak_count_mat), target_pair_mapping_df)
    # }

    na_legend <- Legend(labels = "Filtered Targets", legend_gp = gpar(fill = "grey90"), title = "")
    group_col <- c("#8ECFC9", "#FFBE7A", "#FA7F6F", "#82B0D2", "#BEB8DC")

    ht_opt$ROW_ANNO_PADDING = unit(0.4, "cm")
    ht_opt$COLUMN_ANNO_PADDING = unit(0.4, "cm")

    log_mat <- log2(peak_count_mat + 1)
    log_min <- min(log_mat, na.rm = TRUE)
    log_max <- max(log_mat, na.rm = TRUE)
    log_avg <- (log_min + log_max)/2

    pdf(output_pdf_name, width = 12, height = 12)
    ht <- Heatmap(log_mat, rect_gp = gpar(type = "none"), na_col = "grey90",
                  cluster_rows = FALSE, cluster_columns = FALSE,
                  cell_fun = function(j, i, x, y, w, h, fill) {
                    if (i >= j) {
                      grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                    }
                  },
                  row_labels = TeX(rownames(peak_count_mat)),
                  column_labels = TeX(colnames(peak_count_mat)),
                  row_names_side = "left",
                  column_title = "Target Pair Peak Count",
                  #column_names_rot = 45,
                  col = colorRamp2(c(log_min, log_avg, log_max), c("blue", "white", "red")),
                  heatmap_legend_param = list(
                    title = "log2(peak counts)",
                    color_bar = "continuous",
                    title_gp = gpar(fontsize = 15)
                  )#,row_title = rep("", length(unique(splits)))
    )
    draw(ht, annotation_legend_list = list(na_legend), merge_legend = TRUE)
    dev.off()

  }

}
