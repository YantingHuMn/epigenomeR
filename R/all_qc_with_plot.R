# Post: Perform QC on BAM/BED files and generate a log2 peak count heatmap PDF.
# Parameters: file_path: A character vector of BAM or BED file paths for QC and visualization (required).
#             filtered_path: (Optional) A character vector of file paths for filtered QC matrix.
#             filtered_percentile: Numeric value between 0 and 1 to filter the lowest-read-count samples (default = 0.25).
#             output_path_dir: Directory to store output CSVs and PDF (default = same as first input file's directory).
#             split_crf_by: Delimiter used to split tag pairs into individual TF names (default = "-").
#             save: Logical; if TRUE, write `all_read_counts_*.csv`, `filtered_read_counts_*.csv`, and `total_reads.txt`.
#             group_csv: (Optional) A data frame defining tag groupings and categories for annotation.
#             tag_names: Column name in `group_csv` for tag identifiers (default = "tag").
#             category_names: Column name in `group_csv` for group categories (default = "category").
#             target_pair_mapping_df: (Optional) A data frame for renaming tag names in the heatmap matrix (default = NULL).
# Output: CSVs (if `save = TRUE`):
#            * all_read_counts_bam.csv or all_read_counts_bed.csv
#            * filtered_read_counts_bam.csv or filtered_read_counts_bed.csv
#            * total_reads.txt
#         Heatmap in PDF: AAA_test_qc_heatmap.pdf in `output_path_dir`
all_qc_with_plot <- function(file_path, filtered_path = NULL, filtered_percentile = 0.25, output_path_dir = NULL, split_crf_by = "-", save = TRUE, group_csv = NULL, tag_names = "tag", category_names = "category", target_pair_mapping_df = NULL) {
  if (is.null(output_path_dir)) {
    output_path_dir <- dirname(file_paths[1])
  }
  dir.create(output_path_dir, recursive = TRUE, showWarnings = FALSE)

  result_qc <- qc(file_paths = file_path, filtered_percentile = filtered_percentile, output_path_dir = output_path_dir, save = save)
  all_df <- result_qc$all_df
  filtered_df <- result_qc$filtered_df
  filtered_crf <- result_qc$filtered_crf
  total_reads <- result_qc$total_reads
  if (!is.null(filtered_path) && all(file.exists(filtered_path))) {
    result_qc <- qc(file_paths = filtered_path, filtered_percentile = filtered_percentile, output_path_dir = output_path_dir, save = save)
    filtered_df <- result_qc$filtered_df
    filtered_crf <- result_qc$filtered_crf
  }

  plot_peak_count_heatmap(all_df = all_df, filtered_df = filtered_df, group_csv = group_csv, by = split_crf_by, tag_names = tag_names, category_names = category_names, output_path_dir = output_path_dir, target_pair_mapping_df = target_pair_mapping_df)
}
