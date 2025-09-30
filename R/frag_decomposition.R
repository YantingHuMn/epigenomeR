# Post: Complete fragment length analysis pipeline including quality control, valley detection, and visualization for nucleosome positioning analysis from BAM files.
# Parameter: file_path: Vector of BAM file paths to analyze
#            save_dir: Directory path to save all output files and plots
#            filtered_percentile: Percentile threshold for quality control filtering (default: 0.25, removes bottom 25%)
#            dens_reso: Resolution for kernel density estimation in valley detection (default: 2^15)
#            target_pair_mapping_df: Optional data frame for mapping sample names to display names (default: NULL)
#            frag_decomp_file: Output parameter - path to generated fragment decomposition file
# Output: Generates comprehensive fragment analysis including histograms, valley detection, summary statistics, and bar plots; saves all results to specified directory
frag_decomposition <- function(file_path, save_dir, filtered_percentile = 0.25, dens_reso = 2^15, target_pair_mapping_df = NULL, frag_decomp_file) {
  # load library
  suppressPackageStartupMessages({
    library(GenomicAlignments)
  })

  result_qc <- qc(file_paths = file_path, filtered_percentile = filtered_percentile, save = FALSE)
  target_pairs_remained <- result_qc$filtered_crf

  bam_dir <- dirname(file_path[1])

  sample_dir_filename_list = c()
  for (sample_name in target_pairs_remained) {
    sample_filename = paste0(sample_name, ".bam")
    sample_dir_filename = file.path(bam_dir, sample_filename)
    sample_dir_filename_list = c(sample_dir_filename_list, sample_dir_filename)
  }
  bamWidths <- concatBams(sample_dir_filename_list, save_dir)
  print(paste0("After filter, the number of files remain: ", length(sample_dir_filename_list)))

  result <- reads_hist(reads_vector = bamWidths, save_dir = save_dir, save_name = "premerge_frag_hist.pdf")

  # find two loca minimum
  local_min_val <- find_two_global_valleys(y = bamWidths, dens_reso = dens_reso)

  summary_report <- data.frame(local_min1 = local_min_val[1], local_min2 = local_min_val[2], min_frag_len = min(bamWidths), max_frag_len = max(bamWidths))
  write.csv(summary_report, file.path(save_dir, "summary_report.csv"), row.names = FALSE)

  # plot again hist
  reads_hist(reads_vector = bamWidths, save_dir = save_dir, save_name = "premerge_frag_hist_with_cutoff.pdf", local_min1 = local_min_val[1], local_min2 = local_min_val[2])

  # df for plot
  result2 <- df_for_barplot_with_perc(bam_file_path = sample_dir_filename_list, save_dir = save_dir, valley1 = local_min_val[1], valley2 = local_min_val[2])
  frag_decomp_file <- result2$file_path
  # head(result2$df_for_barplot_with_perc)

  # plot bar
  frag_barplot(frag_decomp_file = frag_decomp_file, target_pair_list = target_pairs_remained, target_pair_mapping_df = target_pair_mapping_df, out_dir = paste0(save_dir, "/barplots"))

  # bam file decomposition in bash
  # finished
}
