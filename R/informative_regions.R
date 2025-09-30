# Select Informative Regions
# Post: Filters informative regions based on log2(mean) and hypervar values,
# Parameter: input_hypervar_path:  HyperVar Matrix path.
#            input_CRF_orig_path: Count Matrix path.
#            output_dir_path: Output directory path.
#            split_num: The number of splits to divide log2(mean) values into for filtering (default is 100).
#            keep_percent_list: List of percentages to keep from each split based on the top informative regions (default is 0.01).
#            log2mean_quantile_thres: Quantile threshold for log2(mean) values to filter the most informative regions (default is 0.99).
#            plot: Whether to generate and save plots (default is FALSE).
# Output: A filtered matrix of informative regions saved as a .feather file.
#         Two scatter plots saved as .PNG (log2(mean) vs. hypervar, log2(mean) vs. log2(var))
informative_regions <- function(input_hypervar_path, input_CRF_orig_path, output_dir_path = NULL, split_num = 100, keep_percent_list = c(0.01), log2mean_quantile_thres = 0.99, plot=FALSE) {
  # load library
  suppressPackageStartupMessages({
    library(feather)
    library(ggplot2)
    library(dplyr)
    library(reshape2)
    library(cowplot)
    library(arrow)
    library(tibble)
  })

  if (is.null(output_dir_path) || output_dir_path == "") {
    output_dir <- dirname(input_hypervar_path)
  } else {
    output_dir <- output_dir_path
  }

  hypervar_summary <- arrow::read_feather(input_hypervar_path)

  hypervar_summary$log2_mean <- log2(hypervar_summary$mean)
  hypervar_summary$log2_var <- log2(hypervar_summary$var)
  mean_min <- min(hypervar_summary$log2_mean)
  mean_max <- max(hypervar_summary$log2_mean)
  log2mean_thres <- quantile(hypervar_summary$log2_mean, log2mean_quantile_thres)
  print(log2mean_thres)

  mean_range_list <- seq(mean_min, mean_max, length.out = split_num + 1)

  # Select regions based on thresholds
  for (keep_percent in keep_percent_list) {
    keep_num <- nrow(hypervar_summary) * keep_percent
    keep_num_split <- floor(keep_num / split_num)

    envelop_select_list <- list()
    for (split_idx in 1:split_num) {
      mean_ub <- mean_range_list[split_idx]
      mean_lb <- mean_range_list[split_idx + 1]

      if (split_idx == split_num) {
        hypervar_summary_snap <- subset(hypervar_summary, log2_mean >= mean_ub & log2_mean <= mean_lb & hypervar > 1)
      } else {
        hypervar_summary_snap <- subset(hypervar_summary, log2_mean >= mean_ub & log2_mean < mean_lb & hypervar > 1)
      }

      hypervar_summary_snap_sort <- hypervar_summary_snap %>% arrange(desc(hypervar))
      rownames(hypervar_summary_snap_sort) <- NULL

      hypervar_summary_snap_sort_keep <- hypervar_summary_snap_sort[1:keep_num_split, ]
      envelop_select_list[[split_idx]] <- hypervar_summary_snap_sort_keep
    }
  }
  # Combine selected data
  envelop_select <- bind_rows(envelop_select_list)
  evenlop_select_keep <- envelop_select %>% filter(log2_mean >= log2mean_thres)

  # Load transformed count matrix file
  wgc_raw <- arrow::read_feather(input_CRF_orig_path)
  wgc <- wgc_raw %>% tibble::column_to_rownames("pos")
  wgc_evenlop_select_keep <- wgc[rownames(wgc) %in% evenlop_select_keep$pos, ]

  # save
  input_name <- basename(tools::file_path_sans_ext(input_hypervar_path))
  input_dir_dir <- dirname(input_hypervar_path)

  wgc_evenlop_select_keep_for_save <- wgc_evenlop_select_keep %>% tibble::rownames_to_column("pos")
  wgc_evenlop_select_keep_filename <- paste0(input_name, "_filtered_regions.feather")
  wgc_evenlop_select_keep_dir_filename <- file.path(output_dir, wgc_evenlop_select_keep_filename)
  write_feather(wgc_evenlop_select_keep_for_save, wgc_evenlop_select_keep_dir_filename)

  # Plots
  if (plot == TRUE) {
    # Plot scatter plot
    ggplot() +
      geom_point(data = hypervar_summary, aes(x = log2_mean, y = hypervar),
                 alpha = 0.3, color = "blue") +
      geom_point(data = envelop_select, aes(x = log2_mean, y = hypervar),
                 alpha = 0.3, color = "red") +
      geom_point(data = evenlop_select_keep, aes(x = log2_mean, y = hypervar),
                 size = 3, color = "green", alpha = 0.8) +
      labs(x = "log2(mean)", y = "log2(var)") +
      theme_minimal()

    # Saving the plot
    fig_wgc_evenlop_select_keep_filename <- paste0(input_name, "_filtered_regions_plot1.png")
    fig_wgc_evenlop_select_keep_dir_filename <- file.path(output_dir, fig_wgc_evenlop_select_keep_filename)
    ggsave(fig_wgc_evenlop_select_keep_dir_filename, dpi = 600)  # Saving the plot with 600 DPI

    # Plot joint plot
    ggplot(hypervar_summary, aes(x = log2_mean, y = log2_var)) +
      geom_point(alpha = 0.3, color = "#D4E9D9") +
      geom_point(data = evenlop_select_keep, aes(x = log2_mean, y = log2_var),
                 size = 3, color = "#F6E0D3", alpha = 0.8) +
      labs(x = "log2(mean)", y = "log2(var)") +
      theme_minimal()

    # Saving the plot
    fig_wgc_evenlop_select_keep_filename <- paste0(input_name, "_filtered_regions_plot2.png")
    fig_wgc_evenlop_select_keep_dir_filename <- file.path(output_dir, fig_wgc_evenlop_select_keep_filename)
    ggsave(fig_wgc_evenlop_select_keep_dir_filename, dpi = 600)  # Saving the plot with 600 DPI

    # Plot heatmap (clustermap)
    # cluster_plot_clip <- Heatmap(wgc_evenlop_select_keep, name = "cluster_plot",
    #                             cluster_columns = TRUE, cluster_rows = TRUE,
    #                             show_column_names = TRUE, show_row_names = FALSE,
    #                             heatmap_legend_param = list(title = "Correlation"))

  }
}
