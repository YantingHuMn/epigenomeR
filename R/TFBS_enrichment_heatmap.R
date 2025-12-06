# Peak Enrichment Heatmap
# Post: Generate heatmaps of motif enrichment (log-odds ratio) from *_result.tsv files.
# Parameter: motif_result_dir: tsv file's directory
#            out_dir: outpur directory
#            top_n: Number of top enriched motifs (per column) to display in the heatmap.
#            selected_tfs: Vector of transcription factor names or substrings to filter rows (motifs).
#            transformations: List of transformations to apply in order: log2, t, Z, cluster
# Output: Heatmap PDF (full matrix): peak_no_perturb_<transformations>.pdf
#         Heatmap PDF (top_n motifs): reorganized_heatmap_filtered_rna_subtraction_<top_n>_<transformations>.pdf
#         log2 odds ratio matrix: peak_no_perturb.csv
#         FDR matrix:  peak_no_perturb_fdr.csv
motif_enrichment_heatmap <- function(motif_result_dir, out_dir, top_n = NULL, selected_tfs = NULL, transformations = c("log2")) {
  # load library
  suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(plyranges)
    library(rtracklayer)
    library(matrixTests)
    library(latex2exp)
    library(glue)
    library(circlize)
  })

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  motif_result_files <- Sys.glob(glue("{motif_result_dir}/*_result.tsv"))
  motif_results <- list()
  folders <- c()

  for (f in motif_result_files) {
    folder <- tools::file_path_sans_ext(basename(f))
    folder <- gsub("_result$", "", folder)
    df <- read.table(f, sep="\t", header=TRUE, row.names=1)
    motif_results[[folder]] <- df
    folders <- c(folders, folder)
  }

  # arrange folders
  suppressWarnings({
    test_numeric <- suppressWarnings(as.numeric(folders))
    if (all(!is.na(test_numeric))) {
      folders <- as.character(sort(test_numeric))
    } else {
      folders <- sort(folders)
    }
  })
  motif_results <- motif_results[folders]

  motif_ids <- rownames(motif_results[[1]])
  odds_ratio_df <- as.data.frame(lapply(motif_results, function(x) x[motif_ids, "odds_ratio"]), check.names = FALSE)
  FDR_df <- as.data.frame(lapply(motif_results, function(x) x[motif_ids, "FDR"]), check.names = FALSE)
  target_hit_df <- as.data.frame(lapply(motif_results, function(x) x[motif_ids, "target_hit"]), check.names = FALSE)
  rownames(odds_ratio_df) <- rownames(FDR_df) <- rownames(target_hit_df) <- motif_ids

  target_hit_thresh <- quantile(unlist(target_hit_df), 0.1)
  keep <- Reduce(intersect, list(
    rownames(odds_ratio_df[rowSums(odds_ratio_df >= 2) >= 1, ]),
    rownames(FDR_df[rowSums(FDR_df <= 0.05) >= 1, ]),
    rownames(target_hit_df[rowSums(target_hit_df >= target_hit_thresh) >= 1, ])
  ))
  odds_ratio_df <- odds_ratio_df[keep, ]
  FDR_df <- FDR_df[keep, ]

  non_perturb <- !grepl("(", rownames(odds_ratio_df), fixed = TRUE)
  odds_ratio_df <- odds_ratio_df[non_perturb, ]
  FDR_df <- FDR_df[rownames(odds_ratio_df), ]

  if (!is.null(selected_tfs)) {
    selected_rows <- unique(unlist(lapply(selected_tfs, function(tf) {
      grep(tf, rownames(odds_ratio_df), value = TRUE)
    })))
    odds_ratio_df <- odds_ratio_df[selected_rows, , drop = FALSE]
    FDR_df <- FDR_df[selected_rows, , drop = FALSE]
  }

  # heatmap
  # all
  # transformation
  transformed_df <- transform_df(odds_ratio_df, transformations = transformations)
  col_fun <- colorRamp2(c(-max(abs(transformed_df)), 0, max(abs(transformed_df))), c("#3155C3", "white", "#AF0525"))

  # name
  preferred_order <- c("log2", "t", "Z", "cluster")
  transformations_sorted <- preferred_order[preferred_order %in% transformations]
  names <- paste(transformations_sorted, collapse = "_")

  if ("cluster" %in% transformations) {
    h <- Heatmap(transformed_df, name = names, col = col_fun, cluster_columns = TRUE)
    ht_drawn <- draw(h)
    col_order <- column_order(ht_drawn)
    ordered_colnames <- rev(colnames(transformed_df)[col_order])
    h <- Heatmap(transformed_df[,ordered_colnames], name = names, col = col_fun, cluster_columns = FALSE)
    V1V2_size <- calc_ht_size_heatmap((transformed_df[,ordered_colnames]))
    pdf(glue("{out_dir}/peak_no_perturb_{names}.pdf"), width = V1V2_size[1], height = V1V2_size[2])
    draw(h)
    dev.off()
  } else {
    h <- Heatmap(transformed_df, name = names, col = col_fun, cluster_columns = FALSE)
    V1V2_size <- calc_ht_size_heatmap(transformed_df)
    pdf(glue("{out_dir}/peak_no_perturb_{names}.pdf"), width = V1V2_size[1], height = V1V2_size[2])
    draw(h)
    dev.off()
  }

  # top n
  if (!is.null(top_n) && is.numeric(top_n) && top_n > 0) {
    selected_motif_list <- list()
    selected_motifs <- c()
    for (folder in folders) {
      motif_odds_ratio_mat_sub <- odds_ratio_df
      odds_ratio_i <- motif_odds_ratio_mat_sub[, folder, drop=FALSE]
      odds_ratio_i <- odds_ratio_i[order(odds_ratio_i[[folder]], decreasing = TRUE),,drop=FALSE]
      selected_motifs <- append(selected_motifs, rownames(odds_ratio_i)[1:top_n])
      selected_motif_list[[folder]] <- rownames(odds_ratio_i)[1:top_n]
    }
    selected_motifs <- unique(selected_motifs)
    odds_ratio_df <- odds_ratio_df[selected_motifs,]

    # transformation
    transformed_df <- transform_df(odds_ratio_df, transformations = transformations)
    col_fun <- colorRamp2(c(-max(abs(transformed_df)), 0, max(abs(transformed_df))), c("#3155C3", "white", "#AF0525"))

    # name
    preferred_order <- c("log2", "t", "Z", "cluster")
    transformations_sorted <- preferred_order[preferred_order %in% transformations]
    names <- paste(transformations_sorted, collapse = "_")

    if ("cluster" %in% transformations) {
      h <- Heatmap(transformed_df, name = names, col = col_fun, cluster_columns = TRUE)
      ht_drawn <- draw(h)
      col_order <- column_order(ht_drawn)
      ordered_colnames <- rev(colnames(transformed_df)[col_order])
      h <- ComplexHeatmap::Heatmap(transformed_df[,ordered_colnames],
                                   column_names_gp = grid::gpar(fontsize = 25),
                                   row_names_gp = grid::gpar(fontsize = 25),
                                   cluster_columns = F,
                                   cluster_rows = F,
                                   # column_names_rot = 45,
                                   row_names_side = "left",
                                   col = col_fun,
                                   heatmap_legend_param = list(
                                     title = "log(Odds Ratio)",
                                     color_bar = "continuous",
                                     title_gp = gpar(fontsize = 24),
                                     legend_direction = "horizontal",
                                     labels_gp = gpar(fontsize = 20) ,
                                     legend_width = unit(6, "cm")
                                   ))
      V1V2_size <- calc_ht_size_heatmap(transformed_df[,ordered_colnames])
      pdf(glue("{out_dir}/reorganized_heatmap_filtered_rna_subtraction_{top_n}_{names}.pdf"), width = V1V2_size[1], height = V1V2_size[2])
      draw(h, heatmap_legend_side="bottom", padding = unit(c(8, 8, 8, 8), "mm"))
      dev.off()
    } else {
      h <- ComplexHeatmap::Heatmap(transformed_df,
                                   column_names_gp = grid::gpar(fontsize = 25),
                                   row_names_gp = grid::gpar(fontsize = 25),
                                   cluster_columns = F,
                                   cluster_rows = F,
                                   # column_names_rot = 45,
                                   row_names_side = "left",
                                   col = col_fun,
                                   heatmap_legend_param = list(
                                     title = "log(Odds Ratio)",
                                     color_bar = "continuous",
                                     title_gp = gpar(fontsize = 24),
                                     legend_direction = "horizontal",
                                     labels_gp = gpar(fontsize = 20) ,
                                     legend_width = unit(6, "cm")
                                   ))
      V1V2_size <- calc_ht_size_heatmap(transformed_df)
      pdf(glue("{out_dir}/reorganized_heatmap_filtered_rna_subtraction_{top_n}_{names}.pdf"), width = V1V2_size[1], height = V1V2_size[2])
      draw(h, heatmap_legend_side="bottom", padding = unit(c(8, 8, 8, 8), "mm"))
      dev.off()
    }
  }
  # save
  write.csv(log2_df, glue("{out_dir}/peak_no_perturb.csv"))
  write.csv(FDR_df, glue("{out_dir}/peak_no_perturb_fdr.csv"))
}
