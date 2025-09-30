map_target_names <- function(target_pair_list, target_pair_mapping_df_path = NULL, from = "targets", to = "shorthand" ) {
  if (is.null(target_pair_mapping_df_path)) {
    target_pair_mapping_df_path <- "/dcs05/hongkai/data/next_cutntag/script/utils/target_pair_short_hand.csv"
  }

  target_pair_mapping_df <- read.table(target_pair_mapping_df_path, header = TRUE, check.names = FALSE)
  cur_names <- target_pair_mapping_df[[from]]
  new_names <- target_pair_mapping_df[[to]]
  result <- target_pair_list
  for (i in seq_along(cur_names)) {
    cur_name <- cur_names[i]
    new_name <- new_names[i]
    result <- gsub(cur_name, new_name, result, fixed = TRUE)
  }
  return(result)
}

# Define a helper function to create a small dataframe for one group.
get_cluster_df <- function(df, prefix) {
  # Get the column names from the given dataframe.
  original_names <- colnames(df)

  # Create a dataframe with:
  # - "feature": prefixed column names, e.g. "V1:colName"
  # - "label": the prefix (group label)
  # - "nonprefix": the original column names (to help with ordering)
  data.frame(
    feature    = paste0(prefix, ":", original_names),
    label      = prefix,
    nonprefix  = original_names,
    stringsAsFactors = FALSE
  )
}

calc_ht_size1 = function(ht, unit = "inch", show_annotation_legend=FALSE, column_title=NULL, column_title_fontsize=14) {
  pdf(NULL)
  if (show_annotation_legend) {
    ht = draw(ht, background = "transparent", column_title=column_title, column_title_gp = gpar(fontsize=column_title_fontsize), merge_legend = TRUE, annotation_legend_side = "top")
  } else {
    ht = draw(ht, background = "transparent", column_title=column_title, column_title_gp = gpar(fontsize=column_title_fontsize), show_annotation_legend = FALSE)
  }
  ht = draw(ht, background = "transparent", column_title=column_title, column_title_gp = gpar(fontsize=column_title_fontsize))
  w = ComplexHeatmap:::width(ht)
  w = convertX(w, unit, valueOnly = TRUE)
  h = ComplexHeatmap:::height(ht)
  h = convertY(h, unit, valueOnly = TRUE)
  dev.off()

  c(w, h)
}

filter_target_pairs <- function(percentage_cutoff = 0.25, target_pairs=NULL, frag_len_num_file="/dcs05/hongkai/data/next_cutntag/bulk/frag_len/frag_split_fastq-demux_num_V_peak-valley.tsv") {
  frag_len_num <- read.table(frag_len_num_file, sep = "\t", header = TRUE, row.names = 1)
  frag_len_num$total_frag <- rowSums(frag_len_num)
  percentile_res <- quantile(frag_len_num$total_frag, type=3, probs = c(percentage_cutoff))
  # print(percentile_res)
  frag_len_num_cutoff <- unname(percentile_res)
  frag_len_num_filtered <- frag_len_num[frag_len_num$total_frag >= frag_len_num_cutoff, ]
  filtered_target_pairs <- rownames(frag_len_num_filtered)
  if (!is.null(target_pairs)) {
    filtered_target_pairs <- filtered_target_pairs[filtered_target_pairs %in% target_pairs]
  }
  return(filtered_target_pairs)
}

list_to_dataframe <- function(dataList) {
  if (is.null(names(dataList)))
    return(do.call('rbind', dataList))

  cn <- lapply(dataList, colnames) %>% unlist %>% unique
  cn <- c('.id', cn)
  dataList2 <- lapply(seq_along(dataList), function(i) {
    data = dataList[[i]]
    data$.id = names(dataList)[i]
    idx <- ! cn %in% colnames(data)
    if (sum(idx) > 0) {
      for (i in cn[idx]) {
        data[, i] <- NA
      }
    }
    return(data[,cn])
  })
  res <- do.call('rbind', dataList2)
  res$.id <- factor(res$.id, levels=rev(names(dataList)))
  return(res)
}
