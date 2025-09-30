# Post: Generate a tag × tag matrix from all_df and filter non-significant lower-triangle pairs based on filtered_df.
# Parameter: all_df: A data frame with 2 columns: pair name (like "YY1-cJun") and numeric value (e.g., peak count)
#            filtered_df: A data frame with at least 1 column: filtered significant pairs
#            group_csv: Optional data frame (not file path) with columns `tag_names` and `category_names`, defining tag groups and split order
#            tag_names: Column name in group_csv for tag identifiers
#            category_names: Column name in group_csv for grouping categories
#            by: Delimiter used in tag pairs (default = "-")
# Output: A list with:
#            - tags_peak_num_mat: tag × tag matrix, lower triangle filled with values or NA
#            - split: split vector (if group_csv is given), else NULL
#            - category_vector: unique group labels from group_csv, else NULL
generate_df <- function(all_df, filtered_df, group_csv = NULL, tag_names = "tag", category_names = "category",  by = "-") {
  # create original matrix
  if (!is.null(group_csv)) {
    group_df <- group_csv

    missing_cols <- setdiff(c(tag_names, category_names), colnames(group_df))
    if (length(missing_cols) > 0) {
      stop(glue::glue("Missing required column(s): {paste(missing_cols, collapse = ', ')} in csv file"))
    }

    csv_df <- group_df
    # sort tags
    csv_df <- csv_df[order(csv_df[[category_names]]), ]
    # add split to category
    csv_df$split <- as.integer(factor(csv_df$category))

    tags <- csv_df[[tag_names]]
    split <- csv_df$split
    category_vector <- unique(csv_df[[category_names]])
    tags_peak_num_mat <- matrix(0, nrow = length(tags), ncol = length(tags))
    rownames(tags_peak_num_mat) <- tags
    colnames(tags_peak_num_mat) <- tags
    colnames(tags_peak_num_mat) <- sub("\\..*$", "", colnames(tags_peak_num_mat))
    rownames(tags_peak_num_mat) <- sub("\\..*$", "", rownames(tags_peak_num_mat))
  } else {
    tags <- split_crf_pairs_to_single_crf(all_df$file, by = by)
    split <- NULL
    category_vector <- NULL
    tags_peak_num_mat <- matrix(0, nrow = length(tags), ncol = length(tags))
    rownames(tags_peak_num_mat) <- tags
    colnames(tags_peak_num_mat) <- tags
    colnames(tags_peak_num_mat) <- sub("\\..*$", "", colnames(tags_peak_num_mat))
    rownames(tags_peak_num_mat) <- sub("\\..*$", "", rownames(tags_peak_num_mat))
  }

  # import data
  for (i in seq_len(nrow(all_df))) {
    pair <- unlist(strsplit(all_df[i, 1], "-"))
    tag1 <- pair[1]
    tag2 <- pair[2]
    value <- df[i, 2]

    if (tag1 %in% tags && tag2 %in% tags) {
      index1 <- which(tags == tag1)
      index2 <- which(tags == tag2)

      if (index1 < index2) {
        tmp <- tag1
        tag1 <- tag2
        tag2 <- tmp
      }

      tags_peak_num_mat[tag1, tag2] <-  value
    }
  }

  valid_pairs <- unique(filtered_df[[1]])
  tags <- rownames(tags_peak_num_mat)
  for (i in seq_len(nrow(tags_peak_num_mat))) {
    for (j in seq_len(ncol(tags_peak_num_mat))) {
      if (i > j) {
        tag1 <- tags[i]
        tag2 <- tags[j]
        pair1 <- paste(tag1, tag2, sep = "-")
        pair2 <- paste(tag2, tag1, sep = "-")

        if (!(pair1 %in% valid_pairs || pair2 %in% valid_pairs)) {
          tags_peak_num_mat[i, j] <- NA
        }
      }
    }
  }

  return (list(tags_peak_num_mat = tags_peak_num_mat, split = split, category_vector = category_vector))
}
