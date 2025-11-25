map_target_names <- function(target_pair_list, target_pair_mapping_df_path = NULL, from = "targets", to = "shorthand" ) {

  if (is.null(target_pair_mapping_df_path)) {
    datasets <- c("target_pair_mapping_df_path")
    for (d in datasets) {
      if (!exists(d)) {
        data(list = d, package = "epigenomeR")
      }
    }
    target_pair_mapping_df_path <- target_pair_mapping_df_path
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
