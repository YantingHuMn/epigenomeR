# Post: Extract unique TF names or components from a vector of pair strings like "YY1-cJun"
# Parameter: pairs: A character vector of strings representing pairs, e.g., "YY1-cJun"
#            by: Delimiter used to separate pairs (default = "-")
# Output: A character vector of unique component names. (exclude values like "is", "0", "200")
split_crf_pairs_to_single_crf <- function(pairs, by = "-") {
  valid_pairs <- pairs[grepl(by, pairs)]
  all_names <- unique(unlist(strsplit(valid_pairs, by)))
  all_names <- all_names[grepl("^[A-Za-z]", all_names)]
  all_names <- setdiff(all_names, c("is", "0", "200"))

  return (all_names)
}
