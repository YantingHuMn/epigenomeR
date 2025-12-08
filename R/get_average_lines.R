get_average_lines <- function(pdp_results, target_pair) {
  all_arrays <- list()=
    for (random_seed in names(pdp_results)) {
      individual_lines <- pdp_results[[random_seed]][["coding_all"]][[target_pair]][["individual"]]
      all_arrays <- append(all_arrays, list(individual_lines[[1]]))
    }
  stacked <- do.call(rbind, all_arrays)
  return(stacked)
}
