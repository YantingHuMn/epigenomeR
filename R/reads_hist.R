# Post: Generate and save histogram plot of read count distribution with optional vertical lines marking local minima for quality control analysis.
# Parameter: reads_vector: Numeric vector of read counts from BAM files
#            save_dir: Directory path to save the histogram plot
#            save_name: Filename for the saved histogram (should include .pdf extension)
#            local_min1: Optional first local minimum value to mark with red dashed line (default: NULL)
#            local_min2: Optional second local minimum value to mark with blue dotted line (default: NULL)
# Output: Saves histogram plot as PDF file to specified directory
reads_hist <- function(reads_vector, save_dir, save_name, local_min1 = NULL, local_min2 = NULL) {
  # load library
  suppressPackageStartupMessages({
    library(GenomicAlignments)
  })

  if (is.null(save_dir)) {
    stop("Error: Need output dir path: 'save_dir' is NULL.")
  }

  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }

  histDirectory <- file.path(save_dir, save_name)
  pdf(histDirectory)
  hist(reads_vector, breaks = 160)

  if (!is.null(local_min1)) {
    abline(v = local_min1, col = "red", lwd = 2, lty = 2)
  }
  if (!is.null(local_min2)) {
    abline(v = local_min2, col = "blue", lwd = 2, lty = 3)
  }
  dev.off()

  hist_data <- hist(reads_vector, breaks = 160, plot = FALSE)
}
