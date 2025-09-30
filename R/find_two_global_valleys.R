# Post: Find two global valleys in fragment length distribution using kernel density estimation to identify nucleosome positioning patterns in ATAC-seq or similar data.
# Parameter: y: Numeric vector of fragment lengths from sequencing data
#            dens_reso: Resolution for kernel density estimation (default: 2^15 = 32768 points)
#            density_kernel: Kernel type for density estimation (default: "gaussian")
#            valley1_range: Range to search for first valley - nucleosome-free region (default: c(73, 221))
#            valley2_range: Range to search for second valley - dinucleosome region (default: c(222, 368))
# Output: Named numeric vector with local_min1 and local_min2 positions (rounded to integers)
find_two_global_valleys <- function(y, dens_reso = 2^15, density_kernel = "gaussian", valley1_range = c(73, 221), valley2_range = c(222, 368)) {

  library(sfsmisc)

  if (!is.numeric(y) || length(y) == 0) {
    stop("Input 'y' must be a non-empty numeric vector")
  }

  if (dens_reso <= 0 || !is.numeric(dens_reso)) {
    stop("dens_reso must be a positive number")
  }

  # Compute kernel density estimation
  density_result <- density(y, kernel = density_kernel, n = dens_reso)
  hist_x <- density_result$x
  hist_y <- density_result$y

  # Extract regions for valley analysis
  # Valley 1 region (typically around nucleosome-free fragments)
  valley1_mask <- hist_x >= valley1_range[1] & hist_x <= valley1_range[2]
  valley1_x_region <- hist_x[valley1_mask]
  valley1_y_region <- hist_y[valley1_mask]

  # Valley 2 region (typically around dinucleosome fragments)
  valley2_mask <- hist_x >= valley2_range[1] & hist_x <= valley2_range[2]
  valley2_x_region <- hist_x[valley2_mask]
  valley2_y_region <- hist_y[valley2_mask]

  # Find valleys (minima) in each region
  if (length(valley1_y_region) == 0) {
    warning("No data points found in valley1_range")
    valley1_x <- NA
  } else {
    valley1_min_idx <- which.min(valley1_y_region)
    valley1_x <- valley1_x_region[valley1_min_idx]
  }

  if (length(valley2_y_region) == 0) {
    warning("No data points found in valley2_range")
    valley2_x <- NA
  } else {
    valley2_min_idx <- which.min(valley2_y_region)
    valley2_x <- valley2_x_region[valley2_min_idx]
  }
  left_rough <- round(valley1_x)
  right_rough <- round(valley2_x)

  # refine_valley

  return(c(local_min1 = left_rough, local_min2 = right_rough))
}
