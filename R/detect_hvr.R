# Filter highly variable regions
# Description:
#   This is a convenience wrapper function that combines mean-variance screening 
#   and highly variable region (HVR) selection into a single step. It automatically:
#   1. Performs MAV (Mean-Variance) screening to calculate hypervariance metrics
#   2. Filters regions based on hypervariance and expression levels
#   3. Returns a filtered count matrix containing only informative regions
# 
# Parameters:
#   transformed_cm_path: Path to transformed count matrix file (.feather format with "pos" column)
#   out_dir: Output directory path (default: "./")
#   plot: Whether to generate diagnostic scatter plots (default: FALSE)
# 
# Output Files:
#   1. <input>_mav_screen.feather - Mean-variance statistics (intermediate)
#   2. <input>_filtered_regions.feather - Final filtered count matrix (main output)
#   3. <input>_mav_screen.png - MAV diagnostic plot (if plot=TRUE)
#   4. <input>_filtered_regions.png - HVR selection plot (if plot=TRUE)
# Returns:
#   Character string: File path to the filtered count matrix (.feather format)

detect_hvr <- function(transformed_cm_path, out_dir = "./", plot = FALSE) {
    mav_path <- mav_screen(transformed_cm_path = transformed_cm_path, out_dir = out_dir, plot = plot)
    filtered_cm_path <- filter_hvr(transformed_cm_path = transformed_cm_path, mav_stats_path = mav_path, out_dir = out_dir, plot = plot)
    return(filtered_cm_path)
}