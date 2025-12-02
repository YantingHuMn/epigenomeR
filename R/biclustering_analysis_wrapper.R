# Biclustering Analysis Wrapper
# Description: 
# 
# Parameters:
#   cm_path: Path to the count matrix `.feather` file or a vector of paths to multiple count matrix files to be merged.
#   out_dir: Directory to save all output files (cluster tables, filtered matrices, annotation plots).
#   apply_filter: Logical. Controls whether to further filter genomic regions.
#       - When the genome was segmented into equal-sized bins
#         (i.e., the count matrix was built using a numeric `regions` argument),
#         this should be TRUE so that low-information or uninformative bins can be removed.
#       - When the user supplied specific genomic intervals of interest
#         (i.e., the count matrix was built using a region file path),
#         this should be FALSE because no additional filtering is needed.
#   row_km: Number of k-means clusters for rows (genomic regions).
#   col_km: Number of k-means clusters for columns (CRF pairs).
#   apply_annotation: Logical. Controls whether to annotate genomic regions to nearby genes.
#       - When the genome was segmented into bins (numeric `regions`), this should be TRUE,
#         since bins lack inherent biological meaning and benefit from gene-level annotation.
#       - When the user provided specific regions of interest (region file path),
#         this should be FALSE, as those regions are already meaningful and do not require annotation.
#   plot: Logical. Whether to generate diagnostic plots during filtering and biclustering steps.
# 
# Output:

biclustering_analysis_wrapper <- function(cm_path, out_dir, apply_filter = TRUE, row_km = 15, col_km = 3, apply_annotation = TRUE, plot = TRUE) {
    # Step1: Merge all the count matrix files
    cat("\n","=="*10,"\n")
    cat("  Merge all count matrix files")
    cat("=="*10,"\n")
    merged_cm_path <- merge_count_matrices(cm_path = cm_path, out_dir = out_dir)

    # Step2: Filter highly variable regions
    cat("\n","=="*10,"\n")
    cat("  Filter highly variable regions")
    cat("=="*10,"\n")
    if (apply_filter) {
        f_cm_path <- detect_hvr(transformed_cm_path = merged_cm_path, out_dir = out_dir, plot = plot)
    } else {
        f_cm_path <- merged_cm_path
    }

    # Step3: Biclustering 
    cat("\n","=="*10,"\n")
    cat("  Biclustering")
    cat("=="*10,"\n")
    cluster_list <- biclustering(cm_path = f_cm_path, row_km = row_km, col_km = col_km, out_dir = out_dir, plot =  plot)

    # Step4: Biclustering annotation
    if (apply_annotation) {
        cat("\n","=="*10,"\n")
        cat("  Annotation")
        cat("=="*10,"\n")
        annotation_ccre_hmm(row_cluster_file_path = cluster_list$row_table, output_dir_path = out_dir)
    }
}