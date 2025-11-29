
biclustering_wrapper <- function(cm_path, out_dir, apply_filter = TRUE, row_km = 15, col_km = 3, apply_annotation = TRUE) {
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
        mav_path <- mav_screen(path = merged_cm_path, out_dir = out_dir, plot=FALSE)
        f_cm_path <- informative_regions(input_hypervar_path=mav_path,input_CRF_orig_path=merged_cm_path,output_dir_path=out_dir, plot=FALSE)
    } else {
        f_cm_path <- merged_cm_path
    }

    # Step3: Biclustering & Annotation
    cat("\n","=="*10,"\n")
    cat("  Biclustering & Annotation")
    cat("=="*10,"\n")
    table_list <- apply_cluster_heatmap_advanced(count_matrix_file_path = f_cm_path, row_km = row_km, col_km = col_km,output_dir_path=out_dir)
    if (apply_annotation) {
        biclustering_annotation_ccre_hmm(row_cluster_file_path=table_list$row_table, output_dir_path=out_dir)
    }
}