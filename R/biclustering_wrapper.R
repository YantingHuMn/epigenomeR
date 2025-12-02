
# Parameters:
#   apply_filter: Logical. Controls whether to further filter genomic regions.
#       - When the genome was segmented into equal-sized 
#         (i.e., binsthe count matrix was built using a numeric `regions` argument ), this should be TRUE 
#         so that low-information or uninformative bins can be removed.
#       - When the user supplied specific genomic intervals of interestt
#         (i.e., he count matrix was built using a region file path ), 
#         this should be FALSE because no additional filtering is needed.


biclustering_wrapper <- function(cm_path, out_dir, apply_filter = TRUE, row_km = 15, col_km = 3, apply_annotation = TRUE, plot = TRUE) {
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
    cat("\n","=="*10,"\n")
    cat("  Annotation")
    cat("=="*10,"\n")
    if (apply_annotation) {
        biclustering_annotation_ccre_hmm(row_cluster_file_path=table_list$row_table, output_dir_path = out_dir)
    }
}