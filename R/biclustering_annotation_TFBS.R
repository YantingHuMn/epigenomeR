

biclustering_annotation_TFBS <- function(row_cluster_file_path, out_dir = "./", ref_genome = "hg38", regions = 800, plot = TRUE) {
    # Load packages
    suppressPackageStartupMessages({
        library(data.table)
        library(GenomicRanges)
        library(ComplexHeatmap)
        library(circlize)
        library(rtracklayer)
    })

    row_cluster <- read.table(row_cluster_file_path, header = TRUE, sep = "\t")
    pos_df <- do.call(rbind, (strsplit(row_cluster$feature, "_")))
    colnames(pos_df) <- c("seqnames", "start", "end")
    row_cluster <- cbind(pos_df, row_cluster)
    row_gr <- makeGRangesFromDataFrame(row_cluster, seqnames.field = "seqnames", start.field = "start", end.field = "end", keep.extra.columns = TRUE)
    style <- seqlevelsStyle(row_gr)[1]
    row_grl <- split(row_gr, row_gr$label)

    control_grl <- lapply(names(row_grl), function(label) {
        cat("Processing cluster:", label, "with", length(row_grl[[label]]), "regions\n")
        target_gr <- row_grl[[lab]]
        control_gr <- get_matched_control(target_gr = target_gr, ref_genome = ref_genome, style = style, n_rep = 1, regions = regions)
        return(control_gr)
    })
    names(control_grl) <- names(row_grl)
    control_gr <- do.call(c, control_grl)

    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }
    export(control_gr, file.path(out_dir, "all_controls.bed"))

    tsv_path <- TFBS_enrichment(target_region = row_gr, control_region = control_gr, regions = regions, out_path = file.path(out_dir, "TFBS_enrichment.tsv"), ref_genome = ref_genome, style = style)

    if (apply_heatmap) {
        TFBS_enrichment_heatmap()
    }
}

