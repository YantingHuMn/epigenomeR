

biclustering_annotation_TFBS <- function(row_cluster_file_path, apply_heatmap = TRUE, out_dir = "./") {
    # Load packages
    suppressPackageStartupMessages({
        library(data.table)
        library(GenomicRanges)
        library(ComplexHeatmap)
        library(circlize)
    })

    row_cluster <- read.table(row_cluster_file_path, header = TRUE, sep = "\t")
    pos_df <- do.call(rbind, (strsplit(row_cluster$feature, "_")))
    colnames(pos_df) <- c("seqnames", "start", "end")
    row_cluster <- cbind(pos_df, row_cluster)
    row_gr <- makeGRangesFromDataFrame(row_cluster, seqnames.field = "seqnames", start.field = "start", end.field = "end", keep.extra.columns = TRUE)

    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }

    if (apply_heatmap) {

    }
}