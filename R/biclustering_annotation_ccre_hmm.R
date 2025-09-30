# Post: Annotate and Plot CCRE and ChromHMM Composition for Row Cluster
# Post: Analyzes genomic annotation composition of clustered regions by overlapping with CCRE and ChromHMM annotations, then generates comparative bar plots showing regulatory element distribution across clusters.
# Parameter: row_cluster_file_path: Path to row cluster assignment .tsv file with feature positions and cluster labels
#            output_dir_path: Directory to save annotation plot outputs
# Output: Generates and saves two bar plot PDF files for CCRE and ChromHMM composition analysis
biclustering_annotation_ccre_hmm <- function(row_cluster_file_path, output_dir_path) {
  # load library
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("The cowplot package is required but not installed.")
  }

  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("The GenomicRanges package is required but not installed.")
  }

  # if (!requireNamespace("RJSONIO", quietly = TRUE)) {
  #     stop("The RJSONIO package is required but not installed.")
  # }

  # if (!requireNamespace("rtracklayer", quietly = TRUE)) {
  #     stop("The rtracklayer package is required but not installed.")
  # }

  suppressPackageStartupMessages({
    library(ggplot2)
    library(cowplot)
    library(GenomicRanges)
    library(glue)
    #library(RJSONIO)
    library(tidyverse)
    #library(rtracklayer)
  })

  source("/dcl02/hongkai/data/kyu/multitag_scripts/scripts_peakAnnotation/ccre_annotation/CCREUtils.R")
  source("/dcs05/hongkai/data/next_cutntag/script/utils/utils.R")
  source("/dcs05/hongkai/data/next_cutntag/script/utils/filter_targets.R")

  # split pos
  biclustering_result <- read.table(row_cluster_file_path, header = TRUE, sep = "\t")
  biclustering_result <- biclustering_result[!biclustering_result$label %in% c("Epitope_Specific", "Background"), ]
  pos_df <- do.call(rbind, (strsplit(biclustering_result$feature, "_")))
  colnames(pos_df) <- c("seqnames", "start", "end")
  biclustering_result <- cbind(pos_df, biclustering_result)
  biclustering_gr <- makeGRangesFromDataFrame(
    biclustering_result,
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE
  )

  cluster_all_levels <- c(
    as.character(1:100),
    "A",
    "B",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "J",
    "K",
    "L",
    "M",
    "N",
    "O",
    "Epitope_Specific",
    "non-target_specific",
    "background",
    "Background"
  )

  unique_clusters <- sort(unique(biclustering_result$label), method = "radix")
  chipseeker_ccre_annotations <- list()
  chipseeker_ccre_celltype_agnostic_annotations <- list()
  repeatmasker_annotations <- list()
  chromhmm_short_annotations <- list()
  chromhmm_full_annotations <- list()
  chromhmm_group_annotations <- list()
  for (cluster_id in unique_clusters) {
    biclustering_result_i <- biclustering_result[biclustering_result$label == cluster_id, ]
    biclustering_grange_i <- makeGRangesFromDataFrame(biclustering_result_i, seqnames.field = "seqnames", start.field = "start", end.field = "end", keep.extra.columns = TRUE)
    #####
    chipseeker_ccre_annotation <- annotatePeakByOverlappingChIPSeekerCCRE(biclustering_grange_i, annotation, categories)
    chipseeker_ccre_celltype_agnostic_annotation <- annotatePeakByOverlappingChIPSeekerCCRE(biclustering_grange_i, annotation_celltype_agnostic, categories, featureColname="V6")
    repeatmasker_annotation <- annotatepeakByOverlappingRepeatMasker(biclustering_grange_i, annotationRepeatMasker, repeatMaskerFeatures)
    #####
    chromhmm_short_annotation <- annotatepeakByOverlappingChromHMM(biclustering_grange_i, annotationChromHMM, categoriesChromHMM, featureColname = "V4")
    chromhmm_full_annotation <- annotatepeakByOverlappingChromHMM(biclustering_grange_i, annotationChromHMM, unique(annotationChromHMM$full_anno), featureColname = "full_anno")
    chromhmm_group_annotation <- annotatepeakByOverlappingChromHMM(biclustering_grange_i, annotationChromHMM, unique(annotationChromHMM$group), featureColname = "group")
    chipseeker_ccre_annotations[[glue("{cluster_id}")]] <- chipseeker_ccre_annotation
    repeatmasker_annotations[[glue("{cluster_id}")]] <- repeatmasker_annotation
    chromhmm_short_annotations[[glue("{cluster_id}")]] <- chromhmm_short_annotation
    chromhmm_full_annotations[[glue("{cluster_id}")]] <- chromhmm_full_annotation
    chromhmm_group_annotations[[glue("{cluster_id}")]] <- chromhmm_group_annotation
    chipseeker_ccre_celltype_agnostic_annotations[[glue("{cluster_id}")]] <- chipseeker_ccre_celltype_agnostic_annotation
  }

  # plotting chipseeker + ccre
  anno_ccre <- lapply(chipseeker_ccre_celltype_agnostic_annotations, getAnnoStatCCRE)
  anno_ccre.df <- list_to_dataframe(anno_ccre)
  anno_ccre.df$Feature <- factor(anno_ccre.df$Feature, levels = c(ChIPSeekerCCRECategoriesOrder))
  categoryColumn <- ".id"

  cluster_levels <- cluster_all_levels[cluster_all_levels %in% unique(biclustering_result$label)]
  anno_ccre.df$.id <- factor(anno_ccre.df$.id, levels = rev(cluster_levels))
  print(anno_ccre.df$.id)
  p_ccre_agnostic <- plotAnnoBar.data.frame.one.target(anno_ccre.df,
                                                       categoryColumn = categoryColumn,
                                                       colorOption = 1,
                                                       features=ChIPSeekerCCRECategoriesOrder) + theme(
                                                         legend.position = "right",
                                                         plot.title = element_text(size = 18, ),
                                                         axis.text.x = element_text(size = 12, ),
                                                         axis.text.y = element_text(size = 12, ),
                                                         legend.text = element_text(size = 9),
                                                         legend.key.size = unit(0.5, 'cm'),
                                                         axis.title.x = element_text(size = 20)
                                                       ) +
    ggtitle("Cis-regulatory Elements (cell type agnostic)") +
    guides(fill = guide_legend(title = NULL, ncol = 1, reverse = TRUE))


  # plotting chromhmm
  anno_chromhmm_short <- lapply(chromhmm_short_annotations, getAnnoStatCCRE)
  anno_chromhmm_short.df <- list_to_dataframe(anno_chromhmm_short)
  anno_chromhmm_short.df$Feature <- factor(anno_chromhmm_short.df$Feature, levels = c(unique(annotationChromHMM$V4), "other"))
  categoryColumn <- ".id"
  anno_chromhmm_short.df$.id <- factor(anno_chromhmm_short.df$.id, levels = rev(cluster_levels))
  print(anno_chromhmm_short.df$.id)
  p_chromhmm_short <- plotAnnoBar.data.frame.one.target(anno_chromhmm_short.df,
                                                        categoryColumn = categoryColumn,
                                                        colorOption = 6,
                                                        features=c(unique(annotationChromHMM$V4), "other")) + theme(
                                                          legend.position = "right",
                                                          plot.title = element_text(size = 18, ),
                                                          axis.text.x = element_text(size = 12, ),
                                                          axis.text.y = element_text(size = 12, ),
                                                          legend.text = element_text(size = 9),
                                                          legend.key.size = unit(0.5, 'cm'),
                                                          axis.title.x = element_text(size = 20)
                                                        ) +
    ggtitle("ChromHMM Elements") +
    guides(fill = guide_legend(title = NULL, ncol = 1, reverse = TRUE))

  # save
  if (!dir.exists(output_dir_path)) {
    dir.create(output_dir_path, recursive = TRUE, showWarnings = FALSE)
  }
  ggsave(file.path(output_dir_path, "p_ccre_agnostic.pdf"), plot = p_ccre_agnostic, width = 6, height = 5)
  ggsave(file.path(output_dir_path, "p_chromhmm_short.pdf"), plot = p_chromhmm_short, width = 6, height = 5)

}
