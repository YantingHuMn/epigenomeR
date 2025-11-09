.load_annotation_data <- function() {
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene

  annotationFile <- system.file("extdata", "ENCFF414OGC_ENCFF806YEZ_ENCFF849TDM_ENCFF736UDR.7group.bed", package = "epigenomeR")
  annotation_celltype_agnostic_file <- system.file("extdata", "GRCh38-cCREs.bed", package = "epigenomeR")
  annotationFileChromHMM <- system.file("extdata", "hg38_genome_100_segments.bed", package = "epigenomeR")
  annotationFileRepeatMasker <- system.file("extdata", "hg38.fa.out", package = "epigenomeR")

  annotationDf <- read.table(annotationFile, header = FALSE, sep="\t", stringsAsFactors=FALSE, quote="")
  annotation <- GenomicRanges::makeGRangesFromDataFrame(annotationDf,
                                                        keep.extra.columns=TRUE,
                                                        seqnames.field="V1",
                                                        start.field="V2",
                                                        end.field="V3",
                                                        strand.field="V6")
  categories <- unique(annotation$V10)

  annotation_celltype_agnostic_df <- read.table(annotation_celltype_agnostic_file, header = FALSE, sep="\t", stringsAsFactors=FALSE, quote="")
  annotation_celltype_agnostic <- GenomicRanges::makeGRangesFromDataFrame(annotation_celltype_agnostic_df,
                                                                          keep.extra.columns=TRUE,
                                                                          seqnames.field="V1",
                                                                          start.field="V2",
                                                                          end.field="V3")

  annotationDfChromHMM <- read.table(annotationFileChromHMM, header = FALSE, sep="\t", stringsAsFactors=FALSE, quote="")
  annotationChromHMM <- GenomicRanges::makeGRangesFromDataFrame(annotationDfChromHMM,
                                                                keep.extra.columns=TRUE,
                                                                seqnames.field="V1",
                                                                start.field="V2",
                                                                end.field="V3",
                                                                strand.field="V6")
  cleaned_strings <- gsub("\\d+_", "", annotationChromHMM$V4)

  replace_strings <- function(strings, replacement_df) {
    sapply(strings, function(x) {
      if (x %in% replacement_df$original) {
        replacement_df$replacement[match(x, replacement_df$original)]
      } else {
        x
      }
    })
  }

  path <- system.file("extdata", "state_annotations_processed.csv", package = "epigenomeR")
  ChromHMM_state <- read.csv(path)

  replacement_df <- ChromHMM_state[, c("mneumonics", "Group")]
  colnames(replacement_df) <- c("original", "replacement")
  updated_strings <- replace_strings(cleaned_strings, replacement_df)
  updated_strings <- unname(updated_strings)
  annotationChromHMM$group <- updated_strings

  annotationChromHMM$full_anno <- annotationChromHMM$V4
  annotationChromHMM$full_anno <- gsub("^\\d+", "", annotationChromHMM$full_anno)
  annotationChromHMM$full_anno <- gsub("_", "", annotationChromHMM$full_anno)
  annotationChromHMM$V4 <- gsub("^\\d+|\\d+$", "", annotationChromHMM$V4)
  annotationChromHMM$V4 <- gsub("_", "", annotationChromHMM$V4)
  categoriesChromHMM <- unique(annotationChromHMM$V4)

  path_rm <- system.file("extdata", "repeatmasker.rds", package = "epigenomeR")
  annotationRepeatMasker <- readRDS(path_rm)

  newV10 <- sub(",.*", "", annotation$V10)
  annotation$V10 <- newV10
  annotation_celltype_agnostic$V6 <- sub(",.*", "", annotation_celltype_agnostic$V6)
  categories <- unique(annotation$V10)
  categoriesOrder <- c("dELS", "pELS", "PLS", "CTCF-only", "DNase-H3K4me3", "DNase-only", "Low-DNase")
  CCREFeatures <- c(categoriesOrder, "other")
  CCREChIPSeekerCategoriesOrder <- c(categoriesOrder, c("5' UTR", "3' UTR", "Exon", "Intron", "other"))
  ChIPSeekerCCRECategoriesOrder <- c("dELS", "pELS", "PLS", "5' UTR", "Exon", "Intron", "3' UTR", "DNase-H3K4me3", "Low-DNase", "DNase-only", "CTCF-only", "other")
  CCREChIPSeekerFeatures <- CCREChIPSeekerCategoriesOrder
  repeatMaskerFeatures <- unique(annotationRepeatMasker$X11)
  repeatMaskerCategoriesOrder <- c(repeatMaskerFeatures, "other")

  allFeatures <- c("Distal Intergenic", "Promoter", "5' UTR", "1st Exon", "1st Intron", "Other Exon", "Other Intron", "3' UTR", "Downstream (<=300)")

  list(
    txdb = txdb,
    annotation = annotation,
    annotation_celltype_agnostic = annotation_celltype_agnostic,
    annotationChromHMM = annotationChromHMM,
    annotationRepeatMasker = annotationRepeatMasker,
    categories = categories,
    categoriesChromHMM = categoriesChromHMM,
    repeatMaskerFeatures = repeatMaskerFeatures,
    ChIPSeekerCCRECategoriesOrder = ChIPSeekerCCRECategoriesOrder,
    allFeatures = allFeatures
  )
}

.annotation_env <- new.env()

.get_annotations <- function() {
  if (is.null(.annotation_env$data)) {
    .annotation_env$data <- .load_annotation_data()
  }
  return(.annotation_env$data)
}

setClass("csCCREAnno",
         representation=representation(
           annoStat="data.frame",
           peakNum="numeric",
           anno="GRanges"
         ))

getAnnoStatCCRE <- function(x) {
  if (!is(x, "csCCREAnno"))
    stop("not supported...")
  return(x@annoStat)
}

getAnnotateStatPeakByOverlappingClosestFeatureHelper <- function(peak, annotation, categories, featureColname = "V10") {
  allPeakCtsList <- list()
  for (category in categories) {
    allPeakCtsList[[category]] <- 0
  }
  overlapRes <- GenomicRanges::findOverlaps(peak, annotation)
  peakHitsUnique <- unique(S4Vectors::queryHits(overlapRes))
  peakHits <- S4Vectors::queryHits(overlapRes)
  peakHitsGR <- peak[peakHits]
  annoHits <- S4Vectors::subjectHits(overlapRes)
  annoHitsGR <- annotation[annoHits]
  peakCenter <- ceiling((GenomicRanges::end(peakHitsGR)-GenomicRanges::start(peakHitsGR))/2 + GenomicRanges::start(peakHitsGR))
  annoStarts <- GenomicRanges::start(annoHitsGR)
  annoEnds <- GenomicRanges::end(annoHitsGR)
  startsDistance <- abs(peakCenter - annoStarts)
  endsDistance <- abs(peakCenter - annoEnds)
  minDistance <- pmin(startsDistance, endsDistance)

  annoHitsDT <- data.table::as.data.table(annoHitsGR)
  annoHitsDT$peakHits <- peakHits
  annoHitsDT$minDistance <- minDistance

  annoHitsDTGroup <- annoHitsDT[annoHitsDT[ , .I[which.min(minDistance)], by = peakHits]$V1]

  peakCategories <- annoHitsDTGroup[[featureColname]]

  peakCategoriesTable <- table(peakCategories)
  return(peakCategoriesTable)
}

annotatePeakByOverlappingClosestFeature <- function(peak, annotation, categories, featureColname = "V10") {
  peakCategoriesTable <- getAnnotateStatPeakByOverlappingClosestFeatureHelper(peak, annotation, categories, featureColname)
  otherLength <- length(peak)-sum(peakCategoriesTable)
  peakFreq <- c(unname(peakCategoriesTable), otherLength) / length(peak) * 100
  res <- data.frame(Feature=c(names(peakCategoriesTable), "other"), Frequency = peakFreq)
  x <- new("csCCREAnno", annoStat = res, peakNum=length(peak))
  return(x)
}

annotatePeakByOverlappingClosestFeatureHelper <- function(peak, annotation, categories, featureColname = "V10") {
  allPeakCtsList <- list()
  for (category in categories) {
    allPeakCtsList[[category]] <- 0
  }
  overlapRes <- GenomicRanges::findOverlaps(peak, annotation)
  peakHitsUnique <- unique(S4Vectors::queryHits(overlapRes))
  peakHits <- S4Vectors::queryHits(overlapRes)
  peakHitsGR <- peak[peakHits]
  annoHits <- S4Vectors::subjectHits(overlapRes)
  annoHitsGR <- annotation[annoHits]
  peakCenter <- ceiling((GenomicRanges::end(peakHitsGR)-GenomicRanges::start(peakHitsGR))/2 + GenomicRanges::start(peakHitsGR))
  annoStarts <- GenomicRanges::start(annoHitsGR)
  annoEnds <- GenomicRanges::end(annoHitsGR)
  startsDistance <- abs(peakCenter - annoStarts)
  endsDistance <- abs(peakCenter - annoEnds)
  minDistance <- pmin(startsDistance, endsDistance)

  annoHitsDT <- data.table::as.data.table(annoHitsGR)
  annoHitsDT$peakHits <- peakHits
  annoHitsDT$minDistance <- minDistance
  annoHitsDT$queryStart <- GenomicRanges::start(peakHitsGR)
  annoHitsDT$queryEnd <- GenomicRanges::end(peakHitsGR)
  annoHitsDT$queryIndex <- peakHitsGR$index

  annoHitsDTGroup <- annoHitsDT[annoHitsDT[ , .I[which.min(minDistance)], by = peakHits]$V1]
  peakCategories <- annoHitsDTGroup[[featureColname]]

  peakCategoriesTable <- table(peakCategories)
  return(annoHitsDTGroup)
}

annotatePeakByOverlappingChIPSeekerCCRE <- function(peak, annotation, categories, featureColname="V10") {
  anno_data <- .get_annotations()
  txdb <- anno_data$txdb
  allFeatures <- anno_data$allFeatures

  chipSeekerAnno <- ChIPseeker::annotatePeak(peak, TxDb=txdb, tssRegion=c(-1000, 1000), verbose=FALSE)
  chipSeekerAnnoFreq <-  c()
  chipSeekerAnnoFreq[allFeatures] <- 0
  chipSeekerAnnoFreq[as.character(chipSeekerAnno@annoStat$Feature)] <- chipSeekerAnno@annoStat$Frequency
  condensedFeatures <- c("5' UTR", "3' UTR", "Exon", "Intron")
  chipSeekerAnnoFreqCondensed <- c()
  chipSeekerAnnoFreqCondensed["5' UTR"] <- chipSeekerAnnoFreq["5' UTR"]
  chipSeekerAnnoFreqCondensed["3' UTR"] <- chipSeekerAnnoFreq["3' UTR"]
  chipSeekerAnnoFreqCondensed["Exon"] <- chipSeekerAnnoFreq["1st Exon"] + chipSeekerAnnoFreq["Other Exon"]
  chipSeekerAnnoFreqCondensed["Intron"] <- chipSeekerAnnoFreq["1st Intron"] + chipSeekerAnnoFreq["Other Intron"]
  chipSeekerAnnoFreqCondensed <- chipSeekerAnnoFreqCondensed * chipSeekerAnno@peakNum / 100

  chipSeekerAnno@anno$annotation <- gsub("\\s*\\([^\\)]+\\)","",chipSeekerAnno@anno$annotation)
  finalAnno <- rep("other", length(peak))
  finalAnno[chipSeekerAnno@anno$annotation %in% condensedFeatures] <- chipSeekerAnno@anno$annotation[chipSeekerAnno@anno$annotation %in% condensedFeatures]
  chipSeekerAnno@anno$index <- 1: length(peak)
  disposeFeatures <- c("Promoter", "Downstream (<=300)", "Distal Intergenic")
  unmappedPeakGR <- chipSeekerAnno@anno[chipSeekerAnno@anno$annotation %in% disposeFeatures]
  if (length(unmappedPeakGR) > 0) {
    annoResult <- annotatePeakByOverlappingClosestFeatureHelper(unmappedPeakGR, annotation, categories, featureColname)
    peakCategories <- annoResult[[featureColname]]
    peakCategoriesTable <- table(peakCategories)
    finalAnno[annoResult$queryIndex] <- annoResult[[featureColname]]
  }

  if (length(unmappedPeakGR)==length(peak)) {
    otherCt <- length(peak) - sum(peakCategoriesTable)
    peakFeature <- c(names(peakCategoriesTable), "other")
    peakFreq <- c(unname(peakCategoriesTable), otherCt) / length(peak) * 100
  } else if (length(unmappedPeakGR)==0) {
    otherCt <- length(peak) - sum(chipSeekerAnnoFreqCondensed)
    peakFeature <- c(names(chipSeekerAnnoFreqCondensed), "other")
    peakFreq <- c(unname(chipSeekerAnnoFreqCondensed), otherCt) / length(peak) * 100
  } else {
    otherCt <- length(peak) - sum(peakCategoriesTable) - sum(chipSeekerAnnoFreqCondensed)
    peakFeature <- c(names(peakCategoriesTable), names(chipSeekerAnnoFreqCondensed), "other")
    peakFreq <- c(unname(peakCategoriesTable), unname(chipSeekerAnnoFreqCondensed), otherCt) / length(peak) * 100
  }
  peak$annotation <- finalAnno
  res <- data.frame(Feature=peakFeature, Frequency = peakFreq)
  x <- new("csCCREAnno", annoStat = res, peakNum=length(peak), anno=peak)
  return(x)
}

annotatepeakByOverlappingRepeatMasker <- function(peak, annotationRepeatMasker, repeatMaskerFeatures) {
  return(annotatePeakByOverlappingClosestFeature(peak, annotationRepeatMasker, repeatMaskerFeatures, featureColname = "X11"))
}

annotatepeakByOverlappingChromHMM <- function(peak, annotationChromHMM, categoriesChromHMM, featureColname = "V4") {
  return(annotatePeakByOverlappingClosestFeature(peak, annotationChromHMM, categoriesChromHMM, featureColname = featureColname))
}
