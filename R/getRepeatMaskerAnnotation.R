.pkgEnv <- new.env(parent = emptyenv())

getRepeatMaskerAnnotation <- function(cache_dir = tools::R_user_dir("epigenomeR", "cache")) {

  if (is.null(.pkgEnv$annotationRepeatMasker)) {

    cache_file <- file.path(cache_dir, "annotationRepeatMasker.rds")

    if (file.exists(cache_file)) {
      message("Loading cached RepeatMasker annotation...")
      .pkgEnv$annotationRepeatMasker <- readRDS(cache_file)
    } else {
      message("Downloading RepeatMasker annotation (first time only, ~185MB)...")

      # downlaod
      url <- "https://www.repeatmasker.org/genomes/hg38/rmsk4.0.5_rb20140131/hg38.fa.out.gz"
      temp_file <- tempfile(fileext = ".out.gz")
      download.file(url, temp_file, mode = "wb")

      # 用你原来的处理方法
      lines <- readLines(gzfile(temp_file))
      lines <- stringr::str_trim(lines[4:length(lines)], "left")
      annotationRepeatMasker <- data.frame(stringr::str_split(lines, "\\s+", simplify = TRUE))
      annotationRepeatMasker$X11 <- gsub("\\/.*", "", annotationRepeatMasker$X11)
      annotationRepeatMasker <- annotationRepeatMasker[!grepl("?", annotationRepeatMasker$X11, fixed = TRUE), ]
      annotationRepeatMasker <- GenomicRanges::makeGRangesFromDataFrame(
        annotationRepeatMasker,
        keep.extra.columns = TRUE,
        seqnames.field = "X5",
        start.field = "X6",
        end.field = "X7"
      )

      dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(annotationRepeatMasker, cache_file)
      unlink(temp_file)

      .pkgEnv$annotationRepeatMasker <- annotationRepeatMasker
      message("Done! Cached for future use.")
    }
  }
  return(.pkgEnv$annotationRepeatMasker)
}
