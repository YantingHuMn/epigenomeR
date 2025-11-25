.pkgEnv <- new.env(parent = emptyenv())

#' @export
getRepeatMaskerAnnotation <- function(cache_dir = tools::R_user_dir("epigenomeR", "cache")) {
  if (is.null(.pkgEnv$annotationRepeatMasker)) {
    cache_file <- file.path(cache_dir, "annotationRepeatMasker.rds")

    if (file.exists(cache_file)) {
      message("Loading cached RepeatMasker annotation...")
      .pkgEnv$annotationRepeatMasker <- readRDS(cache_file)
    } else {
      message("Downloading RepeatMasker annotation (first time only, ~185MB)...")
      message("This may take 5-10 minutes...")

      url <- "https://www.repeatmasker.org/genomes/hg38/rmsk4.0.5_rb20140131/hg38.fa.out.gz"
      temp_file <- tempfile(fileext = ".out.gz")

      curl::curl_download(url, temp_file, quiet = FALSE)

      lines <- readLines(gzfile(temp_file))
      lines <- stringr::str_trim(lines[4:length(lines)], "left")
      df <- data.frame(stringr::str_split(lines, "\\s+", simplify = TRUE))
      df$X11 <- gsub("\\/.*", "", df$X11)
      df <- df[!grepl("?", df$X11, fixed = TRUE), ]

      annotationRepeatMasker <- GenomicRanges::makeGRangesFromDataFrame(
        df, keep.extra.columns = TRUE,
        seqnames.field = "X5", start.field = "X6", end.field = "X7"
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
