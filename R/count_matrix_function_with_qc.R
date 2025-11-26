# Count Matrix Function With QC
# Post: build count matrix from bam file with pre-specified genomic regions.
# Parameter: bam_path: A vector of bam file path.
#            regions: Regions should be either an integer or a (vector of) file path ending with .tsv, .txt, .csv, or .bed.
#            libnorm_type: Default to "libnorm" for 1E6 normalization.
#            transformation: Default to "remove0", "libnorm", "log2p1", "qnorm" in order.
#            save_dir: Folder path for saving output files.
#            datasetName_full: the file name, default to "Count_matrix" + region type.
#            save_each_step: Whether save each step, default to TRUE.
#            do_qc: Whether to perform quality control filtering on BAM files.
#            qc_filtered_percentile: Percentile threshold for QC filtering.
# Output: None (saves count matrix and transformed data to files).
count_matrix_function_with_qc <- function(bam_path, regions, save_dir, ref = "hg38", libnorm_type = "libnorm", apply_transformation = TRUE, transformations = NULL, save_each_step = TRUE, datasetName_full = NULL, do_qc = FALSE, qc_filtered_percentile = 0.25) {
    start_time <- Sys.time()

    # Create folder
    if (!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE)
    }

    # initiate packages
    suppressPackageStartupMessages({
        library(R.utils)
        library(GenomicAlignments)
        library(GenomicRanges)
        library(Biostrings)
        library(BSgenome.Hsapiens.NCBI.GRCh38)
        library(BSgenome.Mmusculus.UCSC.mm10)
        library(plyranges)
        library(arrow)
        library(preprocessCore)
        library(tibble)
        library(matrixStats)
        library(BiocParallel)
        library(data.table)
        library(rtracklayer)
    })

    num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))
    register(MulticoreParam(workers = num_cores))

    # Define Regions
    if (is.numeric(regions)) {
        BINSIZE <- regions
        use_custom_region <- FALSE
    } else if (is.character(regions) && all(file.exists(regions))) {
        region_path <- regions
        ext <- tools::file_ext(region_path)
        if (all(ext == "tsv") || all(ext == "txt")) {   
            region_df <- data.table::fread(region_path)
        } else if (all(ext == "csv")){
            region_df <- read.csv(region_path, row.names = NULL)
        } else if (all(ext == "bed")) {
            region_df <- NULL
        } else {
            stop("Error: Invalid region file format. Only support .csv, .tsv, .txt, .bed")
        }
        use_custom_region <- TRUE
    } else {
        stop("Error: Custom regions must be provided. 'regions' argument is missing or invalid.")
    }

    # Variables set
    pos_colname = "pos"
    # qc
    if (do_qc == TRUE) {
        result <- qc(file_paths = bam_path, filtered_percentile = qc_filtered_percentile, save = FALSE)
        vector_crf <- result$filtered_crf
        bamFiles <- bam_path[tools::file_path_sans_ext(basename(bam_path)) %in% vector_crf]
    } else {
        bamFiles <- bam_path
    }
    
    if (ref == "hg38") {
        refGenome <- BSgenome.Hsapiens.NCBI.GRCh38
        # chrSizes <- seqlengths(refGenome)[1:24]
        # chr_list = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")
    } else if (ref == "mm10") {
        refGenome <- BSgenome.Mmusculus.UCSC.mm10
        # chrSizes <- seqlengths(refGenome)[1:21]
        # chr_list = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", 
        #                 "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", 
        #                 "chr16", "chr17", "chr18", "chr19", "chrX", "chrY")
    }
    chr_list <- GenomeInfoDb::standardChromosomes(refGenome)
    chr_list <- chr_list[!tolower(chr_list) %in% c("mt", "chrm", "m", "mito")]
    chrSizes <- seqlengths(refGenome)[chr_list]


    # Post: Use bplapply for parallel chromosome processing
    #       Generate Count Matrix for each chr
    # Process custom regions (no chromosome loop needed)
    if (use_custom_region) {
        ext <- tools::file_ext(region_path)[1]
        
        # Handle CSV/GTF files
        if (!is.null(region_df)) {
            region_df <- fix_region_colnames(region_df)
            bin <- GRanges(seqnames = region_df$seqnames, ranges = IRanges(start = region_df$start, end = region_df$end))
            binChriDataframe <- as.data.frame(bin)[, c("seqnames", "start", "end")]
            colnames(binChriDataframe)[1] <- "CHR"
            
            if (!all(grepl("^chr", seqlevels(bin)))) {
                binChriDataframe$CHR <- paste0("chr", binChriDataframe$CHR)
            }
            if ("gene_id" %in% colnames(region_df)) {
                binChriDataframe$gene_id <- region_df$gene_id
            }
        } 
        # Handle BED files
        else {
            gr_list <- lapply(region_path, function(p) {
                bed_data <- read.table(p, sep = "\t", stringsAsFactors = FALSE)
                GRanges(seqnames = bed_data$V1, ranges = IRanges(start = bed_data$V2 + 1, end = bed_data$V3))
            })
            bin <- reduce(do.call("c", gr_list))
            binChriDataframe <- as.data.frame(bin)[, c("seqnames", "start", "end")]
            colnames(binChriDataframe)[1] <- "CHR"
            
            if (!all(grepl("^chr", seqlevels(bin)))) {
                binChriDataframe$CHR <- paste0("chr", binChriDataframe$CHR)
            }
        }
        
        # Process BAM files for custom regions
        for (k in seq_along(bamFiles)) {
            bamFile <- bamFiles[k]
            temp <- readGAlignmentPairs(bamFile)
            locus <- data.frame(first_start = start(temp@first), 
                                first_end = end(temp@first),
                                last_start = start(temp@last), 
                                last_end = end(temp@last))
            strand <- '*'

            if (nrow(locus) == 0) {
                bamContent <- GRanges()
            } else {
                start <- rowMin(as.matrix(locus))
                end <- rowMax(as.matrix(locus))
                mid <- ceiling((start + end) / 2)
                bamContent <- makeGRangesFromDataFrame(data.frame(
                    seqnames = as.vector(seqnames(temp)), strand = strand, 
                    start = mid, end = mid))
            }

            overlapCount <- countOverlaps(bin, bamContent)
            bamName <- tools::file_path_sans_ext(basename(bamFile))
            binChriDataframe[[bamName]] <- overlapCount
        }
        
        # Format output
        if (length(region_path) == 1 && ("gene_id" %in% colnames(binChriDataframe))) {
            tmp_pos <- binChriDataframe[, c("CHR", "start", "end", "gene_id")]
            tmp_pos$pos <- tmp_pos$gene_id
            pos_df <- data.frame(pos = tmp_pos$pos)
            tmp_wgc <- binChriDataframe[, !(names(binChriDataframe) %in% c("CHR", "start", "end", "gene_id"))]
            binChriDataframe_full <- cbind(pos_df, tmp_wgc)
        } else {
            tmp_pos <- binChriDataframe[, c("CHR", "start", "end")]
            tmp_pos$pos <- paste0(tmp_pos$CHR, "_", tmp_pos$start, "_", tmp_pos$end)
            pos_df <- data.frame(pos = tmp_pos$pos)
            tmp_wgc <- binChriDataframe[, !(names(binChriDataframe) %in% c("CHR", "start", "end"))]
            binChriDataframe_full <- cbind(pos_df, tmp_wgc)
        }
        
    } else {
        # Process fixed bins with parallel chromosome processing
        binChriDataframe_list <- bplapply(chr_list, function(chr_i) {
            chrSizei <- chrSizes[chr_i]
            bin <- tileGenome(chrSizei, tilewidth=BINSIZE, cut.last.tile.in.chrom=TRUE)
            binChriDataframe <- as.data.frame(bin)[, c("start", "end")]
            seqlevels(bin) <- paste0("chr", seqlevels(bin))

            chr_df <- data.frame(CHR = paste0("chr", names(chrSizei)), stringsAsFactors = FALSE)
            binChriDataframe <- cbind(chr_df, binChriDataframe)

            for (k in seq_along(bamFiles)) {
                bamFile <- bamFiles[k]
                temp <- readGAlignmentPairs(bamFile)
                locus <- data.frame(first_start = start(temp@first), 
                                    first_end = end(temp@first),
                                    last_start = start(temp@last), 
                                    last_end = end(temp@last))
                strand <- '*'

                if (nrow(locus) == 0) {
                    bamContent <- GRanges()
                } else {
                    start <- rowMin(as.matrix(locus))
                    end <- rowMax(as.matrix(locus))
                    mid <- ceiling((start + end) / 2)
                    bamContent <- makeGRangesFromDataFrame(data.frame(
                        seqnames = as.vector(seqnames(temp)), strand = strand, 
                        start = mid, end = mid))
                }

                overlapCount <- countOverlaps(bin, bamContent)
                bamName <- tools::file_path_sans_ext(basename(bamFile))
                binChriDataframe[[bamName]] <- overlapCount
            }

            tmp_pos <- binChriDataframe[, c("CHR", "start", "end")]
            tmp_pos$pos <- paste0(tmp_pos$CHR, "_", tmp_pos$start, "_", tmp_pos$end)
            pos_df <- data.frame(pos = tmp_pos$pos)
            tmp_wgc <- binChriDataframe[, !(names(binChriDataframe) %in% c("CHR", "start", "end"))]
            binChriDataframe_final <- cbind(pos_df, tmp_wgc)

            return(binChriDataframe_final)
        })

        binChriDataframe_full <- as.data.frame(do.call(rbind, binChriDataframe_list))
    }

    if (is.null(datasetName_full)) {
        if (is.numeric(regions)) {
            datasetName_full <- paste0("Count_Matrix_", BINSIZE)
        } else if (is.character(regions)) {
            ext <- tools::file_ext(regions)
            add <- ext[1]
            datasetName_full <- paste0("Count_Matrix_", add)
        }
    }

  # Report
    print(warnings())
    print(datasetName_full)

    preprocess_time <- Sys.time()
    preprocess_time_taken <- round(preprocess_time - start_time, 2)
    print(c("prepocess time taken: ", preprocess_time_taken))

    datasetName_full_filename = paste0(datasetName_full, "_orig.feather")
    datasetName_full_dir_filename = file.path(save_dir, datasetName_full_filename)
    write_feather(binChriDataframe_full, datasetName_full_dir_filename)

    saving_time_0 <- Sys.time()
    saving_time_taken_0 <- round(saving_time_0 - preprocess_time, 2)
    print(c("saving time original taken: ", saving_time_taken_0))

    col2idx_time <- Sys.time()
    col2idx_time_taken <- round(col2idx_time - saving_time_0, 2)
    print(c("col2idx time taken: ", col2idx_time_taken))

  # Transformation
    if (apply_transformation == TRUE) {
        rownames(binChriDataframe_full) <- NULL
        binChriDataframe_full = column_to_rownames(binChriDataframe_full, var=pos_colname)
        binChriDataframe_full <- binChriDataframe_full

        apply_transformations(df = binChriDataframe_full, libnorm_type1 = libnorm_type, transformations = transformations, save_each_step = save_each_step, save_dir = save_dir, datasetName_full = datasetName_full)
    }

}