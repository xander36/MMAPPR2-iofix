#' @title Generate candidate mutations and consequences in peak regions
#'
#' @name generateCandidates
#'
#' @usage Follows the \code{\link{peakRefinement}} step and produces a
#' \code{\linkS4class{md}} object ready for
#' \code{\link{outputmd}}.
#'
#' @param md The \code{\linkS4class{md}} object to be analyzed.
#'
#' @return A \code{\linkS4class{md}} object with the \code{candidates}
#'   slot filled with a \code{\link[GenomicRanges]{GRanges}} object for each
#'   peak chromosome containing variants and predicted consequences.
#' @export
#'
#' @examples
#' if (requireNamespace('MMAPPR2data', quietly=TRUE)) {
#'     mmappr_param <- mmapprParam(wtFiles = MMAPPR2data::exampleWTbam(),
#'                                 mutFiles = MMAPPR2data::exampleMutBam(),
#'                                 refFasta = MMAPPR2data::goldenFasta(),
#'                                 gtf = MMAPPR2data::gtf(),
#'                                 outputFolder = tempOutputFolder())
#' }
#'
#' \dontrun{
#' md <- mmapprData(mmappr_param)
#' postCalcDistMD <- calculateDistance(md)
#' postLoessMD <- loessFit(postCalcDistMD)
#' postPrePeakMD <- prePeak(postLoessMD)
#' postPeakRefMD <- peakRefinement(postPrePeakMD)
#'
#' postCandidatesMD <- generateCandidates(postPeakRefMD)
#' }
#'
NULL

# ---- helpers (species-agnostic) ------------------------------------------------

# pick ONE canonical seqlevel style from Seqinfo (avoid vector-of-styles warnings)
.choose_target_style <- function(si) {
  st <- unique(unlist(GenomeInfoDb::seqlevelsStyle(si)))
  if (length(st)) st[1] else NA_character_
}

# Seqinfo from FaFile (no hardcoded genome string)
.faSeqinfo <- function(fa) {
  faidx   <- Rsamtools::scanFaIndex(fa)
  seqlens <- GenomeInfoDb::seqlengths(faidx)
  GenomeInfoDb::Seqinfo(seqnames = names(seqlens), seqlengths = seqlens)
}

# ---- main API ------------------------------------------------------------------

generateCandidates <- function(md) {
  # get GRanges representation of peak
  .messageAndLog("Getting Variants in Peak", outputFolder(param(md)))
  peakGRanges <- lapply(md@peaks, .getPeakRange)

  .messageAndLog("Peak report generation", outputFolder(param(md)))

  #Write peak info to log
  .messageAndLog("Peak Summary:", outputFolder(param(md)))
  
  for (seqname in names(md@peaks)){
    peak <- md@peaks[[seqname]]
    peak_pos <- peak$peakPosition

    density_apex_pos <- NA
    density_apex_val <- NA
    if (!is.null(peak$densityData) && all(c("x","y") %in% names(peak$densityData))) {
      max_index <- which.max(peak$densityData$y)
      density_apex_pos = peak$densityData$x[max_index]
      density_apex_val = peak$densityData$y[max_index]
    }

    log_text <- paste0(
      "seqname: ", seqname, "\n",
      "start-end positions:\n",
      as.integer(peak$start), "-", as.integer(peak$end), "\n",
      "Position of curve peak: ", as.integer(round(peak_pos)), "\n",
      "Position of maximum density: ", as.integer(round(density_apex_pos)), "\n",
      "Value of maximum density: ", density_apex_val, "\n\n"
    )

    .messageAndLog(log_text, outputFolder(param(md)))
  }
  

  # call variants in peaks
  md@candidates$snps <- lapply(peakGRanges, FUN = .getVariantsForRange, param = md@param)
  
  # remove NULL peak results
  allPeaks <- length(md@candidates$snps)
  md@candidates$snps <- md@candidates$snps[!vapply(md@candidates$snps, is.null, logical(1))]
  noNullPeaks <- length(md@candidates$snps)
  
  if (allPeaks != noNullPeaks) {
    peaksRemoved <- allPeaks - noNullPeaks
    oF <- outputFolder(param(md))
    .messageAndLog("Warning: peaks without valid snps were called but then removed", oF)
    .messageAndLog(paste("Number of peaks removed:", peaksRemoved), oF)
  }
  
  # predict effects of variants
  .messageAndLog("Predicting Variant Effects", outputFolder(param(md)))
  md@candidates$effects <- lapply(md@candidates$snps, FUN = .predictEffects, param = md@param)
  
  # add differentially expressed genes
  .messageAndLog("Identifying Differentially Expressed Genes in Peak", outputFolder(param(md)))
  md@candidates$diff <- lapply(peakGRanges, FUN = .addDiff, param = md@param)
  
  # density score and order variants
  .messageAndLog("Ordering Candidates by Position", outputFolder(param(md)))
  md@candidates <- .scoreVariants(md@candidates, md@peaks)
  
  md
}

# ---- internals -----------------------------------------------------------------

.getPeakRange <- function(peakList) {
  ir <- IRanges::IRanges(start = as.numeric(peakList$start),
                         end   = as.numeric(peakList$end),
                         names = peakList$seqname)
  GenomicRanges::GRanges(seqnames = names(ir), ranges = ir)
}

.getVariantsForRange <- function(inputRange, param) {
  mergedBam <- file.path(outputFolder(param), "merged.tmp.bam")

  # choose the BAM (BamFile if single; character path if merged)
  if (length(param@mutFiles) < 2) {
    mutBam <- param@mutFiles[[1]]  # BamFile
  } else {
    mutBam <- Rsamtools::mergeBam(param@mutFiles, destination = mergedBam, region = inputRange) 
  }

  # index whatever we have
  if (inherits(mutBam, "BamFile")) {
    p <- BiocGenerics::path(mutBam)
    if (!file.exists(paste0(p, ".bai"))) Rsamtools::indexBam(p)

  } else if (is.character(mutBam)) {
    if (!file.exists(paste0(mutBam, ".bai"))) Rsamtools::indexBam(mutBam)
  } else {
    stop("mutBam must be a BamFile or a character path; got: ", class(mutBam)[1])
  }

  # always hand downstream a BamFileList
  mutBamList <- if (inherits(mutBam, "BamFile")) {
    Rsamtools::BamFileList(mutBam)
  } else {
    Rsamtools::BamFileList(Rsamtools::BamFile(mutBam))
  }

  # (optional) keep old name if later code expects `mutBam`
  mutBam <- mutBamList
  
  # FASTA seqinfo & a single canonical style
  fa     <- param@refGenome                     # Rsamtools::FaFile
  fa_si  <- .faSeqinfo(fa)
  tstyle <- .choose_target_style(fa_si)
  
  # Make the query window use the FASTA style so pileup sees contigs
  if (!is.na(tstyle)) GenomeInfoDb::seqlevelsStyle(inputRange) <- tstyle
  
  # Build tallies via Rsamtools::pileup (no callbacks, no GMAP)
  pile <- .pileupVariants(
    bams           = mutBam,
    genome         = fa,
    which          = inputRange,
    minAltDepth    = 1L,
    baseOnly       = TRUE,
    minBaseQuality = as.integer(minBaseQuality(param)),
    minMapQuality  = as.integer(minMapQuality(param))
  )
  
  # Anchor VRanges to FASTA seqinfo and style
  common <- intersect(GenomeInfoDb::seqlevels(pile), GenomeInfoDb::seqnames(fa_si))
  pile   <- GenomeInfoDb::keepSeqlevels(pile, common, pruning.mode = "coarse")
  GenomeInfoDb::seqinfo(pile) <- fa_si[GenomeInfoDb::seqlevels(pile)]
  if (!is.na(tstyle)) GenomeInfoDb::seqlevelsStyle(pile) <- tstyle
  
  # Call variants
  resultVr <- VariantTools::callVariants(pile)
  resultVr <- resultVr[ VariantAnnotation::altDepth(resultVr) / VariantAnnotation::totalDepth(resultVr) > 0.8 ]
  
  # cleanup merged BAM + index if created
  if (file.exists(mergedBam)) {
    file.remove(mergedBam)
    bai <- paste0(mergedBam, ".bai")
    if (file.exists(bai)) file.remove(bai)
  }
  
  if (length(resultVr) > 0) {
    Biobase::sampleNames(resultVr) <- paste0(names(param@mutFiles), collapse = " -- ")
    S4Vectors::mcols(resultVr) <- NULL
    resultVr
  } else {
    NULL
  }
}

.predictEffects <- function(inputVariants, param) {
  # Build TxDb with contig lengths from FASTA; species-agnostic
  fa     <- Rsamtools::FaFile(refFasta(param))
  fa_si  <- .faSeqinfo(fa)                     # Seqinfo from FASTA (helper defined above)
  tstyle <- .choose_target_style(fa_si)        # single canonical style (helper)
  
  # Use FASTA Seqinfo *directly* as chrominfo (valid type for txdbmaker/rtracklayer)
  chrominfo <- fa_si
  
  if (requireNamespace("txdbmaker", quietly = TRUE)) {
    txdb <- txdbmaker::makeTxDbFromGFF(
      file      = gtf(param),
      format    = "gtf",
      chrominfo = chrominfo                    # <-- Seqinfo (not S4 DataFrame)
    )
  } else {
    warning("Package 'txdbmaker' not installed; using GenomicFeatures::makeTxDbFromGFF() (deprecated).")
    txdb <- GenomicFeatures::makeTxDbFromGFF(
      file      = gtf(param),
      format    = "gtf",
      chrominfo = chrominfo                    # <-- Seqinfo works here too
    )
  }
  
  # Keep only contigs present in FASTA and use the FASTA's single style
  keep <- intersect(GenomeInfoDb::seqlevels(txdb), GenomeInfoDb::seqnames(fa_si))
  txdb <- GenomeInfoDb::keepSeqlevels(txdb, keep, pruning.mode = "coarse")
  if (!is.na(tstyle)) GenomeInfoDb::seqlevelsStyle(txdb) <- tstyle
  
  inputVariants <- GenomeInfoDb::keepSeqlevels(inputVariants, keep, pruning.mode = "coarse")
  if (!is.na(tstyle)) GenomeInfoDb::seqlevelsStyle(inputVariants) <- tstyle
  
  VariantAnnotation::predictCoding(
    query     = inputVariants,
    subject   = txdb,
    seqSource = fa,
    varAllele = Biostrings::DNAStringSet(VariantAnnotation::alt(inputVariants))
  )
}

.addDiff <- function(peakGRange, param) {
  # prep gene models from GTF (genes only)
  suppressMessages(
    genes <- data.table::fread(gsub(".bgz", "", gtf(param)), header = FALSE)
  )
  genes <- genes[V3 == "gene"
  ][, gene_id   := gsub(".*gene_id \"(.*?)\";.*",  "\\1", V9)
  ][, gene_name := gsub(".*gene_name \"(.*?)\";.*","\\1", V9)
  ][, .(seqnames = V1, start = V4, end = V5, strand = V7, gene_id, gene_name)]
  genes <- GenomicRanges::makeGRangesFromDataFrame(genes, keep.extra.columns = TRUE)
  
  # Harmonize seqlevels/style with FASTA (through the peak's style)
  fa     <- param@refGenome
  fa_si  <- .faSeqinfo(fa)
  tstyle <- .choose_target_style(fa_si)
  if (!is.na(tstyle)) {
    GenomeInfoDb::seqlevelsStyle(genes)      <- tstyle
    GenomeInfoDb::seqlevelsStyle(peakGRange) <- tstyle
  }
  common <- intersect(GenomeInfoDb::seqlevels(genes), GenomeInfoDb::seqlevels(peakGRange))
  genes  <- GenomeInfoDb::keepSeqlevels(genes, common, pruning.mode = "coarse")
  genes  <- IRanges::subsetByOverlaps(x = genes, ranges = peakGRange)
  
  # summarize counts and compute simple log2FC
  readfiles <- c(param@wtFiles, param@mutFiles)
  counts <- GenomicAlignments::summarizeOverlaps(
    features = genes,
    reads    = readfiles,
    param    = Rsamtools::ScanBamParam(which = peakGRange)
  )
  num_wt <- length(param@wtFiles)
  countDF <- data.table::as.data.table(SummarizedExperiment::assays(counts)$counts)
  countDF[, ave_wt := rowMeans(countDF[, 1:(num_wt + 1)])
  ][,   ave_mt := rowMeans(countDF[, (num_wt + 1):length(countDF)])
  ][,   log2FC := round(log2(ave_mt / ave_wt), 3)]
  
  mcols(genes) <- cbind(mcols(genes), countDF)
  
  # filter for "differentially expressed" using simple thresholds
  genes[(abs(genes$log2FC) > 1 | is.na(genes$log2FC)) &
          (genes$ave_wt > 10 | genes$ave_mt > 10)]
}

.scoreVariants <- function(candList, peaks) {
  for (GRname in names(candList)) {
    GR <- candList[[GRname]]
    for (seqname in names(GR)) {
      # density calculation
      densityFunc <- peaks[[seqname]]$densityFunction
      stopifnot(!is.null(densityFunc))
      positions <- BiocGenerics::start(GR[[seqname]]) +
        ((BiocGenerics::width(GR[[seqname]]) - 1) / 2)
      densityCol <- vapply(positions, densityFunc, FUN.VALUE = numeric(1))
      S4Vectors::mcols(candList[[GRname]][[seqname]])$peakDensity <- densityCol
      
      # re-order
      candList[[GRname]][[seqname]] <- .orderVariants(candList[[GRname]][[seqname]])
    }
  }
  candList
}

.orderVariants <- function(candidateGRanges) {
  if (!(class(candidateGRanges) %in% c("VRanges", "GRanges"))) {
    .messageAndLog("invalid data type for sorting. (orderVariants)", outputFolder(param(md)))
    return(candidateGRanges)
  }
  if (!is.null(candidateGRanges$CONSEQUENCE)) {
    impactLevels <- c("synonymous", "nonsynonymous", "frameshift", "nonsense")
    orderVec <- order(match(candidateGRanges$CONSEQUENCE, impactLevels),
                      candidateGRanges$peakDensity, decreasing = TRUE)
  } else {
    orderVec <- order(candidateGRanges$peakDensity, decreasing = TRUE)
  }
  candidateGRanges[orderVec]
}

# GMAP-free pileup → VRanges (substitutions only)
.pileupVariants <- function(
    bams,
    genome,                      # Rsamtools::FaFile (indexed)
    which = NULL,                # GRanges (e.g., inputRange)
    minAltDepth   = 1L,
    baseOnly      = TRUE,
    minBaseQuality = 0L,
    minMapQuality  = 0L
) {
  # normalize BAMs and sample names
  if (is.character(bams)) {
    bams <- Rsamtools::BamFileList(bams)
  }
  if (is.null(names(bams))) {
    names(bams) <- basename(Rsamtools::path(bams))
  }
  
  # pileup config (no callbacks)
  pupar <- Rsamtools::PileupParam(
    distinguish_nucleotides = TRUE,
    distinguish_strands     = FALSE,
    min_base_quality        = as.integer(minBaseQuality),
    min_mapq                = as.integer(minMapQuality),
    include_insertions      = FALSE,
    include_deletions       = FALSE
  )
  sbpar <- if (!is.null(which)) Rsamtools::ScanBamParam(which = which) else Rsamtools::ScanBamParam()
  
  vr_list <- vector("list", length(bams))
  for (i in seq_along(bams)) {
    smp <- names(bams)[i]
    df  <- Rsamtools::pileup(bams[[i]], scanBamParam = sbpar, pileupParam = pupar)
    
    if (NROW(df) == 0L) {
      vr_list[[i]] <- VariantAnnotation::VRanges()
      next
    }
    
    # keep standard bases only
    df <- df[df$nucleotide %in% c("A","C","G","T","N"), c("seqnames","pos","nucleotide","count")]
    if (NROW(df) == 0L) {
      vr_list[[i]] <- VariantAnnotation::VRanges()
      next
    }
    
    # reference base per unique position from FASTA
    pos_df <- unique(df[, c("seqnames","pos")])
    gr     <- GenomicRanges::GRanges(seqnames = pos_df$seqnames,
                                     ranges   = IRanges::IRanges(pos_df$pos, width = 1))
    pos_df$ref <- as.character(Rsamtools::getSeq(genome, gr))
    
    # total depth per pos (A/C/G/T only)
    tot <- stats::aggregate(count ~ seqnames + pos,
                            df[df$nucleotide %in% c("A","C","G","T"), ],
                            sum)
    names(tot)[3] <- "totalDepth"
    
    # ref counts per pos
    ref_counts <- merge(
      pos_df[, c("seqnames","pos","ref")],
      df,
      by.x = c("seqnames","pos","ref"),
      by.y = c("seqnames","pos","nucleotide"),
      all.x = TRUE
    )
    ref_counts$count[is.na(ref_counts$count)] <- 0L
    names(ref_counts)[names(ref_counts) == "count"] <- "refDepth"
    
    # alt rows: nucleotide != ref
    alt_df <- merge(df, pos_df[, c("seqnames","pos","ref")], by = c("seqnames","pos"))
    alt_df <- alt_df[alt_df$nucleotide != alt_df$ref, ]
    if (baseOnly) {
      alt_df <- alt_df[alt_df$ref != "N" & alt_df$nucleotide != "N", ]
    }
    if (NROW(alt_df) == 0L) {
      vr_list[[i]] <- VariantAnnotation::VRanges()
      next
    }
    
    # join totalDepth & refDepth
    alt_df <- merge(alt_df, tot, by = c("seqnames","pos"), all.x = TRUE)
    alt_df$totalDepth[is.na(alt_df$totalDepth)] <- 0L
    alt_df <- merge(alt_df, ref_counts[, c("seqnames","pos","refDepth")],
                    by = c("seqnames","pos"), all.x = TRUE)
    alt_df$refDepth[is.na(alt_df$refDepth)] <- 0L
    
    # build VRanges rows (one per alt allele) with altDepth filter
    keep <- alt_df$count >= minAltDepth
    alt_df <- alt_df[keep, ]
    if (NROW(alt_df) == 0L) {
      vr_list[[i]] <- VariantAnnotation::VRanges()
      next
    }
    
    vr <- VariantAnnotation::VRanges(
      seqnames    = S4Vectors::Rle(alt_df$seqnames),
      ranges      = IRanges::IRanges(start = alt_df$pos, width = 1L),
      ref         = alt_df$ref,
      alt         = alt_df$nucleotide,
      refDepth    = S4Vectors::Rle(as.integer(alt_df$refDepth)),
      altDepth    = S4Vectors::Rle(as.integer(alt_df$count)),
      totalDepth  = S4Vectors::Rle(as.integer(alt_df$totalDepth)),
      sampleNames = S4Vectors::Rle(smp, nrow(alt_df))
    )
    
    vr_list[[i]] <- vr
  }
  
  if (length(vr_list) == 0L) VariantAnnotation::VRanges() else do.call(c, vr_list)
}
