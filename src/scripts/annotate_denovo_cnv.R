# de novo CNV annotation pipeline
#
# Annotate de novo CNVs from a gCNV callset.
# Usage: Rscript annotate_denovo_cnv.R CALLSET BINS PED DCRS NPROC OUTPUT
# * CALLSET - Final callset produced by gCNV pipeline
# * BINS    - Genomic bins file used by the gCNV pipeline
# * PED     - Pedigree file
# * DCRS    - Paths to the dCR matrices in a format accepted by gelpers::read_dcr_list()
# * NPROC   - Number of processors to use
# * OUTPUT  - Path to the output

suppressPackageStartupMessages(library(gelpers))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(parallel))

# Min and max expected average chrX dCR values for XX and XY samples
XX_MIN <- 1.9
XX_MAX <- 2.1
XY_MIN <- 0.9
XY_MAX <- 1.1

# Maximum sequence length supported by tabix. A hack to get an entire
# chromosome from a dCR matrix through Rsamtools.
TABIX_MAX_SEQLEN <- 536870912L

DCR_EVIDENCE <- data.table(M = NA_real_,
                           MF = NA_real_,
                           MM = NA_real_,
                           MD = NA_real_,
                           bins_fail = NA_integer_,
                           bins = NA_integer_)

format_sample_ids <- function(x) {
    is_long <- nchar(x) > 47L
    if (!any(is_long)) {
        return(x)
    }
    replace(x, is_long, paste0(substr(x[is_long], 1, 47), "..."))
}

# Compute the M, MF, MM, etc. table for a trio
dcr_evidence <- function(x, dcr_map) {
    region <- gregion(x$chr, x$start, x$end)
    trio_dcr <- tryCatch({
        get_trio_dcr(region,
                     c(x$sample, x$paternal_id, x$maternal_id),
                     c(x$batch, x$paternal_batch, x$maternal_batch),
                     dcr_map,
                     include_bg = TRUE,
                     keep_all_ranges = TRUE,
                     squeeze = TRUE,
                     reduce = TRUE)

    },
    error = function(cnd) NULL)

    is_chrx <- grepl("X", x$chr, fixed = TRUE)
    default <- DCR_EVIDENCE
    if (is.null(trio_dcr) || nrow(trio_dcr) == 0) {
        return(default)
    }

    dcr <- trio_dcr[setdiff(colnames(trio_dcr), c("chr", "start", "end"))]

    cutoff <- if (is_chrx) 1 else 0.5
    mads <- apply(dcr, 1, mad)
    mads_fail <- sum(mads >= cutoff, na.rm = TRUE)
    if (mads_fail == length(mads)) {
        default$bins_fail <- mads_fail
        default$bins <- length(mads)
        return(default)
    }

    means <- apply(dcr[mads < cutoff, 1:3, drop = FALSE], 2, mean, na.rm = TRUE)
    md <- min(c(abs(means[[1]] - means[[2]]), abs(means[[1]] - means[[3]])))
    rg <- data.table(
        M = means[[1]],
        MF = means[[2]],
        MM = means[[3]],
        MD = md,
        bins_fail = mads_fail,
        bins = length(mads)
    )

    rg
}

sex_from_dcr <- function(samples, dcr) {
    if (is.null(dcr) || nrow(dcr) == 0) {
        return(data.table(sample = samples, chrX_CN = NA_real_, predicted_sex = NA_integer_))
    }
    mean_x <- dcr[, lapply(.SD, mean), .SDcols = !c("chr", "start", "end")]
    mean_x <- melt(mean_x,
                   measure.vars = colnames(mean_x),
                   variable.name = "sample",
                   value.name = "chrX_CN",
                   variable.factor = FALSE)
    mean_x <- mean_x[, predicted_sex := fcase(
        chrX_CN >= XY_MIN & chrX_CN <= XY_MAX, 1L,
        chrX_CN >= XX_MIN & chrX_CN <= XX_MAX, 2L,
        default = NA_integer_
    )]

    mean_x
}

predict_sex <- function(x, dcr_map, nproc = 1L) {
    x <- unique(x)
    hg38_region <- gregion("chrX", 1, TABIX_MAX_SEQLEN)
    hg19_region <- gregion("X", 1, TABIX_MAX_SEQLEN)

    dcr_groups <- split(x$sample, x$batch)
    f <- function(samples, batch) {
        dcr_paths <- gethash(dcr_map, batch)
        if (is.null(dcr_paths)) {
            return(sex_from_dcr(samples, NULL))
        }

        dcr <- tryCatch(
            get_samples_dcr(hg38_region, dcr_paths, samples),
            dcr_parse_error = function(cnd) cnd,
            error = function(cnd) cnd
        )

        if (inherits(dcr, "dcr_parse_error")) (
            dcr <- tryCatch(
                get_samples_dcr(hg19_region, dcr_paths, samples),
                error = function(cnd) NULL
            )
        )

        if (is.null(dcr) || inherits(dcr, "error")) {
            return(sex_from_dcr(samples, NULL))
        } else {
            return(sex_from_dcr(samples, as.data.table(dcr)))
        }
    }
    sexes <- mcmapply(f, dcr_groups, names(dcr_groups), SIMPLIFY = FALSE, mc.cores = nproc)

    rbindlist(sexes)
}

denovo_from_ovp <- function(child, father, mother, bins) {
    gr_c <- df_to_gr(child, cnv = TRUE)
    gr_p <- df_to_gr(father, cnv = TRUE)
    gr_m <- df_to_gr(mother, cnv = TRUE)

    child[, cov_p := cnvCoverage(addPrefix(gr_c, mcols(gr_c)$paternal_id),
                                 addPrefix(gr_p, mcols(gr_p)$sample))]
    child[, cov_m := cnvCoverage(addPrefix(gr_c, mcols(gr_c)$maternal_id),
                                 addPrefix(gr_m, mcols(gr_m)$sample))]

    gr_c_bs <- toBinSpace(gr_c, bins)
    gr_p_bs <- toBinSpace(gr_p, bins)
    gr_m_bs <- toBinSpace(gr_m, bins)

    child[, cov_p_bs := cnvCoverage(addPrefix(gr_c_bs, mcols(gr_c_bs)$paternal_id),
                                    addPrefix(gr_p_bs, mcols(gr_p_bs)$sample))]
    child[, cov_m_bs := cnvCoverage(addPrefix(gr_c_bs, mcols(gr_c_bs)$maternal_id),
                                    addPrefix(gr_m_bs, mcols(gr_m_bs)$sample))]

    dn <- child[cov_m_bs < 0.3 & cov_m < 0.3 & cov_p_bs < 0.3 & cov_p < 0.3]
    dn[, inheritance := "denovo"]

    dn
}

add_parent_batch <- function(child, all_calls) {
    batches <- unique(all_calls[, list(sample, batch)])
    colnames(batches) <- c("paternal_id", "paternal_batch")
    child <- batches[child, nomatch = NULL, on = "paternal_id",]
    colnames(batches) <- c("maternal_id", "maternal_batch")
    child <- batches[child, nomatch = NULL, on = "maternal_id",]

    child
}

recal_cnv_freq <- function(calls, ped) {
    parent_ids <- unique(c(ped$paternal_id, ped$maternal_id))
    parent_cnvs <- calls[sample %in% parent_ids, ]
    parent_cnvs <- parent_cnvs[, list(sc = .N), by = "variant_name"]
    parent_cnvs[, sf := sc / length(..parent_ids)]
    setkey(parent_cnvs, variant_name)

    recal <- calls[, .SD, .SDcols = !c("sc", "sf")]
    recal <- parent_cnvs[recal, on = "variant_name"]
    recal[is.na(sc), sc := 0]
    recal[is.na(sf), sf := 0]

    recal
}

autosome_denovo <- function(calls, bins, ped, dcrs, nproc) {
    calls <- calls[!grepl("X|Y", chr)]
    ped <- ped[, list(sample_id, paternal_id, maternal_id)]

    call_samples <- calls[, list(sample)]
    setkey(call_samples)
    call_samples <- unique(call_samples)

    setkey(ped, sample_id)
    ped <- ped[call_samples, on = c(sample_id = "sample"), nomatch = NULL]
    setkey(ped, paternal_id)
    ped <- ped[call_samples, on = c(paternal_id = "sample"), nomatch = NULL]
    setkey(ped, maternal_id)
    ped <- ped[call_samples, on = c(maternal_id = "sample"), nomatch = NULL]
    setkey(ped, NULL)

    # Recalibrate offspring CNV frequency to parents' -------------------------
    calls <- recal_cnv_freq(calls, ped)
    setkey(calls, sample)

    # Filter to high quality calls and merge in pedigree ----------------------
    child_calls <- calls[PASS_SAMPLE & PASS_QS & sf < 0.01]
    child_calls <- child_calls[ped, on = c(sample = "sample_id"), nomatch = NULL]

    # Get preliminary de novo calls based on CNV overlap ----------------------
    log_info("making initial de novo predictions based on CNV overlap")
    dn <- denovo_from_ovp(child_calls,
                          calls[sample %in% ped$paternal_id],
                          calls[sample %in% ped$maternal_id],
                          bins)
    log_info(sprintf("found %d de novo events based on overlap", nrow(dn)))

    # Gather dCR evidence -----------------------------------------------------
    log_info("gathering dCR evidence")
    dn <- add_parent_batch(dn, calls)
    old_dtthreads <- getDTthreads()
    # Ensure that datat.table runs single-threaded inside a mclapply call.
    setDTthreads(1)
    rg <- mclapply(seq_len(nrow(dn)),
                   \(i) dcr_evidence(as.list(dn[i, ]), dcrs),
                   mc.cores = nproc) |>
        rbindlist()
    setDTthreads(old_dtthreads)

    # Regenotype --------------------------------------------------------------
    log_info("regenotyping")
    rg <- cbind(rg, dn[, list(CN, chr, inheritance, svtype)])

    rg[is.na(MF) | is.na(MM), inheritance := "fail_miss_parents"]
    rg[abs(M - CN) > 0.175 & CN <= 3, inheritance := "fail_M"]
    rg[MD < 0.7, inheritance := "fail_MD"]
    rg[abs(M - CN) > 0.175 & CN <= 3 & MD < 0.7, inheritance := "fail_M_MD"]
    rg[bins_fail / bins >= 0.5 & bins < 100, inheritance := "fail_bad_bins"]
    rg[CN > 2 & inheritance == "denovo" & (MF > 2.5 | MM > 2.5), inheritance := "inherited"]
    rg[CN < 2 & inheritance == "denovo" & (MF < 1.5 | MM < 1.5), inheritance := "inherited"]

    rg[, mos_dup_thresh := NA_real_]
    rg[inheritance == "denovo" & svtype == "DUP",
       mos_dup_thresh := sapply(qnorm(0.98, 0, 0.5 / sqrt(bins)) + 2, max, 2.1)]

    rg[, mos_del_thresh := NA_real_]
    rg[inheritance == "denovo" & svtype == "DEL",
       mos_del_thresh := sapply(qnorm(0.98, 0, 0.5 / sqrt(bins)) * -1 + 2, min, 1.9)]

    rg[MF > mos_dup_thresh, inheritance := "mosaic_father"]
    rg[MM > mos_dup_thresh, inheritance := "mosaic_mother"]
    rg[MF < mos_del_thresh, inheritance := "mosaic_father"]
    rg[MM < mos_del_thresh, inheritance := "mosaic_mother"]
    rg[MF > mos_dup_thresh & MM > mos_dup_thresh, inheritance := "mosaic_both"]
    rg[MF < mos_del_thresh & MM < mos_del_thresh, inheritance := "mosaic_both"]

    dn[, inheritance := rg$inheritance]

    dn
}

chrx_denovo <- function(calls, bins, ped, dcrs, nproc) {
    # Need the sex column of every sample
    ped <- ped[sample_id %in% unique(calls$sample) & !is.na(sex)]

    # Predict the sex of all samples ------------------------------------------
    log_info("validating sex of samples")
    predicted_sex <- predict_sex(calls[, list(sample, batch)], dcrs, nproc)
    fail_dcr <- predicted_sex[is.na(chrX_CN)]
    if (nrow(fail_dcr) > 0) {
        log_warn(sprintf("%s samples failed chrX CN calling", nrow(fail_dcr)))
    }

    # Report samples with chrX copy number outside range ----------------------
    fail_sex <- predicted_sex[!is.na(chrX_CN) & is.na(predicted_sex)]
    if (nrow(fail_sex) > 0) {
        sample_ids <- format_sample_ids(fail_sex$sample)
        max_width <- max(nchar(sample_ids))
        header_fmt <- paste0("%", max_width, "s  chrX_CN")
        body_fmt <- paste0("%-", max_width, "s  %7.3f")
        log_warn(c(
            sprintf("%s samples have out-of-range chrX CN", nrow(fail_sex)),
            sprintf(header_fmt, "sample"),
            sprintf(body_fmt, sample_ids, fail_sex$chrX_CN)
        ))
    }

    # Remove callset samples without a pedigree -------------------------------
    sex_match_df <- predicted_sex[!is.na(predicted_sex)]
    sex_match_df <- ped[, list(sample_id, sex)][sex_match_df, on = c(sample_id = "sample")]
    setnames(sex_match_df, "sample_id", "sample")
    missing_sex_cnt <- nrow(sex_match_df[is.na(sex)])
    if (missing_sex_cnt > 0) {
        log_warn(
            paste0("discarding ", missing_sex_cnt, " samples without pedigree")
        )
    }

    # Remove samples whose predicted sex doesn't match pedigree sex -----------
    sex_match_df <- sex_match_df[!is.na(sex)]
    mismatch_sex <- sex_match_df[predicted_sex != sex]
    if (nrow(mismatch_sex) > 0) {
        sample_ids <- format_sample_ids(mismatch_sex$sample)
        max_width <- max(nchar(sample_ids))
        header_fmt <- paste0("%", max_width, "s  chrX_CN  predicted_sex  pedigree_sex")
        body_fmt <- paste0("%-", max_width, "s  %7.3f  %13s  %12s")
        log_warn(c(
            paste0("discarding ", nrow(mismatch_sex), " samples with mismatched sex"),
            sprintf(header_fmt, "sample"),
            sprintf(body_fmt,
                    sample_ids,
                    mismatch_sex$chrX_CN,
                    ifelse(mismatch_sex$predicted_sex == 1, "male", "female"),
                    ifelse(mismatch_sex$sex == 1, "male", "female"))
        ))
    }

    pass_sex <- unique(sex_match_df[predicted_sex == sex], by = "sample")
    log_info(sprintf("continuing analysis with %d samples passing sex check", nrow(pass_sex)))

    # Filter callset and pedigree ---------------------------------------------
    calls <- calls[sample %in% pass_sex$sample]
    ped <- ped[sample_id %in% calls$sample]
    ped <- ped[paternal_id %in% calls$sample]
    ped <- ped[maternal_id %in% calls$sample]

    # Recalibrate offspring CNV frequency to parents' -------------------------
    calls <- recal_cnv_freq(calls, ped)

    # Filter to high quality calls and merge in pedigree ----------------------
    child_calls <- calls[PASS_SAMPLE & PASS_QS & sf < 0.01 & grepl("X", chr)]
    child_calls <- child_calls[ped, on = c(sample = "sample_id"), nomatch = NULL]

    # Get preliminary de novo calls based on CNV overlap ----------------------
    log_info("making initial de novo predictions based on CNV overlap")
    dn <- denovo_from_ovp(child_calls,
                          calls[sample %in% ped$paternal_id & grepl("X", chr)],
                          calls[sample %in% ped$maternal_id & grepl("X", chr)],
                          bins)
    log_info(sprintf("found %d de novo events based on overlap", nrow(dn)))

    # Gather dCR evidence -----------------------------------------------------
    log_info("gathering dCR evidence")
    dn <- add_parent_batch(dn, calls)
    old_dtthreads <- getDTthreads()
    setDTthreads(1)
    rg <- mclapply(seq_len(nrow(dn)),
                   \(i) dcr_evidence(as.list(dn[i, ]), dcrs),
                   mc.cores = nproc) |>
        rbindlist()
    setDTthreads(old_dtthreads)

    # Regenotype --------------------------------------------------------------
    log_info("regenotyping")
    rg <- cbind(rg, dn[, list(sample, chr, CN, inheritance, svtype, paternal_id, maternal_id)])
    rg <- pass_sex[, list(sample, sex)][rg, on = "sample"]
    rg <- pass_sex[, list(paternal_id = sample, paternal_sex = sex)][rg, on = "paternal_id"]
    rg <- pass_sex[, list(maternal_id = sample, maternal_sex = sex)][rg, on = "maternal_id"]

    # Check sex of parents ----------------------------------------------------
    rg[is.na(MF) | is.na(MM), inheritance := "fail_miss_parents"]
    rg[paternal_sex != 1 | maternal_sex != 2, inheritance := "fail_parental_sex"]
    rg[is.na(M), inheritance := "fail_miss_child"]

    # Check offspring is consistent with SV type and sex ----------------------
    # sex  svtype         expected_M
    #  XX     DEL  0.75 <= M <= 1.25
    #  XY     DEL  0.00 <= M <= 0.25
    #  XX     DUP  2.75 <= M
    #  XY     DUP  1.75 <= M
    rg[sex == 2L & svtype == "DEL" & abs(M - 1) > 0.25, inheritance := "fail_M"]
    rg[sex == 1L & svtype == "DEL" & M > 0.25, inheritance := "fail_M"]
    rg[sex == 2L & svtype == "DUP" & M < 2.75, inheritance := "fail_M"]
    rg[sex == 1L & svtype == "DUP" & M < 1.75, inheritance := "fail_M"]

    # Check for CNV evidence in parents ---------------------------------------
    # We check if a parent has a possible CNV in the same region as a child CNV
    # by comparing the parent chrX ploidy to the parent dCR value for the region.
    # For example, if a parent has 1 chrX, then it is expected that the dCR value
    # for any region on chrX will be about 1. If either parent has a suspected CNV
    # that matches the SV type of the CNV in the child, then the CNV in the child
    # is likely inherited. A threshold of max deviation 0.15 was determined by
    # plotting deviations for random sized subsets of some dCR matrices.
    rg[svtype == "DEL" & paternal_sex - MF > 0.15, inheritance := "fail_MF"]
    rg[svtype == "DEL" & maternal_sex - MM > 0.15, inheritance := "fail_MF"]
    rg[svtype == "DEL" & paternal_sex - MF > 0.15 & maternal_sex - MM > 0.15, inheritance := "fail_MF_MM"]
    rg[svtype == "DUP" & MF - paternal_sex > 0.15, inheritance := "fail_MF"]
    rg[svtype == "DUP" & MM - maternal_sex > 0.15, inheritance := "fail_MM"]
    rg[svtype == "DUP" & MF - paternal_sex > 0.15 & MM - maternal_sex > 0.15, inheritance := "fail_MF_MM"]
    rg[bins_fail / bins > 0.5 & bins < 100, inheritance := "fail_bad_bins"]

    dn$inheritance <- rg$inheritance

    dn
}

# Parse command-line arguments ------------------------------------------------
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) != 6) {
    stop("Incorrect number of arguments", call. = FALSE)
}
callset_path <- argv[[1]]
bins_path <- argv[[2]]
ped_path <- argv[[3]]
dcrs_path <- argv[[4]]
nproc <- as.integer(argv[[5]])
output <- argv[[6]]

# Read inputs -----------------------------------------------------------------
log_info("reading callset")
raw_calls <- read_callset(callset_path)
is_hg19 <- any(c(as.character(1:22), "X", "Y") %in% raw_calls$chr)

log_info("reading bins")
bins <- read_gcnv_bins(bins_path, reduce = is_hg19)

log_info("reading pedigree")
ped <- suppressWarnings(read_pedigree(ped_path))

log_info("reading dCR paths")
dcrs <- read_dcr_list(dcrs_path)

log_info("calling autosome de novo CNVs")
dn_auto <- autosome_denovo(raw_calls, bins, ped, dcrs, nproc)
log_info("completed calling autosome de novo CNVs")

log_info("calling chrX de novo CNVs")
dn_chrx <- chrx_denovo(raw_calls, bins, ped, dcrs, nproc)
log_info("completed calling chrX de novo CNVs")

# Write output to file --------------------------------------------------------
common_cols <- intersect(colnames(dn_auto), colnames(dn_chrx))
dn_all <- rbind(dn_auto[, ..common_cols], dn_chrx[, ..common_cols])
dn_all <- ped[, c("family_id", "sample_id")][dn_all,
                                             on = c(sample_id = "sample"),
                                             mult = "first"]
setnames(dn_all, "sample_id", "sample")

log_info("writing output")
write.table(dn_all,
            file = output,
            sep = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)
log_info("done")
