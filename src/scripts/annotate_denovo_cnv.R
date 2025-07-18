# de novo CNV annotation pipeline

DEFAULT_OPTS <- list(recal_freq = TRUE,
                     hq_cols = c("PASS_SAMPLE", "PASS_QS"),
                     max_freq = 0.01,
                     cpus = 1L,
                     allosomes = TRUE)

DEFAULT_ARGS <- list(CALLSET = NULL,
                     BINS = NULL,
                     PED = NULL,
                     DCRS = NULL,
                     OUTPUT = NULL)


usage <- function() {
'Usage: Rscript annotate_denovo_cnv.R [OPTIONS] CALLSET BINS PED DCRS OUTPUT

CALLSET  Final callset produced by gCNV pipeline
BINS     Genomic intervals file used by gCNV pipeline
PED      Pedigree file
DCRS     A file with either the paths, one per line, to the dCR matrices or
         two tab-separated columns with the batch ID in the first column and
         path to the dCR matrix for that batch in the second column
OUTPUT   Where to write the output

Options:
  -h,--help               Print this message and exit.
  -n,--no-recal-freq      Do not recalibrate the variant frequency. If this
                          option is given, there must be a column named "sf" in
                          CALLSET with the frequency of each variant. [FALSE]
  -c,--hq-cols[=]<cols>   Comma-separated list of logical columns in CALLSET
                          that must all be TRUE for a call to be considered
                          high-quality. [PASS_SAMPLE,PASS_QS]
  -f,--max-freq[=]<freq>  Exclusive maximum variant frequency for a CNV to be
                          considered rare. Only rare de novo CNVs are annotated.
                          [0.01]
  -p,--cpus[=]<cpus>      Number of processors to use. [1]
  -s,--skip-allosomes     Do not call de novo variants on allosomes. [FALSE]
'
}

#' Write to stderr and exit
msg_and_exit <- function(msg, status = 1) {
    cat(msg, sep = "", file = stderr())
    quit(save = "no", status)
}

#' Parse a command-line option requiring an argument
#'
#' Three option types are supported:
#' * short option (e.g. -o <arg>)
#' * long option (e.g. --opt <arg>)
#' * long option with '=' (e.g. --opt=<arg>)
#'
#' @param x A character vector of all the command-line arguments i.e. as
#'   returned by `commandArgs(trailingOnly = TRUE)`.
#' @param m The match object returned by `regexec` from comparing the
#'   option against the option-specific regular expression.
#' @param i The index in `x` of the current option.
#' @returns A list of length two in which the first element is the argument of
#'   the parsed option and the second is the index of the next element of `x`
#'   parse. If an option does not provide an argument, the script will exit
#'   with an error message.
opt_arg <- function(x, m, i) {
    matches <- regmatches(x[[i]], m)[[1]]
    can_inc <- i + 1 <= length(x)
    opt <- ""
    if (length(matches) == 1) { # ^-a$
        opt <- matches[[1]]
    } else if (length(matches) == 3) { # ^(--long)(=(.*)?)?$
        if (!nzchar(matches[[3]])) {
            opt <- matches[[2]]
        } else {
            return(list(matches[[3]], i + 1))
        }
    } else { # ^((-a)|((--long)(=(.*)?)?))$
        if (!nzchar(matches[[7]])) {
            opt <- matches[[5]]
        } else {
            return(list(matches[[7]], i + 1))
        }
    }

    if (!can_inc) {
        msg_and_exit(sprintf("option %s requires an argument", opt))
    }

    list(x[[i + 1]], i + 2)
}

parse_args <- function() {
    argv <- commandArgs(trailingOnly = TRUE)
    if (length(argv) == 0) {
        msg_and_exit(usage(), 0)
    }

    opts <- DEFAULT_OPTS
    args <- DEFAULT_ARGS

    j <- 1
    i <- 1
    while (i <= length(argv)) {
        if (grepl("^((-h)|(--help))$", argv[[i]])) {
            msg_and_exit(usage(), 0)
        }
        if (grepl("^((-n)|(--no-recal-freq))$", argv[[i]])) {
            opts$recal_freq <- FALSE
            i <- i + 1
            next
        }

        if (all((m <- regexec("^((-c)|((--hq-cols)(=(.*)?)?))$", argv[[i]]))[[1]] != -1)) {
            x <- opt_arg(argv, m, i)
            opts$hq_cols <- trimws(strsplit(x[[1]], split = ",", fixed = TRUE)[[1]])
            if (any(!nzchar(opts$hq_cols))) {
                msg_and_exit("names of HQ columns cannot be empty\n")
            }
            i <- x[[2]]
            next
        }

        if (all((m <- regexec("^((-f)|((--max-freq)(=(.*)?)?))$", argv[[i]]))[[1]] != -1)) {
            x <- opt_arg(argv, m, i)
            opts$max_freq <- tryCatch(
                as.double(x[[1]]),
                warning = \(cnd) msg_and_exit(sprintf("could not convert '%s' to double\n", x[[1]])),
                error = \(cnd) msg_and_exit(sprintf("could not convert '%s' to double\n", x[[1]]))
            )
            if (opts$max_freq < 0 || opts$max_freq > 1) {
                msg_and_exit("maximum frequency must be between 0 and 1 inclusive\n")
            }
            i <- x[[2]]
            next
        }

        if (all((m <- regexec("^((-p)|((--cpus)(=(.*)?)?))$", argv[[i]]))[[1]] != -1)) {
            x <- opt_arg(argv, m, i)
            opts$cpus <- tryCatch(
                as.integer(x[[1]]),
                warning = \(cnd) msg_and_exit(sprintf("could not convert '%s' to integer\n", x[[1]])),
                error = \(cnd) msg_and_exit(sprintf("could not convert '%s' to integer\n", x[[1]]))
            )
            if (opts$cpus < 1) {
                msg_and_exit("number of processors must be greater than 0\n")
            }
            i <- x[[2]]
            next
        }

        if (grepl("^((-s)|(--skip-allosomes))$", argv[[i]])) {
            opts$allosomes <- FALSE
            i <- i + 1
            next
        }

        if (grepl("^-", argv[[i]])) {
            msg_and_exit(sprintf("unknown option %s\n", argv[[i]]))
        } else {
            if (j <= length(args)) {
                args[[j]] <- argv[[i]]
                j <- j + 1
                i <- i + 1
            } else {
                msg_and_exit("too many arguments\n")
            }
        }
    }

    missing_args <- Filter(is.null, args)
    if (length(missing_args) > 0) {
        msg_and_exit(sprintf("required argument(s) %s missing\n",
                             paste0(names(missing_args), collapse = ", ")))
    }

    append(args, opts)
}

args <- parse_args()

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
                           bins = NA_integer_,
                           miss_parents = FALSE)

format_sample_ids <- function(x) {
    is_long <- nchar(x) > 47L
    if (!any(is_long)) {
        return(x)
    }
    replace(x, is_long, paste0(substr(x[is_long], 1, 47), "..."))
}

# Check if a call to mclapply produced any errors.
# This is useful because mclapply will not propogate errors that occur in the
# jobs that it runs.
check_mclapply_errors <- function(x) {
    errs <- Filter(\(y) !is.null(attr(y, "condition")) && is(attr(y, "condition"), "error"), x)

    if (length(errs) > 0) {
        stop(attr(errs[[1]], "condition"), call. = FALSE)
    }

    x
}

# Compute the M, MF, MM, etc. table for a trio
dcr_evidence <- function(x, dcr_map) {
    region <- gregion(x$chr, x$start, x$end)
    trio_dcr <- get_trio_dcr(region,
                             c(x$sample, x$paternal_id, x$maternal_id),
                             c(x$batch, x$paternal_batch, x$maternal_batch),
                             dcr_map,
                             include_bg = TRUE,
                             keep_all_ranges = TRUE,
                             squeeze = TRUE,
                             reduce = TRUE)
    if (nrow(trio_dcr) == 0 || all(is.na(trio_dcr[[x$sample]]))) {
        stop(sprintf("CNV in '%s' at '%d' does not have supporting dCR",
                     x$sample,
                     to_string(region)))
    }

    is_chrx <- grepl("X", x$chr, fixed = TRUE)
    default <- DCR_EVIDENCE
    if (all(is.na(trio_dcr[[x$paternal_id]])) || all(is.na(trio_dcr[[x$maternal_id]]))) {
        default$miss_parents <- TRUE
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
        bins = length(mads),
        miss_parents = FALSE
    )

    rg
}

#' Predict the sex of samples from their chrX dCR values
#'
#' Prediction of sex is done based on the mean dCR of chrX for each sample.
#'
#' @param samples Sample IDs.
#' @param dcr `NULL` or a `data.table` of the dCR matrix of the samples in
#'   `samples` subset to the intervals on chromosome X.
#' @returns A `data.table` giving the mean dCR value of chrX and predicted
#'   sex of each sample. Samples that cannot be classified will have an
#'   `NA` for predicted sex. If `dcr` is `NULL`, every sample will be
#'   assigned an `NA` mean chrX dCR and predicted sex.
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

#' Predict the sex of samples
#'
#' Use chrX dCR values to predict chrX copy number of sample sex.
#'
#' @param x A `data.table` of sample and batch IDs.
#' @param dcr_map A hash map of batch ID to dCR matrix paths.
#' @param nproc Number of cpus to use.
#' @returns A `data.table` of all sample chrX copy number and predicted sex.
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
    old_dtthreads <- getDTthreads()
    setDTthreads(1)
    sexes <- mcmapply(f, dcr_groups, names(dcr_groups), SIMPLIFY = FALSE, mc.cores = nproc)
    setDTthreads(old_dtthreads)

    rbindlist(sexes)
}

#' Predict de novo CNVs based on missing overlap between parents and offspring
#'
#' Check all CNVs in a child for overlap with a CNV in either parent. The
#' overlap check is done with both genomic coordinates and bin coordinates. Any
#' CNV in the child that does not have at least 0.3 reciprocal overlap with a
#' CNV in a parent is considered a potential de novo.
#'
#' @param child `data.table` of calls in the child.
#' @param father `data.table` of calls in the father.
#' @param mother `data.table` of calls in the mother.
#' @param bins `GBins` object of the genomic bins used by gCNV.
#' @returns `data.table` of putative de novo calls.
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

#' Add parental batch information to a `data.table` of child calls
add_parent_batch <- function(child, all_calls) {
    batches <- unique(all_calls[, list(sample, batch)])
    colnames(batches) <- c("paternal_id", "paternal_batch")
    child <- batches[child, nomatch = NULL, on = "paternal_id",]
    colnames(batches) <- c("maternal_id", "maternal_batch")
    child <- batches[child, nomatch = NULL, on = "maternal_id",]

    child
}

#' Recalibrate child CNV frequency to match frequency in parents
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

#' Filter a `data.table` callset to high-quality calls
#'
#' @param calls `data.table`.
#' @param cols character vector of columns names in `calls` that must all be
#'   `TRUE` for a CNV call to be considered high-quality.
#' @returns The high-quality calls.
filter_hq_calls <- function(calls, cols) {
    missing_cols <- cols[!cols %in% colnames(calls)]
    if (length(missing_cols) > 0) {
        stop(sprintf("HQ columns not found in callset: %s",
                     paste0(missing_cols, collapse = ", ")),
             call. = FALSE)
    }

    x <- paste0(cols, collapse = " & ")

    if (length(cols) == 1) {
        filter <- call("==", as.name(x), TRUE)
        return(calls[eval(filter)])
    }

    calls[eval(parse(text = x))]
}

#' Autosomal de novo CNV workflow
autosome_denovo <- function(calls, bins, ped, dcrs, recal_freq, hq_cols, max_freq, nproc) {
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
    if (recal_freq) {
        calls <- recal_cnv_freq(calls, ped)
    } else if (!"sf" %in% colnames(calls)) {
        stop("not recomputing variant frequency, but no 'sf' column in callset",
             call. = FALSE)
    }
    setkey(calls, sample)

    # Filter to high quality calls and merge in pedigree ----------------------
    child_calls <- filter_hq_calls(calls, hq_cols)
    child_calls <- child_calls[sf < get("max_freq", envir = parent.frame(3))]
    child_calls <- child_calls[ped, on = c(sample = "sample_id"), nomatch = NULL]

    # Get preliminary de novo calls based on CNV overlap ----------------------
    log_info("making initial de novo predictions based on CNV overlap")
    dn <- denovo_from_ovp(child_calls,
                          calls[sample %in% ped$paternal_id],
                          calls[sample %in% ped$maternal_id],
                          bins)
    log_info(sprintf("found %d de novo events based on overlap", nrow(dn)))
    if (nrow(dn) == 0) {
        return(dn)
    }

    # Gather dCR evidence -----------------------------------------------------
    log_info("gathering dCR evidence")
    dn <- add_parent_batch(dn, calls)
    old_dtthreads <- getDTthreads()
    # Ensure that data.table runs single-threaded inside a mclapply call.
    setDTthreads(1)
    rg <- mclapply(seq_len(nrow(dn)),
                   \(i) dcr_evidence(as.list(dn[i, ]), dcrs),
                   mc.cores = nproc) |>
        suppressWarnings() |>
        check_mclapply_errors() |>
        rbindlist()
    setDTthreads(old_dtthreads)

    # Regenotype --------------------------------------------------------------
    log_info("regenotyping")
    rg <- cbind(rg, dn[, list(CN, chr, inheritance, svtype)])

    rg[miss_parents == TRUE, inheritance := "fail_miss_parents"]
    rg[abs(M - CN) > 0.175 & CN <= 3, inheritance := "fail_M"]
    rg[MD < 0.7, inheritance := "fail_MD"]
    rg[abs(M - CN) > 0.175 & CN <= 3 & MD < 0.7, inheritance := "fail_M_MD"]
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
    rg[bins_fail / bins >= 0.5 & bins < 100, inheritance := "fail_bad_bins"]

    dn[, inheritance := rg$inheritance]

    dn
}

#' chrX de novo CNV workflow
chrx_denovo <- function(calls, bins, ped, dcrs, recal_freq, hq_cols, max_freq, nproc) {
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
    if (recal_freq) {
        calls <- recal_cnv_freq(calls, ped)
    } else if (!"sf" %in% colnames(calls)) {
        stop("not recomputing variant frequency, but no 'sf' column in callset",
             call. = FALSE)
    }

    # Filter to high quality calls and merge in pedigree ----------------------
    child_calls <- filter_hq_calls(calls, hq_cols)
    child_calls <- child_calls[sf < get("max_freq", envir = parent.frame(3)) & grepl("X", chr)]
    child_calls <- child_calls[ped, on = c(sample = "sample_id"), nomatch = NULL]

    # Get preliminary de novo calls based on CNV overlap ----------------------
    log_info("making initial de novo predictions based on CNV overlap")
    dn <- denovo_from_ovp(child_calls,
                          calls[sample %in% ped$paternal_id & grepl("X", chr)],
                          calls[sample %in% ped$maternal_id & grepl("X", chr)],
                          bins)
    log_info(sprintf("found %d de novo events based on overlap", nrow(dn)))
    if (nrow(dn) == 0) {
        return(dn)
    }

    # Gather dCR evidence -----------------------------------------------------
    log_info("gathering dCR evidence")
    dn <- add_parent_batch(dn, calls)
    old_dtthreads <- getDTthreads()
    setDTthreads(1)
    rg <- mclapply(seq_len(nrow(dn)),
                   \(i) dcr_evidence(as.list(dn[i, ]), dcrs),
                   mc.cores = nproc) |>
        suppressWarnings() |>
        check_mclapply_errors() |>
        rbindlist()
    setDTthreads(old_dtthreads)

    # Regenotype --------------------------------------------------------------
    log_info("regenotyping")
    rg <- cbind(rg, dn[, list(sample, chr, CN, inheritance, svtype, paternal_id, maternal_id)])
    rg <- pass_sex[, list(sample, sex)][rg, on = "sample"]
    rg <- pass_sex[, list(paternal_id = sample, paternal_sex = sex)][rg, on = "paternal_id"]
    rg <- pass_sex[, list(maternal_id = sample, maternal_sex = sex)][rg, on = "maternal_id"]

    # Check sex of parents ----------------------------------------------------
    rg[miss_parents == TRUE, inheritance := "fail_miss_parents"]
    rg[paternal_sex != 1 | maternal_sex != 2, inheritance := "fail_parental_sex"]

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

# Read inputs -----------------------------------------------------------------
setDTthreads(args$cpus)

log_info("reading callset")
raw_calls <- read_callset(args$CALLSET)
is_hg19 <- any(c(as.character(1:22), "X", "Y") %in% raw_calls$chr)

log_info("reading bins")
bins <- read_gcnv_bins(args$BINS, reduce = is_hg19)

log_info("reading pedigree")
ped <- suppressWarnings(read_pedigree(args$PED))

log_info("reading dCR paths")
dcrs <- read_dcr_list(args$DCRS)

log_info("calling autosome de novo CNVs")
dn_auto <- autosome_denovo(raw_calls, bins, ped, dcrs,
                           recal_freq = args$recal_freq,
                           hq_cols = args$hq_cols,
                           max_freq = args$max_freq,
                           nproc = args$cpus)
log_info("completed calling autosome de novo CNVs")

if (args$allosomes) {
    log_info("calling chrX de novo CNVs")
    dn_chrx <- chrx_denovo(raw_calls, bins, ped, dcrs,
                           recal_freq = args$recal_freq,
                           hq_cols = args$hq_cols,
                           max_freq = args$max_freq,
                           nproc = args$cpus)
    log_info("completed calling chrX de novo CNVs")
    common_cols <- intersect(colnames(dn_auto), colnames(dn_chrx))
    dn_all <- rbind(dn_auto[, ..common_cols], dn_chrx[, ..common_cols])
} else {
    log_info("skipping de novo CNVs on allosomes")
    dn_all <- dn_auto
}

# Write output to file --------------------------------------------------------
dn_all <- ped[, c("family_id", "sample_id")][dn_all,
                                             on = c(sample_id = "sample"),
                                             mult = "first"]
setnames(dn_all, "sample_id", "sample")

log_info("writing output")
write.table(dn_all,
            file = args$OUTPUT,
            sep = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)
log_info("done")
