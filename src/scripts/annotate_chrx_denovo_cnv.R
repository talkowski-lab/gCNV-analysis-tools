# chrX de novo CNV annotation pipeline
#
# Annotate de novo CNVs from a gCNV callset.
# Usage: Rscript annotate_chrx_denovo_cnv.R CALLSET BINS PED DCRS NPROC OUTPUT
# * CALLSET - Final callset produced by gCNV pipeline
# * BINS    - Genomic bins file used by the gCNV pipeline
# * PED     - Pedigree file
# * DCRS    - List of paths, one per line, to the dCR matrices
# * NPROC   - Number of processors to use
# * OUTPUT  - Where to write the output

XX_MIN <- 1.9
XX_MAX <- 2.1
XY_MIN <- 0.9
XY_MAX <- 1.1
TABIX_MAX_SEQLEN <- 536870912L

get_coverage <- function(x, y, relation = c("paternal", "maternal")) {
    relation <- match.arg(relation)
    ol <- findOverlaps(x, y)
    x_mcols <- mcols(x[queryHits(ol)])
    y_mcols <- mcols(y[subjectHits(ol)])
    if (relation == "paternal") {
        ol <- ol[x_mcols[["paternal_id"]] == y_mcols[["sample"]]]
    } else {
        ol <- ol[x_mcols[["maternal_id"]] == y_mcols[["sample"]]]
    }
    int <- pintersect(x[queryHits(ol)], y[subjectHits(ol)])
    cov_hit <- by(width(int) / width(x[queryHits(ol)]), queryHits(ol), sum)
    out <- double(length = length(x))
    out[as.integer(names(cov_hit))] <- cov_hit

    out
}

get_denovo_evidence <- function(x, dcr_map) {
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

    if (is.null(trio_dcr) || nrow(trio_dcr) == 0) {
        return(regenotype(NULL))
    }

    regenotype(trio_dcr[!colnames(trio_dcr) %in% c("chr", "start", "end")])
}

regenotype <- function(dcr) {
    default_rg <- tibble(
        M = NA_real_,
        MF = NA_real_,
        MM = NA_real_,
        MD = NA_real_,
        bins_fail = NA_integer_,
        bins = NA_integer_
    )
    if (is.null(dcr)) {
        return(default_rg)
    }

    cutoff <- 1
    mads <- apply(dcr, 1, mad)
    mads_fail <- sum(mads >= cutoff, na.rm = TRUE)
    if (mads_fail == length(mads)) {
        default_rg$bins_fail <- mads_fail
        default_rg$bins <- length(mads)
        return(default_rg)
    }

    means <- apply(dcr[mads < cutoff, , drop = FALSE], 2, mean, na.rm = TRUE)
    # Assume child in first column, father in second, and mother in third
    md <- min(c(abs(means[[1]] - means[[2]]), abs(means[[1]] - means[[3]])))
    rg <- tibble(
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
        return(tibble(sample = samples, chrX_CN = NA_real_, predicted_sex = NA_integer_))
    }
    mean_x <- select(dcr, !c(chr, start, end)) |>
        apply(2, mean)

    sexes <- tibble(sample = samples, chrX_CN = mean_x[samples]) |>
        mutate(
            predicted_sex = case_when(
                chrX_CN >= XY_MIN & chrX_CN <= XY_MAX ~ 1L,
                chrX_CN >= XX_MIN & chrX_CN <= XX_MAX ~ 2L,
                .default = NA_integer_
            )
        )

    sexes
}

predict_sex <- function(x, dcr_map, nproc = 1L) {
    x <- distinct(x)
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
            return(sex_from_dcr(samples, dcr))
        }
    }
    sexes <- mcmapply(f, dcr_groups, names(dcr_groups), SIMPLIFY = FALSE, mc.cores = nproc)

    bind_rows(sexes)
}

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

suppressPackageStartupMessages(library(gelpers))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(dplyr))

# Read inputs -----------------------------------------------------------------
log_info("reading callset")
raw_calls <- read_callset(callset_path) |>
    as_tibble()

is_hg19 <- any(c(as.character(1:22), "X", "Y") %in% raw_calls$chr)
log_info("reading bins")
bins <- read_gcnv_bins(bins_path, reduce = is_hg19)

log_info("reading pedigree")
ped <- suppressWarnings(read_pedigree(ped_path)) |>
    as_tibble() |>
    filter(sample_id %in% raw_calls$sample)

log_info("reading dCR paths")
dcrs <- read_dcr_list(dcrs_path)

# Filter pedigree for missing sex ---------------------------------------------
missing_sex <- sum(is.na(ped$sex))
if (missing_sex > 0) {
    log_warn(
        paste0("discarding ", missing_sex, " samples from pedigree with missing sex")
    )
}
ped <- filter(ped, !is.na(sex))

# Predict the sex of all samples ----------------------------------------------
log_info("validating sex of samples")
predicted_sex <- predict_sex(raw_calls[, c("sample", "batch")], dcrs, nproc)
fail_dcr <- filter(predicted_sex, is.na(chrX_CN))
if (nrow(fail_dcr) > 0) {
    log_warn(sprintf("%s samples failed chrX CN calling", nrow(fail_dcr)))
}

# Report samples with chrX copy number outside range --------------------------
fail_sex <- filter(predicted_sex, !is.na(chrX_CN) & is.na(predicted_sex))
if (nrow(fail_sex) > 0) {
    max_width <- max(nchar(fail_sex$sample))
    header_fmt <- paste0("%", max_width, "s  chrX_CN")
    body_fmt <- paste0("%", max_width, "s  %7.3f")
    log_warn(c(
        sprintf("%s samples have out-of-range chrX CN", nrow(fail_sex)),
        sprintf(header_fmt, "sample"),
        sprintf(body_fmt, fail_sex$sample, fail_sex$chrX_CN)
    ))
}

# Remove callset samples without a pedigree -----------------------------------
sex_match_df <- filter(predicted_sex, !is.na(predicted_sex)) |>
    left_join(ped, by = join_by(sample == sample_id))
missing_sex_cnt <- sum(is.na(sex_match_df$sex))
if (missing_sex_cnt > 0) {
    log_warn(c(
        paste0("discarding ", missing_sex_cnt, " samples from callset without pedigree"),
        sprintf("%s", sex_match_df[is.na(sex_match_df$sex), ]$sample)
    ))
}

# Remove samples whose predicted sex doesn't match pedigree sex ---------------
sex_match_df <- filter(sex_match_df, !is.na(sex))
mismatch_sex <- filter(sex_match_df, predicted_sex != sex)
if (nrow(mismatch_sex) > 0) {
    max_width <- max(nchar(mismatch_sex$sample))
    header_fmt <- paste0("%", max_width, "s  chrX_CN  predicted_sex  pedigree_sex")
    body_fmt <- paste0("%", max_width, "s  %7.3f  %13s  %12s")
    log_warn(c(
        paste0("discarding ", nrow(mismatch_sex), " samples with mismatched sex"),
        sprintf(header_fmt, "sample"),
        sprintf(body_fmt,
                mismatch_sex$sample,
                mismatch_sex$chrX_CN,
                ifelse(mismatch_sex$predicted_sex == 1, "male", "female"),
                ifelse(mismatch_sex$sex == 1, "male", "female"))
    ))
}

pass_sex <- filter(sex_match_df, predicted_sex == sex)
log_info(sprintf("continuing analysis with %d samples passing sex check", nrow(pass_sex)))

# Filter the callset and pedigree ---------------------------------------------
candidate_calls <- filter(raw_calls, sample %in% pass_sex$sample)
trio_ped <- filter(ped, sample_id %in% raw_calls$sample) |>
    filter(paternal_id %in% raw_calls$sample) |>
    filter(maternal_id %in% raw_calls$sample)

# Recalibrate offspring CNV frequency to parents' -----------------------------
log_info("making initial de novo predictions based on CNV overlap")
parent_ids <- unique(c(ped$paternal_id, ped$maternal_id))
parent_cnvs <- candidate_calls |>
    filter(sample %in% parent_ids) |>
    group_by(variant_name) |>
    summarize(sc = n()) |>
    mutate(sf = sc / length(parent_ids))

candidate_calls <- select(candidate_calls, !c(sc, sf)) |>
    left_join(parent_cnvs, by = "variant_name", relationship = "many-to-one") |>
    mutate(sc = replace(sc, is.na(sc), 0), sf = replace(sf, is.na(sf), 0)) |>
    arrange(sample, chr, start)

# Filter calls and merge in pedigree ------------------------------------------
child_calls <- filter(
    candidate_calls,
    PASS_SAMPLE & PASS_QS & sf < 0.01 & grepl("X", chr)
) |>
    inner_join(
        trio_ped,
        by = join_by(sample == sample_id),
        multiple = "first"
    ) |>
    select(!c(family_id, phenotype))

# Get preliminary de novo calls based on overlap ------------------------------
gr_c <- df_to_gr(child_calls, cnv = TRUE)

paternal_calls <- filter(candidate_calls, sample %in% ped$paternal_id, grepl("X", chr))
maternal_calls <- filter(candidate_calls, sample %in% ped$maternal_id, grepl("X", chr))
gr_p <- df_to_gr(paternal_calls, cnv = TRUE)
gr_m <- df_to_gr(maternal_calls, cnv = TRUE)

child_calls <- mutate(
    child_calls,
    cov_p = get_coverage(gr_c, gr_p, "paternal"),
    cov_m = get_coverage(gr_c, gr_m, "maternal")
)

gr_c_bs <- toBinSpace(gr_c, bins)
gr_p_bs <- toBinSpace(gr_p, bins)
gr_m_bs <- toBinSpace(gr_m, bins)

child_calls <- mutate(
    child_calls,
    cov_p_bs = get_coverage(gr_c_bs, gr_p_bs, "paternal"),
    cov_m_bs = get_coverage(gr_c_bs, gr_m_bs, "maternal")
)

child_calls <- child_calls |>
    mutate(
        inheritance = case_when(
            (cov_m_bs >= 0.3 | cov_m >= 0.3) & (cov_p_bs >= 0.3 | cov_p >= 0.3) ~ "biparental",
            cov_p_bs >= 0.3 | cov_p >= 0.3 ~ "paternal",
            cov_m_bs >= 0.3 | cov_m >= 0.3 ~ "maternal",
            cov_m_bs < 0.3 & cov_m < 0.3 & cov_p_bs < 0.3 & cov_p < 0.3 ~ "denovo",
            .default = "inherited"
        )
    )
dn <- filter(child_calls, inheritance == "denovo")
log_info(sprintf("found %d putative denovo CNVs based on overlap", nrow(dn)))

# Regenotype ------------------------------------------------------------------
batch_tbl <- select(candidate_calls, sample, batch) |>
    distinct() |>
    mutate(
        paternal_id = sample,
        maternal_id = sample,
        paternal_batch = batch,
        maternal_batch = batch
    )
dn <- select(batch_tbl, paternal_id, paternal_batch) |>
    inner_join(dn, by = "paternal_id", relationship = "one-to-many")
dn <- select(batch_tbl, maternal_id, maternal_batch) |>
    inner_join(dn, by = "maternal_id", relationship = "one-to-many")

log_info("gathering dCR evidence")
rg_info <- mclapply(
    seq_len(nrow(dn)),
    \(i) get_denovo_evidence(as.list(dn[i, ]), dcrs),
    mc.cores = nproc
) |>
    bind_rows()

log_info("regenotyping")
rg_info <- rg_info |>
    mutate(
        sample = dn$sample,
        chr = dn$chr,
        CN = dn$CN,
        inheritance = dn$inheritance,
        svtype = dn$svtype,
        paternal_id = dn$paternal_id,
        maternal_id = dn$maternal_id
    ) |>
    inner_join(select(pass_sex, sample, sex), by = "sample") |>
    left_join(rename(select(pass_sex, sample, sex),
                     paternal_id = "sample",
                     paternal_sex = "sex"),
              by = "paternal_id") |>
    left_join(rename(select(pass_sex, sample, sex),
                     maternal_id = "sample",
                     maternal_sex = "sex"),
              by = "maternal_id")

# Check sex of parents --------------------------------------------------------
idx <- rg_info$paternal_sex != 1 | rg_info$maternal_sex != 2
rg_info[idx, "inheritance"] <- "fail_parental_sex"

# Check offspring is consistent with SV type and sex --------------------------
# sex  svtype         expected_M
#  XX     DEL  0.75 <= M <= 1.25
#  XY     DEL  0.00 <= M <= 0.25
#  XX     DUP  2.75 <= M
#  XY     DUP  1.75 <= M
idx <- with(rg_info, sex == 2L & svtype == "DEL" & abs(M - 1) > 0.25)
rg_info[idx, "inheritance"] <- "fail_M"
idx <- with(rg_info, sex == 1L & svtype == "DEL" & M > 0.25)
rg_info[idx, "inheritance"] <- "fail_M"
idx <- with(rg_info, sex == 2L & svtype == "DUP" & M < 2.75)
rg_info[idx, "inheritance"] <- "fail_M"
idx <- with(rg_info, sex == 1L & svtype == "DUP" & M < 1.75)
rg_info[idx, "inheritance"] <- "fail_M"

# Check for CNV evidence in parents -------------------------------------------
# We check if a parent has a possible CNV in the same region as a child CNV
# by comparing the parent chrX ploidy to the parent dCR value for the region.
# For example, if a parent has 1 chrX, then it is expected that the dCR value
# for any region on chrX will be about 1. If either parent has a suspected CNV
# that matches the SV type of the CNV in the child, then the CNV in the child
# is likely inherited. A threshold of max deviation 0.15 was determined by
# plotting deviations for random sized subsets of some dCR matrices.
max_chrx_dev <- 0.15
idx <- with(rg_info, svtype == "DEL" & paternal_sex - MF > max_chrx_dev)
rg_info[idx, "inheritance"] <- "fail_MF"
idx <- with(rg_info, svtype == "DEL" & maternal_sex - MM > max_chrx_dev)
rg_info[idx, "inheritance"] <- "fail_MM"
idx <- with(rg_info,
            svtype == "DEL" & paternal_sex - MF > max_chrx_dev & maternal_sex - MM > max_chrx_dev)
rg_info[idx, "inheritance"] <- "fail_MF_MM"

idx <- with(rg_info, svtype == "DUP" & MF - paternal_sex > max_chrx_dev)
rg_info[idx, "inheritance"] <- "fail_MF"
idx <- with(rg_info, svtype == "DUP" & MM - maternal_sex > max_chrx_dev)
rg_info[idx, "inheritance"] <- "fail_MM"
idx <- with(rg_info,
            svtype == "DUP" & MF - paternal_sex > max_chrx_dev & MM - maternal_sex > max_chrx_dev)
rg_info[idx, "inheritance"] <- "fail_MF_MM"

# Write annotations to file ---------------------------------------------------
dn$inheritance <- rg_info$inheritance
out <- inner_join(
    dn, select(trio_ped, sample_id, family_id),
    by = join_by(sample == sample_id),
    keep = FALSE,
    multiple = "first"
)

log_info("writing output")
write.table(
    out, file = output, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)
log_info("done")
