# chrX/Y de novo CNV annotation pipeline
#
# Annotate de novo CNVs from a gCNV callset.
# Usage: Rscript annotate_denovo_cnv.R CALLSET BINS PED DCRS NPROC OUTPUT
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

get_denovo_evidence <- function(call_info, dcr_map) {
    cid <- call_info$sample
    pid <- call_info$paternal_id
    mid <- call_info$maternal_id
    cbatch <- call_info$batch
    pbatch <- call_info$paternal_batch
    mbatch <- call_info$maternal_batch
    region <- gregion(call_info$chr, call_info$start, call_info$end)

    coord_cols <- c("chr", "start", "end")
    trio_ids <- c(cid, pid, mid)
    batches <- c(cbatch, pbatch, mbatch)
    dcr_groups <- split(trio_ids, batches)
    dcrs <- vector(mode = "list", length = length(dcr_groups))
    for (i in seq_along(dcr_groups)) {
        batch_i <- names(dcr_groups)[[i]]
        samples_i <- dcr_groups[[i]]
        dcr_paths <- gethash(dcr_map, batch_i)
        if (is.null(dcr_paths)) {
            return(regenotype(NULL))
        }
        dcr <- tryCatch(
            get_samples_dcr(
                region, dcr_paths, samples_i, include_bg = cid %in% samples_i
            ),
            error = function(cnd) NULL
        )
        if (is.null(dcr) || nrow(dcr) == 0) {
            return(regenotype(NULL))
        }

        cols <- setdiff(colnames(dcr), c(trio_ids[!trio_ids %in% samples_i], coord_cols))
        dcr_mat <- data.matrix(dcr[cols])
        rownames(dcr_mat) <- paste0(dcr$chr, ":", dcr$start, "-", dcr$end)
        dcrs[[i]] <- dcr_mat
    }

    bins <- unique(unlist(lapply(dcrs, rownames)))
    if (length(bins) == 0) {
        return(regenotype(NULL))
    }
    sample_ids <- lapply(dcrs, colnames)
    trio_dcr <- matrix(nrow = length(bins), ncol = sum(lengths(sample_ids)))
    rownames(trio_dcr) <- bins
    colnames(trio_dcr) <- unlist(sample_ids)
    for (mat in dcrs) {
        trio_dcr[rownames(mat), colnames(mat)] <- mat
    }

    bg_samples <- colnames(trio_dcr)[!colnames(trio_dcr) %in% trio_ids]
    trio_dcr <- trio_dcr[, c(trio_ids, bg_samples), drop = FALSE]

    regenotype(trio_dcr)
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

predict_sex <- function(x, dcr_map) {
    x <- distinct(x)
    hg38_region <- gregion("chrX", 1, TABIX_MAX_SEQLEN)
    hg19_region <- gregion("X", 1, TABIX_MAX_SEQLEN)

    dcr_groups <- split(x$sample, x$batch)
    sexes <- vector(mode = "list", length = length(dcr_groups))
    for (i in seq_along(dcrs)) {
        batch_i <- names(dcr_groups)[[i]]
        samples_i <- dcr_groups[[i]]
        dcr_paths <- gethash(dcr_map, batch_i)
        if (is.null(dcr_paths)) {
            sexes[[i]] <- sex_from_dcr(samples_i, NULL)
            next
        }

        dcr <- tryCatch(
            get_samples_dcr(hg38_region, dcr_paths, samples_i),
            dcr_parse_error = function(cnd) cnd,
            error = function(cnd) cnd
        )

        if (inherits(dcr, "dcr_parse_error")) (
            dcr <- tryCatch(
                get_samples_dcr(hg19_region, dcr_paths, samples_i),
                error = function(cnd) NULL
            )
        )

        if (is.null(dcr) || inherits(dcr, "error")) {
            sexes[[i]] <- sex_from_dcr(samples_i, NULL)
        } else {
            sexes[[i]] <- sex_from_dcr(samples_i, dcr)
        }
    }

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
ped <- read_pedigree(ped_path) |>
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

# Check sex of all samples ----------------------------------------------------
log_info("validating sex of samples")
predicted_sex <- predict_sex(raw_calls[, c("sample", "batch")], dcrs)
fail_dcr <- filter(predicted_sex, is.na(chrX_CN))
if (nrow(fail_dcr) > 0) {
    log_warn(sprintf("%s samples failed chrX CN calling", nrow(fail_dcr)))
}

fail_sex <- filter(predicted_sex, !is.na(chrX_CN) & is.na(predicted_sex))
if (nrow(fail_sex) > 0) {
    max_width <- max(nchar(fail_sex$sample))
    header_fmt <- paste0("%", max_width, "s  chrX_CN\n")
    body_fmt <- paste0("%", max_width, "s  %7.3f\n")
    log_warn(c(
        sprintf("%s samples have out-of-range chrX CN", nrow(fail_sex)),
        sprintf(header_fmt, "sample"),
        sprintf(body_fmt, fail_sex$sample, fail_sex$chrX_CN)
    ))
}

sex_match_df <- filter(predicted_sex, !is.na(predicted_sex)) |>
    left_join(ped, by = join_by(sample == sample_id))
missing_sex_cnt <- sum(is.na(sex_match_df$sex))
if (missing_sex_cnt > 0) {
    log_warn(c(
        paste0("discarding ", missing_sex_cnt, " samples from callset without pedigree"),
        sprintf("%s\n", sex_match_df[is.na(sex_match_df$sex), ]$sample)
    ))
}

sex_match_df <- filter(sex_match_df, !is.na(sex))
mismatch_sex <- filter(pass_sex, predicted_sex != sex)
if (nrow(mismatch_sex) > 0) {
    log_warn(
        paste0("discarding ", nrow(mismatch_sex), " samples with mismatched sex")
    )
    cat(
        sprintf("%s\n", pass_sex[pass_sex$sex != pass_sex$predicted_sex, ]$sample),
        file = stderr(),
        sep = ""
    )
}

pass_sex <- filter(pass_sex, predicted_sex == sex)
log_info(sprintf("continuing analysis with %d samples passing sex check", nrow(pass_sex)))

# Filter the callset and pedigree ---------------------------------------------
raw_calls <- filter(raw_calls, sample %in% pass_sex$sample)
trio_ped <- filter(ped, sample_id %in% raw_calls$sample) |>
    filter(paternal_id %in% raw_calls$sample) |>
    filter(maternal_id %in% raw_calls$sample)

# Recalibrate offspring CNV frequency to parents' -----------------------------
log_info("making initial de novo predictions based on CNV overlap")
parent_ids <- unique(c(ped$paternal_id, ped$maternal_id))
parent_cnvs <- raw_calls |>
    filter(sample %in% parent_ids) |>
    group_by(variant_name) |>
    summarize(sc = n()) |>
    mutate(sf = sc / length(parent_ids))

raw_calls <- select(raw_calls, !c(sc, sf)) |>
    left_join(parent_cnvs, by = "variant_name", relationship = "many-to-one") |>
    mutate(sc = replace(sc, is.na(sc), 0), sf = replace(sf, is.na(sf), 0)) |>
    arrange(sample, chr, start)

# Filter calls and merge in pedigree ------------------------------------------
child_calls <- filter(
    raw_calls,
    PASS_SAMPLE & PASS_QS & sf < 0.01 & grepl("X|Y", chr)
) |>
    inner_join(
        ped,
        by = join_by(sample == sample_id),
        multiple = "first"
    ) |>
    select(!c(family_id, phenotype))

# Get preliminary de novo calls based on overlap ------------------------------
gr_c <- df_to_gr(child_calls, cnv = TRUE)

paternal_calls <- filter(raw_calls, sample %in% ped$paternal_id, grepl("X|Y", chr))
maternal_calls <- filter(raw_calls, sample %in% ped$maternal_id, grepl("X|Y", chr))
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

# Regenotype
# test CNV overlap between child and parents
