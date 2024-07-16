# de novo CNV annotation pipeline
#
# Annotate de novo CNVs from a gCNV callset.
# Usage: Rscript annotate_denovo_cnv.R CALLSET BINS PED DCRS NPROC OUTPUT
# * CALLSET - Final callset produced by gCNV pipeline
# * BINS    - Genomic bins file used by the gCNV pipeline
# * PED     - Pedigree file
# * DCRS    - List of paths, one per line, to the dCR matrices
# * NPROC   - Number of processors to use
# * OUTPUT  - Where to write the output

# Functions -------------------------------------------------------------------
# Compute the fraction of each CNV in x that is overlapped by the parental CNVs
# in y.
get_coverage <- function(x, y, relation = c("paternal", "maternal")) {
    relation <- match.args(relation)
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

# Get the M, MF, MM, etc. data.frame for a suspected de novo call.
get_denovo_evidence <- function(call_info, dcr_map) {
    cid <- call_info$sample
    pid <- call_info$paternal_id
    mid <- call_info$maternal_id
    cbatch <- call_info$batch
    pbatch <- call_info$paternal_batch
    mbatch <- call_info$maternal_batch
    region <- gregion(call_info$chr, call_info$start, call_info$end)

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
        if (is.null(dcr)) {
            return(regenotype(NULL))
        }

        dcrs[[i]] <- dcr
    }

    coord_cols <- c("chr", "start", "end")
    trio_dcr <- Reduce(\(x, y) merge(x, y, by = coord_cols), dcrs)
    if (nrow(trio_dcr) == 0) {
        return(regenotype(NULL))
    }
    bg_samples <- colnames(trio_dcr)[!colnames(trio_dcr) %in% c(coord_cols, trio_ids)]
    trio_dcr <- trio_dcr[, c(coord_cols, trio_ids, bg_samples), with = FALSE]

    regenotype(trio_dcr)
}

# Compute the M, MF, MM, etc. table from the trio dCR of a de novo call.
regenotype <- function(dcr) {
    default_rg <- data.table(
        M = NA_real_,
        MF = NA_real_,
        MM = NA_real_,
        MD = NA_real_,
        bins_fail = NA_integer_,
        bins = NA_integer_
    )
    if (is.null(dcr)) {
        return(default_rt)
    }

    cutoff <- if (any(grepl("chr[XY]", dcr$chr))) 1 else 0.5
    mat <- data.matrix(dcr[, !c("chr", "start", "end")])
    mads <- apply(mat, 1, mad)
    mads_fail <- sum(mads >= cutoff)
    if (mads_fail == length(mads)) {
        default_rg[1, `:=`(bins_fail = mads_fail, bins = length(mads))]
        return(default_rg)
    }

    means <- apply(mat[mads < cutoff, ,drop = FALSE], 2, mean)
    # Assume child in first column, father in second, and mother in third
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

library(data.table)
library(gelpers)
library(GenomicRanges)
library(parallel)

# Read inputs -----------------------------------------------------------------
raw_calls <- read_callset(callset_path)
bins <- read_gcnv_bins(bins_path)
ped <- read_pedigree(ped_path)
dcrs <- read_dcr_list(dcrs_path)

# Recalibrate offspring CNV frequency to parents' -----------------------------
parent_ids <- unique(c(ped$paternal_id, ped$maternal_id))
parent_cnvs <- raw_calls[sample %in% parent_ids,
                         list(sc = .N),
                         by = "variant_name"]
parent_cnvs[, sf := sc / length(parent_ids)]
raw_calls <- parent_cnvs[raw_calls[, !c("sc", "sf")], on = "variant_name"]
raw_calls[is.na(sc), sc := 0]
raw_calls[is.na(sf), sf := 0]
setkey(raw_calls, sample, chr, start)

# Filter calls and merge in pedigree ------------------------------------------
child_calls <- raw_calls[PASS_SAMPLE & PASS_QS & sf < 0.01]
child_calls <- child_calls[sample %in% ped[["sample_id"]]]
child_calls <- child_calls[ped[, list(sample_id, paternal_id, maternal_id)],
                           on = c(sample = "sample_id"),
                           nomatch = NULL]
child_calls <- child_calls[!(is.na(paternal_id) | is.na(maternal_id)), ]

# Get preliminary de novo calls based on overlap ------------------------------
gr_c <- df_to_gr(child_calls)
mcols(gr_c) <- gr_c[, list(paternal_id, maternal_id)]

paternal_calls <- raw_calls[sample %in% ped$paternal_id, ]
mcols(gr_p)$sample <- paternal_calls$sample
maternal_calls <- raw_calls[sample %in% ped$maternal_id, ]
mcols(gr_m)$sample <- maternal_calls$sample

child_calls[, `:=`(
    cov_p = get_coverage(gr_c, gr_p, "paternal"),
    cov_m = get_coverage(gr_c, gr_m, "maternal")
)]

gr_c_bs <- toBinSpace(gr_c, bins)
gr_p_bs <- toBinSpace(gr_p, bins)
gr_m_bs <- toBinSpace(gr_m, bins)

child_calls[, `:=`(
    cov_p_bs = get_coverage(gr_c_bs, gr_p_bs, "paternal"),
    cov_m_bs = get_coverage(gr_c_bs, gr_m_bs, "maternal")
)]

child_calls[,
  inheritance := fcase(
    (cov_m_bs >= 0.3 | cov_m >= 0.3) & (cov_p_bs >= 0.3 | cov_p >= 0.3), "biparental",
    cov_p_bs >= 0.3 | cov_p >= 0.3, "paternal",
    cov_m_bs >= 0.3 | cov_m >= 0.3, "maternal",
    cov_m_bs < 0.3 & cov_m < 0.3 & cov_p_bs < 0.3 & cov_p < 0.3, "denovo",
    default = "inherited"
  )
]
dn <- child_calls[inheritance == "denovo"]

# Gather dCR evidence ---------------------------------------------------------
rg_info <- mclapply(
    seq_len(nrow(dn)),
    \(i) get_denovo_evidence(as.list(dn[i, ]), dcrs),
    mc.cores = nproc
) |>
    rbindlist()

rg_info[, `:=`(CN = dn$CN, chr = dn$chr, inheritance = dn$inheritance, svtype = dn$svtype)]

# Refine de novo calls with dCR evidence -------------------------------------
rg_info[is.na(MF) | is.na(MM), inheritance := "fail_miss_parents"]
rg_info[abs(M - CN) > 0.175 & !grepl("X|Y", chr) & CN <= 3, inheritance := "fail_M"]
rg_info[MD < 0.7 & !grepl("X|Y", chr), inheritance := "fail_MD"]
rg_info[
  (abs(M - CN) > 0.175 & !grepl("X|Y", chr) & CN <= 3) &
  (MD < 0.7 & !grepl("X|Y", chr)),
  inheritance := "fail_M_MD"
]
rg_info[bins_fail / bins >= 0.5 & bins < 100, inheritance := "fail_bad_bins"]

rg_info[
  grepl("X", chr) & (
    (MF > 1.3 & svtype == "DUP") |
    (MF < 0.7 & svtype == "DEL") |
    (MM > 2.3 & svtype == "DUP") |
    (MM < 1.7 & svtype == "DEL")),
  inheritance := "fail_chrX"
]
rg_info[
  grepl("Y", chr) & (
    (MF > 1.3 & svtype == "DUP") |
    (MF < 0.7 & svtype == "DEL") |
    (MM > 0.3 & svtype == "DUP")),
  inheritance := "fail_chrY"
]
rg_info[
  CN > 2 & !grepl("X|Y", chr) & inheritance == "denovo" & (MF > 2.5 | MM > 2.5),
  inheritance := "inherited"
]
rg_info[
  CN < 2 & !grepl("X|Y", chr) & inheritance == "denovo" & (MF < 1.5 | MM < 1.5),
  inheritance := "inherited"
]

rg_info[
  inheritance == "denovo" & !grepl("X|Y", chr) & svtype == "DUP",
  mos_dup_thresh := sapply(qnorm(0.98, 0, 0.5 / sqrt(bins)) + 2, \(x) max(x, 2.1))
]
rg_info[
  inheritance == "denovo" & !grepl("X|Y", chr) & svtype == "DEL",
  mos_del_thresh := sapply(qnorm(0.98, 0, 0.5 / sqrt(bins)) * -1 + 2, \(x) min(x, 1.9))
]
rg_info[
  MF > mos_dup_thresh,
  inheritance := "mosaic_father"
]
rg_info[
  MM > mos_dup_thresh,
  inheritance := "mosaic_mother"
]
rg_info[
  MF < mos_del_thresh,
  inheritance := "mosaic_father"
]
rg_info[
  MM < mos_del_thresh,
  inheritance := "mosaic_mother"
]
rg_info[
  MF < mos_del_thresh & MM < mos_del_thresh,
  inheritance := "mosaic_both"
]

# Write annotations to file --------------------------------------------------
out <- ped[, list(sample_id, family_id)][g_c2, on = c(sample_id = "sample")]
setnames(out, "sample_id", "sample")

fwrite(out, output, sep = "\t")
