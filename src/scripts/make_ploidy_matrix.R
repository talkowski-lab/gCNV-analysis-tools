# Ploidy matrix generator

usage <- function() {
'Usage: Rscript make_ploidy_matrix.R CALLSET DCRS CONTIGS OUTPUT

CALLSET  Final callset produced by the gCNV pipeline
DCRS     A file with either the paths, one per line, to the dCR matrices or
         two tab-separated columns with the batch ID in the first column and
         path to the dCR matrix for that batch in the second column
CONTIGS  A file listing the names of the contigs to consider
OUTPUT   Where to write the output
'
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
    cat(usage(), sep = "", file = stderr())
    quit(save = "no", 1)
}
args <- as.list(args)
names(args) <- c("CALLSET", "DCRS", "CONTIGS", "OUTPUT")

suppressPackageStartupMessages(library(gelpers))

TABIX_MAX_SEQLEN <- 536870912L

contig_ploidy <- function(contig, dcr_paths, samples) {
    region <- gregion(contig, 1, TABIX_MAX_SEQLEN)
    dcr <- tryCatch(
        as.data.table(get_samples_dcr(region, dcr_paths, samples)),
        dcr_parse_error = function(cnd) NULL
    )

    if (is.null(dcr)) {
        mean_contig <- data.table(sample = samples, contig = NA_real_)
        setnames(mean_contig, "contig", contig)
    } else {
        mean_contig <- dcr[, lapply(.SD, mean), .SDcols = !c("chr", "start", "end")]
        mean_contig <- melt(mean_contig,
                            measure.vars = colnames(mean_contig),
                            variable.name = "sample",
                            value.name = contig,
                            variable.factor = FALSE)
    }

    setkey(mean_contig, sample)

    mean_contig
}

batch_ploidy <- function(batch, samples, dcr_map, contigs) {
    dcr_paths <- gethash(dcr_map, batch)
    if (is.null(dcr_paths)) {
        stop(sprintf("batch %s does not have an associated dCR path", batch),
             call. = FALSE)
    }
    mean_contigs <- lapply(contigs,
                           \(x) contig_ploidy(x, dcr_paths, samples))
    ploidy <- Reduce(\(x, y) merge(x, y, "sample", all= TRUE),
                     mean_contigs)

    ploidy
}

make_ploidy_matrix <- function(callset, dcrs_map, contigs) {
    x <- unique(callset[, list(sample, batch)])

    dcr_groups <- split(x$sample, x$batch)
    batches <- mapply(\(x, y) batch_ploidy(x, y, dcrs_map, contigs),
                      names(dcr_groups),
                      dcr_groups,
                      SIMPLIFY = FALSE)

    rbindlist(batches)
}

log_info("reading callset")
calls <- read_callset(args$CALLSET)

log_info("reading dCR paths")
dcrs <- read_dcr_list(args$DCRS)

log_info("reading contigs list")
contigs <- unique(readLines(args$CONTIGS))
if (length(contigs) == 0) {
    log_error("no contigs given")
    quit("no", 1)
}

log_info("making ploidy matrix")
ploidy <- make_ploidy_matrix(calls, dcrs, contigs)

log_info("writing output")
fwrite(ploidy, file = args$OUTPUT, sep = "\t", na = "NA")
