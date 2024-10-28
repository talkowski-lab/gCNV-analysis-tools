suppressPackageStartupMessages(library(data.table))

matrix_dir <- commandArgs(trailingOnly = TRUE)[[1]]
output <- commandArgs(trailingOnly = TRUE)[[2]]

message(sprintf("searching for matrices in %s", matrix_dir))
mat_paths <- list.files(matrix_dir, pattern = 'ploidy.tsv$', full.names = TRUE)
if (length(mat_paths) == 0) {
    stop("no ploidy matricies found")
}

if (length(mat_paths) == 1) {
    message("found 1 matrix. copying to output")
    file.copy(mat_paths[[1]], output, overwrite = TRUE)
    quit(save = "no")
}

message(sprintf("found %d matrices. merging", length(mat_paths)))
mats <- lapply(mat_paths, \(x) fread(x, header = TRUE, sep = "\t", key = "sample"))
ploidy_mat <- Reduce(\(x, y) merge(x, y, by = "sample", all = TRUE), mats)

message("writing merged matrices to output")
fwrite(ploidy_mat, output, sep = "\t", quote = FALSE)
