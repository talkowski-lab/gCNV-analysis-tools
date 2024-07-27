# test getting a single sample
region <- gregion("chr11", 9808721, 93742239)
dcr <- get_samples_dcr(
    region,
    c("data/dCR/batch00_part0_dcr.bed.gz", "data/dCR/batch00_part1_dcr.bed.gz"),
    "sample00"
)
expect_identical(colnames(dcr), c("chr", "start", "end", "sample00"))
expect_equal(nrow(dcr), 4)

# test getting a single sample with background
dcr <- get_samples_dcr(
    region,
    c("data/dCR/batch00_part0_dcr.bed.gz", "data/dCR/batch00_part1_dcr.bed.gz"),
    "sample00",
    include_bg = TRUE
)
expect_equal(ncol(dcr), 9)
expect_equal(nrow(dcr), 4)
expect_equal(colnames(dcr)[1:3], c("chr", "start", "end"))
expect_equal(colnames(dcr)[[4]], "sample00")
expect_true(all(sprintf("sample%02d", 1:5) %in% colnames(dcr)))

# test getting multiple samples
dcr <- get_samples_dcr(
    region,
    c("data/dCR/batch00_part0_dcr.bed.gz", "data/dCR/batch00_part1_dcr.bed.gz"),
    c("sample00", "sample01", "sample04")
)
expect_equal(ncol(dcr), 6)
expect_equal(nrow(dcr), 4)
expect_equal(colnames(dcr)[1:3], c("chr", "start", "end"))
expect_true(
    all(c("sample00", "sample01", "sample04") %in% colnames(dcr))
)

# test getting multiple samples with background
dcr <- get_samples_dcr(
    region,
    c("data/dCR/batch00_part0_dcr.bed.gz", "data/dCR/batch00_part1_dcr.bed.gz"),
    c("sample00", "sample01", "sample04"),
    include_bg = TRUE
)
expect_equal(ncol(dcr), 9)
expect_equal(nrow(dcr), 4)
expect_equal(colnames(dcr)[1:3], c("chr", "start", "end"))
expect_true(
    all(c("sample00", "sample01", "sample04") %in% colnames(dcr)[4:6])
)

# test getting a sample not in the dCR
expect_error(
    get_samples_dcr(
        region,
        c("data/dCR/batch00_part0_dcr.bed.gz", "data/dCR/batch00_part1_dcr.bed.gz"),
        "notasample"
    )
)
