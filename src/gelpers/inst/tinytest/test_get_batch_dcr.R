# test reading a dCR file with squeezing
region <- gregion("chr5", 35697592, 138201546)
dcr <- get_batch_dcr(
    region,
    c("data/dCR/batch00_part0_dcr.bed.gz", "data/dCR/batch00_part1_dcr.bed.gz")
)
expect_equal(ncol(dcr), 9)
expect_equal(nrow(dcr), 3)
expect_equal(colnames(dcr)[[1]], "chr")
expect_equal(colnames(dcr)[[2]], "start")
expect_equal(colnames(dcr)[[3]], "end")
expect_true(all(dcr[4:9] >= 0))
expect_true(all(dcr[4:9] <= 5))

# test reading a dCR file without squeezing
region <- gregion("chrY", 17514734, 21521525)
dcr <- get_batch_dcr(
    region,
    c("data/dCR/batch00_part0_dcr.bed.gz", "data/dCR/batch00_part1_dcr.bed.gz"),
    squeeze = FALSE
)
expect_equal(ncol(dcr), 9)
expect_equal(nrow(dcr), 3)
expect_equal(colnames(dcr)[[1]], "chr")
expect_equal(colnames(dcr)[[2]], "start")
expect_equal(colnames(dcr)[[3]], "end")
expect_true(any(dcr[4:9] < 0))

region <- gregion("chr17", 11869036, 20814826)
dcr <- get_batch_dcr(
    region,
    c("data/dCR/batch00_part0_dcr.bed.gz", "data/dCR/batch00_part1_dcr.bed.gz"),
    squeeze = FALSE
)
expect_equal(ncol(dcr), 9)
expect_equal(nrow(dcr), 3)
expect_equal(colnames(dcr)[[1]], "chr")
expect_equal(colnames(dcr)[[2]], "start")
expect_equal(colnames(dcr)[[3]], "end")
expect_true(any(dcr[4:9] > 5))

# test reading dCR without header
region <- gregion("chr1", 70249450, 111389364)
expect_error(
    get_batch_dcr(region, "data/dCR/no_header_dcr.bed.gz"),
    "missing header",
    class = "dcr_parse_error"
)

# test getting non-existing sequence
region <- gregion("chrZ", 1e5, 1e8)
expect_error(
    get_batch_dcr(region, "data/dCR/batch00_part0_dcr.bed.gz"),
    "no sequence",
    class = "dcr_parse_error"
)

# test reading a dCR with out-of-order region columns
region <- gregion("chr1", 70249450, 111389364)
expect_error(
    get_batch_dcr(region, "data/dCR/out_of_order_columns_dcr.bed.gz"),
    "first three columns",
    class = "dcr_parse_error"
)

# test reading a batch with unequal number of rows
region <- gregion("chr1", 20713356, 111389364)
expect_error(
    get_batch_dcr(
        region,
        c("data/dCR/mismatched_rows_part0_dcr.bed.gz", "data/dCR/mismatched_rows_part1_dcr.bed.gz"),
    )
)

# test reading a dCR with non-syntatic names
region <- gregion("chr1", 20713356, 111389364)
dcr <- get_batch_dcr(region, "data/dCR/non_syntatic_samples_dcr.bed.gz")
expect_true(all(c("123sample", "+sample", "sample1777") %in% colnames(dcr)))
expect_equal(ncol(dcr), 6)
expect_equal(nrow(dcr), 4)
