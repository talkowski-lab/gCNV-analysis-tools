map <- read_dcr_list("data/dCR_paths/valid_one_column_dcr_paths.tsv")
expect_true(is.hashtab(map))
expect_equal(
    utils::gethash(map, "batch00"),
    c("/foo/bar/batch00_CASE.dcr.bed.gz", "/foo/bar/batch00_COHORT.dcr.bed.gz")
)
expect_equal(
    utils::gethash(map, "batch01"),
    "/foo/bar/batch01_COHORT.dcr.bed.gz"
)

map <- read_dcr_list("data/dCR_paths/valid_two_column_dcr_paths.tsv")
expect_true(is.hashtab(map))
expect_equal(
    utils::gethash(map, "cluster00"),
    c("/foo/bar/batch00_CASE.dcr.bed.gz", "/foo/bar/batch00_COHORT.dcr.bed.gz")
)
expect_equal(
    utils::gethash(map, "cluster01"),
    "/foo/bar/batch01_COHORT.dcr.bed.gz"
)

expect_error(
    read_dcr_list("data/dCR_paths/invalid_three_column_dcr_paths.tsv")
)

expect_error(
    read_dcr_list("data/dCR_paths/invalid_one_column_dcr_paths.tsv")
)

expect_error(
    read_dcr_list("data/dCR_paths/invalid_two_column_dcr_paths.tsv")
)
