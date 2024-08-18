dcr_path_map <- utils::hashtab(type = "identical", 2)
utils::sethash(dcr_path_map, "batch00", "data/trio_dcr/trio_dcr_00.dcr.bed.gz")
utils::sethash(dcr_path_map, "batch01", "data/trio_dcr/trio_dcr_01.dcr.bed.gz")

# Test getting trio from a single batch
region <- gregion("chr2", 6e7, 211e6)
dcr <- get_trio_dcr(region,
                    paste0("sample", c("00", "01", "02")),
                    rep("batch00", 3),
                    dcr_path_map)
expect_equal(c("chr", "start", "end", "sample00", "sample01", "sample02"), colnames(dcr))
expect_equal(4, nrow(dcr))
expect_true(all(dcr$chr == "chr2"))
expect_equal(dcr$start, c(64100300L, 84604241L, 177550899L, 201778343L))
expect_equal(dcr$end, c(64100597L, 84604649L, 177551846L, 201778585L))
expect_false(anyNA(dcr, recursive = TRUE))


# Test getting trio from two batches keeping all ranges
region <- gregion("chr19", 2e6, 39e6)
dcr <- get_trio_dcr(region,
                    paste0("sample", c("00", "04", "05")),
                    paste0("batch", c("00", "01", "01")),
                    dcr_path_map)
expect_equal(c("chr", "start", "end", "sample00", "sample04", "sample05"), colnames(dcr))
expect_equal(5, nrow(dcr))
expect_true(all(dcr$chr == "chr19"))
expect_true(is.na(dcr[dcr$start == 38528251, "sample04", drop = TRUE]))
expect_true(is.na(dcr[dcr$start == 38528251, "sample05", drop = TRUE]))

# Test getting trio from two batches keeping only common ranges
dcr <- get_trio_dcr(region,
                    paste0("sample", c("00", "04", "05")),
                    paste0("batch", c("00", "01", "01")),
                    dcr_path_map,
                    keep_all_ranges = FALSE)
expect_equal(c("chr", "start", "end", "sample00", "sample04", "sample05"), colnames(dcr))
expect_equal(4, nrow(dcr))
expect_true(all(dcr$chr == "chr19"))
expect_false(anyNA(dcr, recursive = TRUE))

# Test getting trio from two batches with background
dcr <- get_trio_dcr(region,
                    paste0("sample", c("00", "04", "05")),
                    paste0("batch", c("00", "01", "01")),
                    dcr_path_map,
                    include_bg = TRUE)
expect_equal(c("chr", "start", "end", sprintf("sample%02d", c(0, 4, 5, 1, 2))), colnames(dcr))
expect_equal(5, nrow(dcr))
expect_true(anyNA(dcr, recursive = TRUE))

# Test getting trio without any ranges in common
region <- gregion("chr19", 38528251L, 38528576L)
dcr <- get_trio_dcr(region,
                    paste0("sample", c("00", "04", "05")),
                    paste0("batch", c("00", "01", "01")),
                    dcr_path_map,
                    include_bg = TRUE,
                    keep_all_ranges = FALSE)
expect_equal(8, ncol(dcr))
expect_equal(0, nrow(dcr))
