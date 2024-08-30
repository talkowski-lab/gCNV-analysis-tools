gbins <- GBins(c("chr19", rep("chr20", 10)),
               IRanges::IRanges(c(77, seq.int(1, 19, 2)),
                                c(100, seq.int(2, 20, 2))))
region <- gregion("chr20", 5, 14)
xreg <- expand_region_by_bins(region, gbins, pad = 0.2)
expect_identical(xreg, gregion("chr20", 3, 16))

# Test expansion doesn't overflow into other sequences
region <- gregion("chr20", 1, 10)
xreg <- expand_region_by_bins(region, gbins, pad = 0.2)
expect_identical(xreg, gregion("chr20", 1, 12))
