# Test creating a gregion object
reg <- gregion("chr21", 400, 700)
expect_inherits(reg, "gregion")
expect_equal(reg$chr, "chr21")
expect_equal(reg$start, 400)
expect_equal(reg$end, 700)

# Test conversion to string
reg_str <- to_string(reg)
expect_equal(reg_str, "chr21:400-700")

# Test printing
expect_stdout(print(reg), "chr21:400-700")

# Test creating invalid gregion objects
expect_error(gregion(0, 1, 1))
expect_error(gregion("chr21", 1:10, 1))
expect_error(gregion("chr21", 1, 1:10))
expect_error(gregion("chr21", -1, 10))
expect_error(gregion("chr21", 1, -10))
expect_error(gregion("chr21", 2, 1))
