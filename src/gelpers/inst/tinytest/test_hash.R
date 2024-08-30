h <- utils::hashtab("identical", 26)
set_hash(h, LETTERS, 1:26)
expect_equal(utils::numhash(h), 26)

values <- get_hash(h, LETTERS)
expect_equal(values, as.list(1:26))
