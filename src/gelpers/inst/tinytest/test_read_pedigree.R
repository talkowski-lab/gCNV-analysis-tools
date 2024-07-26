# test reading a valid pedigree
ped_clean <- utils::read.table(
    "data/pedigree/pedigree_clean.tsv",
    header = TRUE,
    colClasses = c(rep("character", 4), rep("integer", 2))
)

ped <- read_pedigree("data/pedigree/pedigree_valid.tsv")
expect_identical(ped, ped_clean)

# test reading pedigree with missing samples
expect_warning(
    read_pedigree("data/pedigree/pedigree_missing_samples.tsv"),
    pattern = "removing 1 missing sample(s) from pedigree",
    strict = TRUE
)
ped <- suppressWarnings(read_pedigree("data/pedigree/pedigree_missing_samples.tsv"))
expect_identical(ped, ped_clean)

# test reading pedigree with invalid phenotypes
ped <- read_pedigree("data/pedigree/pedigree_invalid_phenotypes.tsv")
expect_identical(ped, ped_clean)

# test reading pedigree with invalid sex
ped <- read_pedigree("data/pedigree/pedigree_invalid_sexes.tsv")
expect_identical(ped, ped_clean)

# test that missing family, paternal, and maternal IDs get converted to NA
ped <- read_pedigree("data/pedigree/pedigree_missing_other_ids.tsv")
expect_identical(ped, ped_clean)

# test that duplicate sample IDs are removed
expect_warning(
    read_pedigree("data/pedigree/pedigree_duplicate_samples.tsv"),
    pattern = "removing 2 duplicate sample(s) from pedigree",
    strict = TRUE
)
ped <- suppressWarnings(read_pedigree("data/pedigree/pedigree_duplicate_samples.tsv"))
expect_identical(ped, ped_clean)
