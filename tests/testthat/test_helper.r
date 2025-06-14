test_that("calc_percent() works", {
  expect_equal(calc_percent(50, 100), 50.0)
  expect_equal(calc_percent(100, 100), 100)
  expect_equal(calc_percent(0, 100), 0)
  expect_length(calc_percent(0, 100), 1)
})

test_that("check_unrecognized() works", {
    expect_error(check_unrecognized(c("A", "T", "G"), c("A", "T")))
    expect_error(check_unrecognized(NA, c("A", "T")))
    expect_error(check_unrecognized(c(NA, "A"), c("A", "T")))
})

test_that("chars() works", {
    expect_vector(chars("ATG"))
    expect_vector(chars(c("ATG")))
    expect_identical(chars("ATG"), c("A", "T", "G"))
    expect_identical(chars(c("ATG", NA)), c("A", "T", "G"))
    expect_type(chars("ATG"), "character")
})

test_that("assert() works", {
    expect_no_error(assert(2 == 2))
    expect_no_error(assert(NULL))
    expect_no_error(assert(c(1 == 1, 2 == 2)))
    expect_error(assert(c(1 == 1, 1 == 2)))
    expect_error(assert())
    expect_error(assert(NA))
    expect_error(assert(1 == 2))
    expect_error(assert(1 == 2, "error"))
    expect_error(assert(c(2 == 2, 2 == 3), "error"), "error")
})

test_that("remap_snps() works", {
    mat <- matrix(c("AA", "AG", "--", "GG"), ncol = 2)
    mat2 <- matrix(c("AA", "AG", "--", "FF"), ncol = 2)
    expect_equal(remap_snps(mat, iupac_map),  matrix(c("A", "R", "-", "G"), ncol = 2))
    expect_equal(remap_snps(mat2, iupac_map),  matrix(c("A", "R", "-", "FF"), ncol = 2))
})