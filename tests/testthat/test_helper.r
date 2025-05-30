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