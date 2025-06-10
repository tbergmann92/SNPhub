test_that("snp_types works", {
    expect_length(allowed_vals_single, 11)
    expect_length(allowed_vals_double, 17)
    expect_length(iupac_map, 17)
})