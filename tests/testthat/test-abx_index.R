testthat::context("Testing abx_indx functionality")

testthat::test_that("Testing abx_idx on test df", {

  idx_outcome <- abx_index(test_df, "Vancomycin")

  testthat::expect_type(idx_outcome, "list")
  testthat::expect_equal(length(idx_outcome), 4)
  testthat::expect_equal(idx_outcome > 0, c(a = TRUE, b = FALSE, c = TRUE, d = FALSE)) ##test that there are two samples have positive idx indices

  ##write tests for the df as well

})
