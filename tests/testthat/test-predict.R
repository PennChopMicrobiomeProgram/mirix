test_that("predict_abundance works for simple case", {
  a <- c(0.2, 0.8)
  r <- c("resistant", "susceptible")
  expect_equal(predict_abundance( 0.0, a, r), c( 1 /  2,  1 /  2))
  expect_equal(predict_abundance( 1.0, a, r), c(10 / 11,  1 / 11))
  expect_equal(predict_abundance(-1.0, a, r), c( 1 / 11, 10 / 11))
})

abundance <- c(0.1, 0.2, 0.1, 0.2, 0.3, 0.1)
susceptibility <- c(
  "resistant", "resistant", "resistant",
  "susceptible", "susceptible", "susceptible")
idx <- antibiotic_index(abundance, susceptibility)
abundance_idx0 <- c(0.125, 0.25, 0.125, 0.1666666667, 0.25, 0.08333333333)

test_that("predict_abundance works for normal case", {
  observed <- predict_abundance(0.0, abundance, susceptibility)
  expect_equal(sum(observed[1:3]), 0.5)
  expect_equal(sum(observed[4:6]), 0.5)
  expect_equal(observed, abundance_idx0, tolerance = 1e-5)
})

test_that("predict_abundance works for multiple index values", {
  observed <- predict_abundance(c(idx, 0.0), abundance, susceptibility)
  expected <- matrix(c(abundance, abundance_idx0), ncol = 2)
  expect_equal(observed, expected, tolerance = 1e-5)
})

test_that("predict_abundance only changes resistant/susceptible taxa", {
  observed <- predict_abundance(
    0.0, c(abundance, 0.3), c(susceptibility, NA_character_))
  expected <- c(abundance_idx0, 0.3)
  expect_equal(observed, expected, tolerance = 1e-5)
})

test_that("predict_abundance preserves total abundance", {
  observed <- predict_abundance(0.0, abundance * 5, susceptibility)
  expected <- abundance_idx0 * 5
  expect_equal(observed, expected, tolerance = 1e-5)
})
