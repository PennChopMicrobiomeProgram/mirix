weiss_sepsis <- subset(weiss2021_data, sample_id %in% "Sepsis.22441.C")
weiss_sepsis$proportion <- round(weiss_sepsis$proportion, 3)

weiss_healthy <- subset(weiss2021_data, sample_id %in% "Healthy.148")
weiss_healthy$proportion <- round(weiss_healthy$proportion, 3)

test_that("vancomycin_susceptibility works on Weiss examples", {
  expect_equal(
    vancomycin_susceptibility(weiss_sepsis$lineage),
    c("resistant", "susceptible", "resistant", "resistant", "resistant"))

  expect_equal(
    vancomycin_susceptibility(weiss_healthy$lineage),
    c("resistant", "susceptible", "resistant", "resistant", "susceptible",
      "susceptible", "resistant", "susceptible", "susceptible", "susceptible",
      "susceptible", "susceptible", "susceptible", "resistant", "susceptible",
      "susceptible", "resistant", "susceptible", "resistant", "susceptible",
      "susceptible", "resistant", "susceptible", "susceptible", "susceptible",
      "resistant"))

})

test_that("vancomycin_index works on Weiss examples", {
  expect_equal(
    vancomycin_index(weiss_sepsis$proportion, weiss_sepsis$lineage),
    0.5670429, tolerance = 1e-5)

  expect_equal(
    vancomycin_index(weiss_healthy$proportion, weiss_healthy$lineage),
    0.3609116, tolerance = 1e-5)
})

test_that("tetracycline_susceptibility works on Weiss examples", {
  expect_equal(
    tetracycline_susceptibility(weiss_sepsis$lineage),
    c("susceptible", "susceptible", "resistant", "susceptible", "susceptible"))

  expect_equal(
    tetracycline_susceptibility(weiss_healthy$lineage),
    c("susceptible", "susceptible", "resistant", "susceptible", "susceptible",
      "susceptible", "susceptible", "susceptible", "susceptible", "susceptible",
      "susceptible", "susceptible", "susceptible", "susceptible", "susceptible",
      "susceptible", "susceptible", "susceptible", "susceptible", "susceptible",
      "susceptible", "susceptible", "susceptible", "susceptible", "susceptible",
      "susceptible"))

})

test_that("tetracycline_index works on Weiss examples", {
  expect_equal(
    tetracycline_index(weiss_sepsis$proportion, weiss_sepsis$lineage),
    -1.879934, tolerance = 1e-5)

  expect_equal(
    tetracycline_index(weiss_healthy$proportion, weiss_healthy$lineage),
    -1.045228, tolerance = 1e-5)
})

test_that("penicillin_susceptibility works on Weiss examples", {
  expect_equal(
    penicillin_susceptibility(weiss_sepsis$lineage),
    c("susceptible", "susceptible", "resistant", "susceptible", "susceptible"))

  expect_equal(
    penicillin_susceptibility(weiss_healthy$lineage),
    c("susceptible", "susceptible", "resistant", "susceptible", "susceptible",
      "susceptible", "susceptible", "susceptible", "susceptible", "susceptible",
      "susceptible", "susceptible", "susceptible", "susceptible", "susceptible",
      "susceptible", "susceptible", "susceptible", "susceptible", "susceptible",
      "susceptible", "susceptible", "susceptible", "susceptible", "susceptible",
      "susceptible"))
})

test_that("penicillin_index works on Weiss examples", {
  expect_equal(
    penicillin_index(weiss_sepsis$proportion, weiss_sepsis$lineage),
    -1.879934, tolerance = 1e-5)

  expect_equal(
    penicillin_index(weiss_healthy$proportion, weiss_healthy$lineage),
    -1.045228, tolerance = 1e-5)
})

test_that("aminoglycoside_susceptibility works on Weiss examples", {
  expect_equal(
    aminoglycoside_susceptibility(weiss_sepsis$lineage),
    c("resistant", "susceptible", "resistant", "resistant", "susceptible"))

  expect_equal(
    aminoglycoside_susceptibility(weiss_healthy$lineage),
    c("resistant", "resistant", "resistant", "resistant", "resistant",
      "resistant", "resistant", "resistant", "resistant", "resistant",
      "resistant", "resistant", "susceptible", "resistant", "resistant",
      "resistant", "resistant", "resistant", "resistant", "resistant",
      NA, NA, "resistant", "resistant", "resistant", "resistant"))
})

test_that("aminoglycoside_index works on Weiss examples", {
  expect_equal(
    aminoglycoside_index(weiss_sepsis$proportion, weiss_sepsis$lineage),
    0.5644559, tolerance = 1e-5)

  expect_equal(
    aminoglycoside_index(weiss_healthy$proportion, weiss_healthy$lineage),
    2.21396, tolerance = 1e-5)
})

test_that("gram_positive_index works on Weiss examples", {
  expect_equal(
    gram_positive_index(weiss_sepsis$proportion, weiss_sepsis$lineage),
    0.5644559, tolerance = 1e-5)

  expect_equal(
    gram_positive_index(weiss_healthy$proportion, weiss_healthy$lineage),
    0.9650513, tolerance = 1e-5)
})

test_that("gram_positive_index and gram_negative_index are inverses", {
  expect_equal(
    gram_positive_index(weiss_sepsis$proportion, weiss_sepsis$lineage),
    -gram_negative_index(weiss_sepsis$proportion, weiss_sepsis$lineage))

  expect_equal(
    gram_positive_index(weiss_healthy$proportion, weiss_healthy$lineage),
    -gram_negative_index(weiss_healthy$proportion, weiss_healthy$lineage))
})

test_that("anaerobes_index works on Weiss examples", {
  expect_equal(
    anaerobes_index(weiss_sepsis$proportion, weiss_sepsis$lineage),
    1.788434, tolerance = 1e-5)

  expect_equal(
    anaerobes_index(weiss_healthy$proportion, weiss_healthy$lineage),
    -2.21396, tolerance = 1e-5)
})

test_that("anaerobes_index and aerobes_index are inverses", {
  expect_equal(
    anaerobes_index(weiss_sepsis$proportion, weiss_sepsis$lineage),
    -aerobes_index(weiss_sepsis$proportion, weiss_sepsis$lineage))

  expect_equal(
    anaerobes_index(weiss_healthy$proportion, weiss_healthy$lineage),
    -aerobes_index(weiss_healthy$proportion, weiss_healthy$lineage))
})

test_that("is_susceptible returns correct values for taxa", {
  expect_equal(
    is_susceptible(
      "k__Bacteria; p__Firmicutes; c__Negativicutes; o__Veillonellales; f__Veillonellaceae; g__Veillonella",
      c("gram_positive", "vancomycin")),
    0.0)

  expect_equal(
    is_susceptible(
      "k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Clostridiaceae; g__Clostridium; s__Clostridium perfringens",
      c("gram_positive", "vancomycin")),
    1.0)
})

testthat::test_that("Testing abxidx on abx_test_df", {

  ##vanco test
  vanco_idx_outcome <- apply(abx_test_df, 2, vancomycin_index, row.names(abx_test_df))

  testthat::expect_type(vanco_idx_outcome, "double")
  testthat::expect_equal(length(vanco_idx_outcome), 4)
  testthat::expect_equal(vanco_idx_outcome > 0, c(a = FALSE, b = TRUE, c = FALSE, d = TRUE)) ##test that there are two samples have positive idx indices

  ##tetracycline test
  tet_idx_outcome <- apply(abx_test_df, 2, tetracycline_index, row.names(abx_test_df))

  testthat::expect_type(tet_idx_outcome, "double")
  testthat::expect_equal(length(tet_idx_outcome), 4)

  ##gram positive test
  gram_pos_idx_outcome <- apply(abx_test_df, 2, gram_positive_index, row.names(abx_test_df))

  testthat::expect_type(gram_pos_idx_outcome, "double")
  testthat::expect_equal(length(gram_pos_idx_outcome), 4)

  ##gram negative test
  gram_neg_idx_outcome <- apply(abx_test_df, 2, gram_negative_index, row.names(abx_test_df))

  testthat::expect_type(gram_neg_idx_outcome, "double")
  testthat::expect_equal(length(gram_neg_idx_outcome), 4)

  ##anaerobes test
  ana_idx_outcome <- apply(abx_test_df, 2, anaerobes_index, row.names(abx_test_df))

  testthat::expect_type(ana_idx_outcome, "double")
  testthat::expect_equal(length(ana_idx_outcome), 4)

  ##aerobes test
  aero_idx_outcome <- apply(abx_test_df, 2, aerobes_index, row.names(abx_test_df))

  testthat::expect_type(aero_idx_outcome, "double")
  testthat::expect_equal(length(aero_idx_outcome), 4)

})

testthat::test_that("Testing abx_idx_df contains the required columns", {

  testthat::expect_equal(paste0(colnames(abx_idx_df), collapse = ", "), "attribute, boo, name, rank, doi")
  testthat::expect_equal(ncol(abx_idx_df), 5)
  testthat::expect_equal(class(abx_idx_df), "data.frame")

})

testthat::test_that("Testing abx_idx_plot to plot antibiotics indices", {

  testthat::expect_type(abx_idx_plot(apply(abx_test_df, 2, vancomycin_index, row.names(abx_test_df))), "double")
  testthat::expect_type(abx_idx_plot(apply(abx_test_df, 2, vancomycin_index, row.names(abx_test_df)), order = TRUE), "double")

})
