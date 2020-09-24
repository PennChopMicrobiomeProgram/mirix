testthat::context("Testing abxindx functionality")

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
  gram_pos_idx_outcome <- apply(abx_test_df, 2, gram_pos_index, row.names(abx_test_df))

  testthat::expect_type(gram_pos_idx_outcome, "double")
  testthat::expect_equal(length(gram_pos_idx_outcome), 4)

  ##gram negative test
  gram_neg_idx_outcome <- apply(abx_test_df, 2, gram_neg_index, row.names(abx_test_df))

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

testthat::context("Testing list of susceptible or resistant bacteria")

testthat::test_that("Testing list of susceptible or resistant bacteria", {

  ##vanco
  vanco_list <- do.call(rbind, apply(abx_test_df, 2, vancomycin_list, row.names(abx_test_df)))
  testthat::expect_equal(paste0(colnames(vanco_list), collapse = ", "), "lineage, abundance, phenotype")
  testthat::expect_equal(ncol(vanco_list), 3)
  testthat::expect_equal(class(vanco_list), "data.frame")

  ##tet
  tet_list <- do.call(rbind, apply(abx_test_df, 2, tetracycline_list, row.names(abx_test_df)))
  testthat::expect_equal(paste0(colnames(tet_list), collapse = ", "), "lineage, abundance, phenotype")
  testthat::expect_equal(ncol(tet_list), 3)
  testthat::expect_equal(class(tet_list), "data.frame")

  ##gram_pos
  gp_list <- do.call(rbind, apply(abx_test_df, 2, gram_pos_list, row.names(abx_test_df)))
  testthat::expect_equal(paste0(colnames(gp_list), collapse = ", "), "lineage, abundance, phenotype")
  testthat::expect_equal(ncol(gp_list), 3)
  testthat::expect_equal(class(gp_list), "data.frame")

  ##gram_neg
  gn_list <- do.call(rbind, apply(abx_test_df, 2, gram_neg_list, row.names(abx_test_df)))
  testthat::expect_equal(paste0(colnames(gn_list), collapse = ", "), "lineage, abundance, phenotype")
  testthat::expect_equal(ncol(gn_list), 3)
  testthat::expect_equal(class(gn_list), "data.frame")

  ##anaerobe
  anaero_list <- do.call(rbind, apply(abx_test_df, 2, anaerobe_list, row.names(abx_test_df)))
  testthat::expect_equal(paste0(colnames(anaero_list), collapse = ", "), "lineage, abundance, phenotype")
  testthat::expect_equal(ncol(anaero_list), 3)
  testthat::expect_equal(class(anaero_list), "data.frame")

  ##aerobe
  aero_list <- do.call(rbind, apply(abx_test_df, 2, aerobe_list, row.names(abx_test_df)))
  testthat::expect_equal(paste0(colnames(aero_list), collapse = ", "), "lineage, abundance, phenotype")
  testthat::expect_equal(ncol(aero_list), 3)
  testthat::expect_equal(class(aero_list), "data.frame")

})


testthat::context("Testing abx_idx_df")

testthat::test_that("Testing abx_idx_df contains the required columns", {

  testthat::expect_equal(paste0(colnames(abx_idx_df), collapse = ", "), "attribute, boo, name, rank, doi")
  testthat::expect_equal(ncol(abx_idx_df), 5)
  testthat::expect_equal(class(abx_idx_df), "data.frame")

})

testthat::context("Testing abx_idx_plot")

testthat::test_that("Testing abx_idx_plot to plot antibiotics indices", {

  testthat::expect_type(abx_idx_plot(apply(abx_test_df, 2, vancomycin_index, row.names(abx_test_df))), "double")
  testthat::expect_type(abx_idx_plot(apply(abx_test_df, 2, vancomycin_index, row.names(abx_test_df)), order = TRUE), "double")

})
