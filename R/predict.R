#' Predict taxon abundances at given values of an index
#'
#' @param index_value Value or values of the index at which to make predictions.
#' @param abundances A vector of taxon abundances in a sample.
#' @param susceptibility A character vector of antibiotic susceptibility, with
#'   values that are "susceptible", "resistant", or \code{NA}.
#' @return A new vector of abundances if \code{index_value} has length 1. If
#'   \code{index_value} is a vector longer than 1, this function returns a
#'   matrix, where each column contains the predicted abundances.
#' @export
predict_abundance <- function (index_value, abundance, susceptibility) {
  is_resistant <- susceptibility %in% "resistant"
  is_susceptible <- susceptibility %in% "susceptible"
  old_total_resistant <- sum(abundance[is_resistant])
  old_total_susceptible <- sum(abundance[is_susceptible])
  total <- old_total_resistant + old_total_susceptible

  adjust_abundance <- function (val) {
    p_resistant <- logistic(val, base = 10)
    p_susceptible <- 1 - p_resistant
    new_total_resistant <- total * p_resistant
    new_total_susceptible <- total * p_susceptible
    multiplier_resistant <- new_total_resistant / old_total_resistant
    multiplier_susceptible <- new_total_susceptible / old_total_susceptible
    result <- abundance
    result[is_resistant] <- result[is_resistant] * multiplier_resistant
    result[is_susceptible] <- result[is_susceptible] * multiplier_susceptible
    result
  }

  if (length(index_value) > 1) {
    vapply(index_value, adjust_abundance, rep(1.0, length(abundance)))
  } else {
    adjust_abundance(index_value)
  }
}

logistic <- function (x, base = exp(1)) {
  ex <- 10 ^ x
  ex / (ex + 1)
}
