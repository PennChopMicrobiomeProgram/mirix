predict_abundances <- function (index_values, abundances, susceptibility) {
  is_resistant <- susceptibility %in% "resistant"
  is_susceptible <- susceptibility %in% "susceptible"
  old_total_resistant <- sum(abundances[is_resistant])
  old_total_susceptible <- sum(abundances[is_susceptible])
  total <- old_total_resistant + old_total_susceptible

  adjust_abundances <- function (val) {
    p_resistant <- logistic(val, base = 10)
    p_susceptible <- 1 - p_resistant
    new_total_resistant <- total * p_resistant
    new_total_susceptible <- total * p_susceptible
    multiplier_resistant <- new_total_resistant / old_total_resistant
    multiplier_susceptible <- new_total_susceptible / old_total_susceptible
    result <- abundances
    result[is_resistant] <- result[is_resistant] * multiplier_resistant
    result[is_susceptible] <- result[is_susceptible] * multiplier_susceptible
    result
  }

  if (length(index_values) > 1) {
    vapply(index_values, adjust_abundances, rep(1.0, length(abundances)))
  } else {
    adjust_abundances(index_values)
  }
}

logistic <- function (x, base = exp(1)) {
  ex <- 10 ^ x
  ex / (ex + 1)
}
