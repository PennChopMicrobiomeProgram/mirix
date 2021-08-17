#' Calculate antibiotic-specific index values
#'
#' @param abundance A vector of taxon abundances in a sample
#' @param lineage A character vector of taxonomic assignments or lineages
#' @param replace_zero Zero-replacement value. If the numerator or denominator
#'   is smaller than this value, they will be replaced with the number here.
#'   This number should generally be slightly less than the minimum measurement
#'   that could be acquired with the method used. For count data, a value of
#'   0.5 is typical. For relative abundances, a number that is slightly lower
#'   than the lowest relative abundance will work.
#' @param antibiotic_db A data frame with columns named "taxon", "rank",
#'   "antibiotic", and "value"
#' @param phenotype_db A data frame with columns named "taxon", "rank",
#'   "gram_stain", and "aerobic_status"
#' @name antibiotic_specific_index
#' @examples
#' h22 <- weiss2021_data[weiss2021_data$sample_id %in% "Healthy.22",]
#' vancomycin_index(h22$proportion, h22$lineage)
NULL


#' @rdname antibiotic_specific_index
#' @export
vancomycin_index <- function(abundance,
                             lineage,
                             replace_zero = 1e-4,
                             antibiotic_db = taxon_susceptibility,
                             phenotype_db = taxon_phenotypes) {
  susceptibility <- vancomycin_susceptibility(
    lineage, antibiotic_db, phenotype_db)
  antibiotic_index(abundance, susceptibility, replace_zero)
}

#' @rdname antibiotic_specific_index
#' @export
tetracycline_index <- function(abundance,
                               lineage,
                               replace_zero = 1e-4,
                               antibiotic_db = taxon_susceptibility) {
  susceptibility <- tetracycline_susceptibility(
    lineage, antibiotic_db)
  antibiotic_index(abundance, susceptibility, replace_zero)
}

#' @rdname antibiotic_specific_index
#' @export
penicillin_index <- function(abundance,
                             lineage,
                             replace_zero = 1e-4,
                             antibiotic_db = taxon_susceptibility) {
  susceptibility <- penicillin_susceptibility(
    lineage, antibiotic_db)
  antibiotic_index(abundance, susceptibility, replace_zero)
}

#' @rdname antibiotic_specific_index
#' @export
gram_positive_index <- function(abundance,
                           lineage,
                           replace_zero = 1e-4,
                           phenotype_db = taxon_phenotypes) {
  susceptibility <- phenotype_susceptibility(
    lineage = lineage,
    phenotype = "gram_stain",
    susceptibility = c(
      "Gram-positive" = "susceptible",
      "Gram-negative" = "resistant"),
    db = phenotype_db)
  antibiotic_index(abundance, susceptibility, replace_zero)
}

#' @rdname antibiotic_specific_index
#' @export
gram_negative_index <- function(abundance,
                           lineage,
                           replace_zero = 1e-4,
                           phenotype_db = taxon_phenotypes) {
  susceptibility <- phenotype_susceptibility(
    lineage = lineage,
    phenotype = "gram_stain",
    susceptibility = c(
      "Gram-positive" = "resistant",
      "Gram-negative" = "susceptible"),
    db = phenotype_db)
  antibiotic_index(abundance, susceptibility, replace_zero)
}

#' @rdname antibiotic_specific_index
#' @export
anaerobes_index <- function(abundance,
                            lineage,
                            replace_zero = 1e-4,
                            phenotype_db = taxon_phenotypes) {
  susceptibility <- phenotype_susceptibility(
    lineage = lineage,
    phenotype = "aerobic_status",
    susceptibility = c(
      "aerobe" = "resistant",
      "facultative anaerobe" = "resistant",
      "obligate anaerobe" = "susceptible"),
    db = phenotype_db)
  antibiotic_index(abundance, susceptibility, replace_zero)
}

#' @rdname antibiotic_specific_index
#' @export
aerobes_index <- function(abundance,
                          lineage,
                          replace_zero = 1e-4,
                          phenotype_db = taxon_phenotypes) {
  susceptibility <- phenotype_susceptibility(
    lineage = lineage,
    phenotype = "aerobic_status",
    susceptibility = c(
      "aerobe" = "susceptible",
      "facultative anaerobe" = "susceptible",
      "obligate anaerobe" = "resistant"),
    db = phenotype_db)
  antibiotic_index(abundance, susceptibility, replace_zero)
}

#' @rdname antibiotic_specific_index
#' @export
aminoglycoside_index <- function(abundance,
                                 lineage,
                                 replace_zero = 1e-4,
                                 antibiotic_db = taxon_susceptibility,
                                 phenotype_db = taxon_phenotypes) {
  susceptibility <- aminoglycoside_susceptibility(
    lineage, antibiotic_db, phenotype_db)
  antibiotic_index(abundance, susceptibility, replace_zero)
}

#' Calculate the antibiotic index
#'
#' @param abundance A vector of taxon abundances in a sample
#' @param susceptibility A character vector of antibiotic susceptibility, with
#'   values that are "susceptible", "resistant", or \code{NA}
#' @param replace_zero Zero-replacement value. If the numerator or denominator
#'   is smaller than this value, they will be replaced with the number here.
#'   This number should generally be slightly less than the minimum measurement
#'   that could be acquired with the method used. For count data, a value of
#'   0.5 is typical. For relative abundances, a number that is slightly lower
#'   than the lowest relative abundance will work.
#'
#' @return The antibiotic index value
#' @export
antibiotic_index <- function (abundance, susceptibility, replace_zero = 1e-4) {
  x_resistant <- sum(abundance[susceptibility %in% "resistant"])
  x_susceptible <- sum(abundance[susceptibility %in% "susceptible"])
  if ((x_resistant < replace_zero) && (x_susceptible < replace_zero)) {
    warning(
      "Numerator and denominator both less than the zero-replacement value.")
  }
  x_resistant <- max(x_resistant, replace_zero)
  x_susceptible <- max(x_susceptible, replace_zero)
  log10(x_resistant / x_susceptible)
}
