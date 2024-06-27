#' Calculate antibiotic-specific MiRIx values
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
#' @name mirix_antibiotic
#' @examples
#' h22 <- weiss2021_data[weiss2021_data$sample_id %in% "Healthy.22",]
#' mirix_vancomycin(h22$proportion, h22$lineage)
NULL


#' @rdname mirix_antibiotic
#' @export
mirix_vancomycin <- function(abundance,
                             lineage,
                             replace_zero = 1e-4,
                             antibiotic_db = whatbacteria::taxon_susceptibility,
                             phenotype_db = whatbacteria::taxon_phenotypes) {
  susceptibility <- antibiotic_susceptibility_vancomycin(
    lineage, antibiotic_db, phenotype_db)
  mirix(abundance, susceptibility, replace_zero)
}

#' @rdname mirix_antibiotic
#' @export
mirix_doxycycline <- function(abundance,
                               lineage,
                               replace_zero = 1e-4,
                               antibiotic_db = whatbacteria::taxon_susceptibility) {
  susceptibility <- antibiotic_susceptibility_tetracycline(
    lineage, antibiotic_db)
  mirix(abundance, susceptibility, replace_zero)
}

#' @rdname mirix_antibiotic
#' @export
mirix_amoxicillin <- function(abundance,
                             lineage,
                             replace_zero = 1e-4,
                             antibiotic_db = whatbacteria::taxon_susceptibility) {
  susceptibility <- antibiotic_susceptibility_penicillin(
    lineage, antibiotic_db)
  mirix(abundance, susceptibility, replace_zero)
}

#' @rdname mirix_antibiotic
#' @export
mirix_metronidazole <- function(abundance,
                            lineage,
                            replace_zero = 1e-4,
                            phenotype_db = whatbacteria::taxon_phenotypes) {
  susceptibility <- phenotype_susceptibility(
    lineage = lineage,
    phenotype = "aerobic_status",
    susceptibility = c(
      "aerobe" = "resistant",
      "facultative anaerobe" = "resistant",
      "obligate anaerobe" = "susceptible"),
    db = phenotype_db)
  mirix(abundance, susceptibility, replace_zero)
}

#' @rdname mirix_antibiotic
#' @export
mirix_ciprofloxacin <- function(abundance,
                          lineage,
                          replace_zero = 1e-4,
                          phenotype_db = whatbacteria::taxon_phenotypes) {
  susceptibility <- phenotype_susceptibility(
    lineage = lineage,
    phenotype = "aerobic_status",
    susceptibility = c(
      "aerobe" = "susceptible",
      "facultative anaerobe" = "susceptible",
      "obligate anaerobe" = "resistant"),
    db = phenotype_db)
  mirix(abundance, susceptibility, replace_zero)
}

#' @rdname mirix_antibiotic
#' @export
mirix_gentamicin <- function(abundance,
                                 lineage,
                                 replace_zero = 1e-4,
                                 antibiotic_db = whatbacteria::taxon_susceptibility,
                                 phenotype_db = whatbacteria::taxon_phenotypes) {
  susceptibility <- antibiotic_susceptibility_aminoglycoside(
    lineage, antibiotic_db, phenotype_db)
  mirix(abundance, susceptibility, replace_zero)
}

#' Calculate Microbiome Response Index (MiRIx) values
#'
#' @param abundance A vector of taxon abundances in a sample
#' @param susceptibility A character vector of taxon susceptibility, with
#'   values that are "susceptible", "resistant", or \code{NA}
#' @param replace_zero Zero-replacement value. If the numerator or denominator
#'   is smaller than this value, they will be replaced with the number here.
#'   This number should generally be slightly less than the minimum measurement
#'   that could be acquired with the method used. For count data, a value of
#'   0.5 is typical. For relative abundances, a number that is slightly lower
#'   than the lowest relative abundance will work.
#'
#' @import whatbacteria
#' @return The MiRIx value
#' @export
mirix <- function (abundance, susceptibility, replace_zero = 1e-4) {
  x_susceptible <- sum(abundance[susceptibility %in% "susceptible"])
  x_resistant <- sum(abundance[susceptibility %in% "resistant"])
  if ((x_resistant < replace_zero) && (x_susceptible < replace_zero)) {
    warning(
      "Numerator and denominator both less than the zero-replacement value.")
  }
  x_susceptible <- max(x_susceptible, replace_zero)
  x_resistant <- max(x_resistant, replace_zero)
  log10(x_susceptible / x_resistant)
}
