#' Determine susceptibility to vancomycin
#'
#' @param lineage A character vector of taxonomic assignments or lineages
#' @param antibiotic_db A data frame with columns named "taxon", "rank",
#'   "antibiotic", and "value"
#' @param phenotype_db A data frame with columns named "taxon", "rank", and
#'   "gram_stain"
#' @return A vector of assigned susceptibility values, which should be either
#'   "susceptible", "resistant", or \code{NA}
#' @name antibiotic_specific_susceptibility
NULL

#' @rdname antibiotic_specific_susceptibility
#' @export
antibiotic_susceptibility_vancomycin <- function (lineage,
                                       antibiotic_db = taxon_susceptibility,
                                       phenotype_db = taxon_phenotypes) {
  # Gram-positive organisms are susceptible to vancomycin
  ph_sus <- phenotype_susceptibility(
    lineage = lineage,
    phenotype = "gram_stain",
    susceptibility = c(
      "Gram-positive" = "susceptible",
      "Gram-negative" = "resistant"),
    db = phenotype_db)
  abx_sus <- antibiotic_susceptibility(
    lineage = lineage,
    antibiotic = "vancomycin",
    db = antibiotic_db)
  ifelse(is.na(abx_sus), ph_sus, abx_sus)
}

#' @rdname antibiotic_specific_susceptibility
#' @export
antibiotic_susceptibility_oxacillin <- function (lineage,
                                                 antibiotic_db = taxon_susceptibility,
                                                 phenotype_db = taxon_phenotypes) {
  ph_sus <- phenotype_susceptibility(
    lineage = lineage,
    phenotype = "gram_stain",
    susceptibility = c(
      "Gram-positive" = "susceptible",
      "Gram-negative" = "resistant"),
    db = phenotype_db)
  abx_sus <- antibiotic_susceptibility(
    lineage = lineage,
    antibiotic = "oxacillin",
    db = antibiotic_db)
  ifelse(is.na(abx_sus), ph_sus, abx_sus)
}

#' @rdname antibiotic_specific_susceptibility
#' @export
antibiotic_susceptibility_tetracycline <- function (lineage,
                                         antibiotic_db = taxon_susceptibility) {
  intrinsic_sus <- rep("susceptible", length(lineage))
  abx_sus <- antibiotic_susceptibility(
    lineage = lineage,
    antibiotic = "tetracycline",
    db = antibiotic_db)
  ifelse(is.na(abx_sus), intrinsic_sus, abx_sus)
}

#' @rdname antibiotic_specific_susceptibility
#' @export
antibiotic_susceptibility_penicillin <- function(lineage,
                                      antibiotic_db = taxon_susceptibility) {
  intrinsic_sus <- rep("susceptible", length(lineage))
  abx_sus <- antibiotic_susceptibility(
    lineage = lineage,
    antibiotic = "penicillin",
    db = antibiotic_db)
  ifelse(is.na(abx_sus), intrinsic_sus, abx_sus)
}

#' @rdname antibiotic_specific_susceptibility
#' @export
antibiotic_susceptibility_aminoglycoside <- function (lineage,
                                           antibiotic_db = taxon_susceptibility,
                                           phenotype_db = taxon_phenotypes) {
  gram_stain_db <- taxon_phenotypes[,c("taxon", "rank", "gram_stain")]
  colnames(gram_stain_db)[3] <- "value"
  gram_stain_phenotype <- match_annotation(lineage, gram_stain_db)

  aerobic_status_db <- taxon_phenotypes[,c("taxon", "rank", "aerobic_status")]
  colnames(aerobic_status_db)[3] <- "value"
  aerobic_status_phenotype <- match_annotation(lineage, aerobic_status_db)

  combined_phenotype <- ifelse(
    is.na(gram_stain_phenotype) | is.na(aerobic_status_phenotype),
    NA_character_,
    paste(gram_stain_phenotype, aerobic_status_phenotype))
  susceptibility <- c(
    "Gram-negative aerobe" = "susceptible",
    "Gram-negative facultative anaerobe" = "susceptible",
    "Gram-negative obligate anaerobe" = "resistant",
    "Gram-positive aerobe" = "resistant",
    "Gram-positive facultative anaerobe" = "resistant",
    "Gram-positive obligate anaerobe" = "resistant")
  ph_sus <- susceptibility[combined_phenotype]

  abx_sus <- antibiotic_susceptibility(
    lineage = lineage,
    antibiotic = "aminoglycoside",
    db = antibiotic_db)

  ifelse(is.na(abx_sus), ph_sus, abx_sus)
}
