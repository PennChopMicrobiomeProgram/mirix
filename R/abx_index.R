#' Calculate the antibiotic index for vancomycin
#'
#' @param abundance A vector of taxon abundances in a sample
#' @param lineage A character vector of taxonomic assignments or lineages
#' @param antibiotic_db A data frame with columns named "taxon", "rank",
#'   "antibiotic", and "value"
#' @param phenotype_db A data frame with columns named "taxon", "rank", and
#'   "gram_stain"
#'
#' @return The vancomycin antibiotic index for the sample
#' @export
#'
#' @examples
#' h22 <- weiss2021_data[weiss2021_data$sample_id %in% "Healthy.22",]
#' vancomycin_index(h22$proportion, h22$lineage)
vancomycin_index <- function(abundance,
                             lineage,
                             antibiotic_db = taxon_susceptibility,
                             phenotype_db = taxon_phenotypes) {
  susceptibility <- vancomycin_susceptibility(
    lineage, antibiotic_db, phenotype_db)
  antibiotic_index(abundance, susceptibility)
}

#' Calculate the antibiotic index for tetracycline
#'
#' @param abundance A vector of taxon abundances in a sample
#' @param lineage A character vector of taxonomic assignments or lineages
#' @param antibiotic_db A data frame with columns named "taxon", "rank",
#'   "antibiotic", and "value"
#'
#' @return The tetracycline antibiotic index for the sample
#' @export
#'
#' @examples
#' h22 <- weiss2021_data[weiss2021_data$sample_id %in% "Healthy.22",]
#' tetracycline_index(h22$proportion, h22$lineage)
tetracycline_index <- function(abundance,
                               lineage,
                               antibiotic_db = taxon_susceptibility) {
  susceptibility <- tetracycline_susceptibility(
    lineage, antibiotic_db)
  antibiotic_index(abundance, susceptibility)
}

#' Calculate the antibiotic index for penicillin-like antibiotics
#'
#' @param abundance A vector of taxon abundances in a sample
#' @param lineage A character vector of taxonomic assignments or lineages
#' @param antibiotic_db A data frame with columns named "taxon", "rank",
#'   "antibiotic", and "value"
#'
#' @return The penicillin antibiotic index for the sample
#' @export
#'
#' @examples
#' h22 <- weiss2021_data[weiss2021_data$sample_id %in% "Healthy.22",]
#' penicillin_index(h22$proportion, h22$lineage)
penicillin_index <- function(abundance,
                             lineage,
                             antibiotic_db = taxon_susceptibility) {
  susceptibility <- penicillin_susceptibility(
    lineage, antibiotic_db)
  antibiotic_index(abundance, susceptibility)
}

#' Calculate the antibiotic index for Gram-positive bacteria
#'
#' Antibiotics such as glycopeptides, macrolides, oxazolidinones, lincosamides,
#' and lipopeptides aside from vancomycin.
#'
#' @param abundance A vector of taxon abundances in a sample
#' @param lineage A character vector of taxonomic assignments or lineages
#' @param phenotype_db A data frame with columns named "taxon", "rank", and
#'   "gram_stain"
#'
#' @return The Gram-positive antibiotic index for the sample
#' @export
#'
#' @examples
#' h22 <- weiss2021_data[weiss2021_data$sample_id %in% "Healthy.22",]
#' gram_positive_index(h22$proportion, h22$lineage)
gram_positive_index <- function(abundance,
                           lineage,
                           phenotype_db = taxon_phenotypes) {
  susceptibility <- phenotype_susceptibility(
    lineage = lineage,
    phenotype = "gram_stain",
    susceptibility = c(
      "Gram-positive" = "susceptible",
      "Gram-negative" = "resistant"),
    db = phenotype_db)
  antibiotic_index(abundance, susceptibility)
}

#' Calculate the antibiotic index for Gram-negative bacteria
#'
#' Antibiotics such as polymyxin and aztreonam.
#'
#' @param abundance A vector of taxon abundances in a sample
#' @param lineage A character vector of taxonomic assignments or lineages
#' @param phenotype_db A data frame with columns named "taxon", "rank", and
#'   "gram_stain"
#'
#' @return The Gram-negative antibiotic index for the sample
#' @export
#'
#' @examples
#' h22 <- weiss2021_data[weiss2021_data$sample_id %in% "Healthy.22",]
#' gram_negative_index(h22$proportion, h22$lineage)
gram_negative_index <- function(abundance,
                           lineage,
                           phenotype_db = taxon_phenotypes) {
  susceptibility <- phenotype_susceptibility(
    lineage = lineage,
    phenotype = "gram_stain",
    susceptibility = c(
      "Gram-positive" = "resistant",
      "Gram-negative" = "susceptible"),
    db = phenotype_db)
  antibiotic_index(abundance, susceptibility)
}

#' Calculate the antibiotic index for obligate anaerboes
#'
#' Antibiotics such as nitroimidazole.
#'
#' @param abundance A vector of taxon abundances in a sample
#' @param lineage A character vector of taxonomic assignments or lineages
#' @param phenotype_db A data frame with columns named "taxon", "rank", and
#'   "aerobic_status"
#'
#' @return The obligate anaerobe antibiotic index for the sample
#' @export
#'
#' @examples
#' h22 <- weiss2021_data[weiss2021_data$sample_id %in% "Healthy.22",]
#' anaerobes_index(h22$proportion, h22$lineage)
anaerobes_index <- function(abundance,
                            lineage,
                            phenotype_db = taxon_phenotypes) {
  susceptibility <- phenotype_susceptibility(
    lineage = lineage,
    phenotype = "aerobic_status",
    susceptibility = c(
      "aerobe" = "resistant",
      "facultative anaerobe" = "resistant",
      "obligate anaerobe" = "susceptible"),
    db = phenotype_db)
  antibiotic_index(abundance, susceptibility)
}

#' Calculate the antibiotic index for aerobes and facultative anaerobes
#'
#' Antibiotics such as fluoroquinolones.
#'
#' @param abundance A vector of taxon abundances in a sample
#' @param lineage A character vector of taxonomic assignments or lineages
#' @param phenotype_db A data frame with columns named "taxon", "rank", and
#'   "aerobic_status"
#'
#' @return The aerobe antibiotic index for the sample
#' @export
#'
#' @examples
#' h22 <- weiss2021_data[weiss2021_data$sample_id %in% "Healthy.22",]
#' aerobes_index(h22$proportion, h22$lineage)
aerobes_index <- function(abundance,
                          lineage,
                          phenotype_db = taxon_phenotypes) {
  susceptibility <- phenotype_susceptibility(
    lineage = lineage,
    phenotype = "aerobic_status",
    susceptibility = c(
      "aerobe" = "susceptible",
      "facultative anaerobe" = "susceptible",
      "obligate anaerobe" = "resistant"),
    db = phenotype_db)
  antibiotic_index(abundance, susceptibility)
}

#' Calculate the antibiotic index for aminoglycoside antibiotics
#'
#' @param abundance A vector of taxon abundances in a sample
#' @param lineage A character vector of taxonomic assignments or lineages
#' @param antibiotic_db A data frame with columns named "taxon", "rank",
#'   "antibiotic", and "value"
#' @param phenotype_db A data frame with columns named "taxon", "rank", and
#'   "gram_stain"
#'
#' @return The aminoglycoside antibiotic index for the sample
#' @export
#'
#' @examples
#' h22 <- weiss2021_data[weiss2021_data$sample_id %in% "Healthy.22",]
#' aminoglycoside_index(h22$proportion, h22$lineage)
aminoglycoside_index <- function(abundance,
                                 lineage,
                                 antibiotic_db = taxon_susceptibility,
                                 phenotype_db = taxon_phenotypes) {
  susceptibility <- aminoglycoside_susceptibility(
    lineage, antibiotic_db, phenotype_db)
  antibiotic_index(abundance, susceptibility)
}

antibiotic_index <- function (abundance, susceptibility) {
  x_resistant <- sum(abundance[susceptibility %in% "resistant"])
  x_susceptible <- sum(abundance[susceptibility %in% "susceptible"])
  log10(x_resistant / x_susceptible)
}
