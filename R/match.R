#' Evaluate susceptibility based on antibiotic-specific annotations
#'
#' This is the main engine for determining antibiotic susceptibility based on
#' a database of antibiotic-specific annotations for various taxa.
#'
#' @param lineage A character vector of taxonomic assignments or lineages
#' @param antibiotic The name of the antibiotic or antibiotic class in \code{db}
#' @param db A data frame with columns named "taxon", "rank", "antibiotic",
#'   and "value"
#' @return A vector of assigned susceptibility values, which should be either
#'   "susceptible", "resistant", or \code{NA}
#' @details
#' To determine susceptibility, the database is first filtered to include only
#' rows relevant to the antibiotic of interest. Then, the filtered database is
#' split into ranks. The susceptibility values are determined for each rank by
#' matching taxa from the rank-specific database to the vector of lineages.
#' If a lineage matches to multiple taxa of different ranks, the value of the
#' taxon with the lowest rank is selected.
#'
#' The taxonomic ranks, in order from highest to lowest, are Kingdom, Phylum,
#' Class, Order, Family, Genus, and Species. The ranks in the database must be
#' capitalized, exactly as they are written here.
#' @examples
#' antibiotic_susceptibility(
#'   c("Enterococcus faecalis", "Lactobacillus", "Lactobacillus delbrueckii"),
#'   "vancomycin")
#' @name antibiotic_susceptibility
#' @export
antibiotic_susceptibility <- function (lineage,
                                       antibiotic,
                                       db = whatbacteria::taxon_susceptibility) {
  whatbacteria::what_antibiotic(lineage, antibiotic, db)
}

#' Evaluate antibiotic susceptibility based on phenotype
#'
#' This is the main engine for determining antibiotic susceptibility based on
#' a database of phenotypes for various taxa. To run this function, you'll need
#' the database as well as a vector that specifies the susceptibility for each
#' value of the phenotype.
#'
#' @param lineage A character vector of taxonomic assignments or lineages
#' @param phenotype The name of the column in \code{db} that contains the
#'   phenotype of interest
#' @param susceptibility A character vector specifying the susceptibility for
#'   different values of the phenotype. The names should correspond to elements
#'   in the column specified with \code{phenotype}. The values should be either
#'   "susceptible" or "resistant"
#' @param db A data frame with columns named "taxon", "rank", and the column
#'   name specified in \code{phenotype}
#' @return A vector of assigned susceptibility values, which should be either
#'   "susceptible", "resistant", or \code{NA}
#' @details
#' This function operates much like \code{antibiotic_susceptibility}, except
#' that it pulls phenotype values from the database instead of susceptibility
#' information. To subsequently determine susceptibility from phenotype, this
#' function uses the named vector provided in the \code{susceptibility}
#' argument.
#'
#' As a reminder, the taxonomic ranks, in order from highest to lowest, are
#' Kingdom, Phylum, Class, Order, Family, Genus, and Species. The ranks in the
#' database must be capitalized, exactly as they are written here.
#' @examples
#' phenotype_susceptibility(
#'   c("Bacteroidetes", "Firmicutes", "Firmicutes; Negativicutes"),
#'   "gram_stain",
#'   c("Gram-positive" = "susceptible", "Gram-negative" = "resistant"))
#' @name phenotype_susceptibility
#' @export
phenotype_susceptibility <- function (lineage,
                                      phenotype,
                                      susceptibility,
                                      db = whatbacteria::taxon_phenotypes) {
  phenotype_values <- whatbacteria::what_phenotype(lineage, phenotype, db)
  susceptibility_values <- susceptibility[phenotype_values]
  susceptibility_values <- unname(susceptibility_values)
  susceptibility_values
}
