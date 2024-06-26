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
                                       db = mirixdb::taxon_susceptibility) {
  what_antibiotic(lineage, antibiotic, db)
}

what_antibiotic <- function (lineage,
                             antibiotic,
                             db = mirixdb::taxon_susceptibility) {
  is_relevant <- db$antibiotic %in% antibiotic
  db <- db[is_relevant, c("taxon", "rank", "value")]

  susceptibility_values <- match_annotation(lineage, db)
  susceptibility_values
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
                                      db = mirixdb::taxon_phenotypes) {
  phenotype_values <- what_phenotype(lineage, phenotype, db)
  susceptibility_values <- susceptibility[phenotype_values]
  susceptibility_values <- unname(susceptibility_values)
  susceptibility_values
}

what_phenotype <- function (lineage,
                            phenotype,
                            db = mirixdb::taxon_phenotypes) {
  db <- db[, c("taxon", "rank", phenotype)]
  # match_annotation() requires a column named "value"
  colnames(db)[3] <- "value"
  match_annotation(lineage, db)
}

# Determine the annotation values for each lineage
#
# @param lineage A vector of taxonomic assignments or lineages
# @param db A data frame with columns named "taxon", "rank", and "value"
# @return A vector of assigned values
match_annotation <- function (lineage, db) {
  get_rank_specific_db <- function (r) {
    rank_is_r <- db[["rank"]] %in% r
    db[rank_is_r,]
  }
  db_ranks <- lapply(rev(taxonomic_ranks), get_rank_specific_db)
  names(db_ranks) <- rev(taxonomic_ranks)

  get_values_by_rank <- function (rank_specific_db) {
    taxa_idx <- match_taxa(lineage, rank_specific_db[["taxon"]])
    rank_specific_db[["value"]][taxa_idx]
  }
  values_by_rank <- vapply(
    db_ranks,
    get_values_by_rank,
    rep("a", length(lineage)))

  if (length(lineage) == 1) {
    assigned_values <- first_non_na_value(values_by_rank)
  } else {
    assigned_values <- apply(values_by_rank, 1, first_non_na_value)
  }
  assigned_values
}

# The 'official' taxonomic ranks supported by this package
taxonomic_ranks <- c(
  "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Return the first value that is not NA. If all values are NA, return NA. The
# resultant vector will not have names.
first_non_na_value <- function (x) {
  unname(x[first_true_idx(!is.na(x))])
}

# For each lineage, return the index of the taxon that is found within the
# lineage. If no taxa are found, return NA for that element. If multiple taxa
# are found, we issue a warning and return the index of the first taxon in the
# vector of taxa.
match_taxa <- function (lineages, taxa) {
  n_lineages <- length(lineages)
  if (length(taxa) == 0) {
    return(rep_len(NA_character_, length(lineages)))
  }

  taxa_patterns <- paste0("(?<=__|\\b)(?:", taxa, ")\\b")
  lineage_matches <- vapply(
    X = taxa_patterns,
    FUN = grepl,
    FUN.VALUE = rep_len(TRUE, n_lineages),
    x = lineages,
    perl = TRUE,
    USE.NAMES = TRUE)

  # If the user passes only one lineage, lineage_matches will be a vector
  # rather than an array. After some trial and error, I found that it's better
  # to deal with this at each stage of the computation, rather than trying to
  # coerce the vector to an array up front.
  if (n_lineages == 1) {
    multi_matches <- sum(lineage_matches) > 1
  } else {
    multi_matches <- rowSums(lineage_matches) > 1
  }
  if (any(multi_matches)) {
    warning(
      "The following lineages match more than one taxon:\n",
      paste(lineages[multi_matches], collapse = "\n"), "\n")
  }

  if (n_lineages == 1) {
    taxon_idx <- first_true_idx(lineage_matches)
  } else {
    taxon_idx <- apply(lineage_matches, 1, first_true_idx)
  }
  taxon_idx
}

# Return the first index of a boolean vector that is TRUE. If all elements of
# the vector are FALSE, return NA. Tempted to call this function minwhich.
first_true_idx <- function (x) {
  if (any(x)) {
    min(which(x == TRUE))
  } else {
    NA_integer_
  }
}
