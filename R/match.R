# Determine the annotation values for each lineage
#
# The database is split into ranks and the values are determined for
# each rank by matching taxa from the rank-specific database to the vector of
# lineages. If a lineage matches to multiple taxa of different ranks, the
# value of the taxon with the lowest rank is selected.
#
# The taxonomic ranks, in order from highest to lowest, are Kingdom, Phylum,
# Class, Order, Family, Genus, and Species. The ranks in the database must be
# capitalized, exactly as they are above.
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
      paste(lineages[multi_matches], collapse = "\n"))
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
