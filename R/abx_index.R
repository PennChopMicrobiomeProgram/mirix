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

#' Function to calculate antibiotics index targeting Gram-positive bacteria
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

#' Function to calculate antibiotics index targeting Gram-negative bacteria
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

#' Function to calculate antibiotics index targeting anaerobes such as Nitroimidazole
#'
#' @param abundance A vector of relative abundances of bacterial taxons for a single sample
#' @param lineage Name of taxonomy lineage for each relative abundance in a sample (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia etc.)
#'
#' @return The calculated antibiotics index for the sample
#' @export
#'
#' @examples
#' apply(abx_test_df, 2, anaerobes_index, row.names(abx_test_df))
#'
anaerobes_index <- function(abundance, lineage) {
  idx <- "anaerobe"
  suscept_vector <- is_susceptible(lineage, idx)
  calc_index(abundance, suscept_vector)
}

#' Function to calculate antibiotics index targeting aerobes such as Fluoroquinolone
#'
#' @param abundance A vector of relative abundances of bacterial taxons for a single sample
#' @param lineage Name of taxonomy lineage for each relative abundance in a sample (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia etc.)
#'
#' @return The calculated antibiotics index for the sample
#' @export
#'
#' @examples
#' apply(abx_test_df, 2, aerobes_index, row.names(abx_test_df))
#'
aerobes_index <- function(abundance, lineage) {
  idx <- "aerobe"
  suscept_vector <- is_susceptible(lineage, idx)
  calc_index(abundance, suscept_vector)
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




#' Calculate the susceptiblity vector from taxonomic lineage based on the given bacterial phenotype
#'
#' @param taxa The lineage name of taxonomic ranks
#' @param idx The bacterial phenotype to calculate the susceptibility
#'
#' @return A vector of 0 and 1, where 0 is resistant and 1 is susceptible
#'
is_susceptible <- function(taxa, idx) {
  suscept_vector <- rep(0, length(taxa))
  if(any(grepl("tetracycline|penicillin", idx))) { ##assume all taxa is susceptible for tetracyclines and penicillins
    suscept_vector <- rep(1, length(taxa))
  }

  for(each_idx in idx) {

    abx_df <- abx_idx_df[abx_idx_df$attribute == each_idx, ]

    for(lvl in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {

      taxon_df <- abx_df[abx_df$rank == lvl, ]

      if (dim(taxon_df)[1] != 0) {
        pos_suscept_vector <- rep(0, length(taxa))
        neg_suscept_vector <- rep(0, length(taxa))
        pos_rank <- taxon_df[taxon_df$boo, ]
        neg_rank <- taxon_df[!taxon_df$boo, ]

        if ((dim(pos_rank)[1] != 0)) {
          pos_suscept_vector <- grepl(paste0(pos_rank$name, collapse = "|"), taxa)
          pos_suscept_vector <- pos_suscept_vector*1
        }

        if ((dim(neg_rank)[1] != 0)) {
          neg_suscept_vector <- grepl(paste0(neg_rank$name, collapse = "|"), taxa)
          neg_suscept_vector <- neg_suscept_vector*-1
        }

        new_vector <- pos_suscept_vector + neg_suscept_vector
        #print(new_vector)
        zero_vector_idx <- which(new_vector == 0)
        old_suscept_vector <- suscept_vector[zero_vector_idx]
        suscept_vector <- new_vector
        suscept_vector[zero_vector_idx] <- old_suscept_vector
      }
    }
  }
  suscept_vector[which(suscept_vector == -1)] <- 0
  suscept_vector
}

antibiotic_index <- function (abundance, susceptibility) {
  x_resistant <- sum(abundance[susceptibility %in% "resistant"])
  x_susceptible <- sum(abundance[susceptibility %in% "susceptible"])
  log10(x_resistant / x_susceptible)
}

#' Calculate the index based on the susceptible vector
#'
#' @param sample_vector The vector of taxonomic abundances for a sample
#' @param suscept_vector The vector of how susceptible a species is to the specified antibiotics
#'
#' @return The calculated antibiotics index for a sample
#'
calc_index <- function(sample_vector, suscept_vector) {

  TF_suscept_vector <- suscept_vector > 0.5
  sum_suscept_taxa <- sum(sample_vector[TF_suscept_vector])

  log10((1-sum_suscept_taxa)/sum_suscept_taxa)
}

#' Function to plot antibiotics indices for each sample
#'
#' @param abx_vector The antibiotics vector generated from any of the abxidx functions
#' @param order Order the samples from highest to lowest antibiotics index before plotting
#'
#' @return A plot for the antibiotics index for each sample
#' @export
#'
#' @examples
#' vanco_idx <- apply(abx_test_df, 2, vancomycin_index, row.names(abx_test_df))
#' abx_idx_plot(vanco_idx)
#'
abx_idx_plot <- function(abx_vector, order = F) {
  ##infinite values are replaced to 10
  show_name <- NULL

  abx_vector[abx_vector > 10] <- 10
  abx_vector[abx_vector < -10] <- -10

  if(order) {
    abx_vector <- sort(abx_vector, decreasing = TRUE)
  }

  vector_cols <- ifelse(abx_vector > 0, "green", "red")

  if(length(abx_vector) > 50) {
    show_name <- FALSE
  }

  barplot(abx_vector, ylab = "Antibiotics index", names.arg = show_name, col = vector_cols, border = NA, space = 0, las = 2)
}
