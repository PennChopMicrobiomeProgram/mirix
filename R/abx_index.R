#' Function to calculate antibiotics index for Vancomycin
#'
#' @param abundance A list of relative abundances of bacterial taxons for a single sample
#' @param lineage Name of taxonomy lineage for each relative abundance in a sample (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia etc.)
#'
#' @return The calculated antibiotics index for the sample
#' @export
#'
#' @examples
#' apply(abx_test_df, 2, vancomycin_index, row.names(abx_test_df))
#'
vancomycin_index <- function(abundance, lineage) {
  idx <- c("gram_positive", "vancomycin")
  suscept_vector <- is_susceptible(lineage, idx)
  calc_index(abundance, suscept_vector)
}

#' Function to calculate antibiotics index for Tetracycline
#'
#' @param abundance A list of relative abundances of bacterial taxons for a single sample
#' @param lineage Name of taxonomy lineage for each relative abundance in a sample (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia etc.)
#'
#' @return The calculated antibiotics index for the sample
#' @export
#'
#' @examples
#' apply(abx_test_df, 2, tetracycline_index, row.names(abx_test_df))
#'
tetracycline_index <- function(abundance, lineage) {
  idx <- "tetracycline"
  suscept_vector <- is_susceptible(lineage, idx)
  calc_index(abundance, suscept_vector)
}

#' Function to calculate antibiotics index targeting gram positive bacteria such as Glycopeptides, Macrolides, Oxazolidinones, Lincosamides, Lipopeptides, and Amoxicillin aside from Vancomycin (see \code{vancomycin_index})
#'
#' @param abundance A list of relative abundances of bacterial taxons for a single sample
#' @param lineage Name of taxonomy lineage for each relative abundance in a sample (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia etc.)
#'
#' @return The calculated antibiotics index for the sample
#' @export
#'
#' @examples
#' apply(abx_test_df, 2, gram_pos_index, row.names(abx_test_df))
#'
gram_pos_index <- function(abundance, lineage) {
  idx <- "gram_positive"
  suscept_vector <- is_susceptible(lineage, idx)
  calc_index(abundance, suscept_vector)
}

#' Function to calculate antibiotics index targeting anaerobes such as Polymyxin and Aztreonam
#'
#' @param abundance A list of relative abundances of bacterial taxons for a single sample
#' @param lineage Name of taxonomy lineage for each relative abundance in a sample (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia etc.)
#'
#' @return The calculated antibiotics index for the sample
#' @export
#'
#' @examples
#' apply(abx_test_df, 2, gram_neg_index, row.names(abx_test_df))
#'
gram_neg_index <- function(abundance, lineage) {
  idx <- "gram_negative"
  suscept_vector <- is_susceptible(lineage, idx)
  calc_index(abundance, suscept_vector)
}

#' Function to calculate antibiotics index targeting anaerobes such as Nitroimidazole
#'
#' @param abundance A list of relative abundances of bacterial taxons for a single sample
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

#' Function to calculate antibiotics index targeting anaerobes such as Fluoroquinolone
#'
#' @param abundance A list of relative abundances of bacterial taxons for a single sample
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

#' Calculate the susceptiblity vector from taxonomic lineage based on the given bacterial phenotype
#'
#' @param taxa The lineage name of taxonomic levels
#' @param idx The bacterial phenotype to calculate the susceptibility
#'
#' @return A vector of 0 and 1, where 0 is resistant and 1 is susceptible
#'
is_susceptible <- function(taxa, idx) {
  suscept_vector <- rep(0, length(taxa))

  for(each_idx in idx) {
    abx_df <- filter(grep_df, grepl(each_idx, attribute))

    for(lvl in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
      taxon_df <- abx_df %>%
        filter(grepl(lvl, level))

      if (dim(taxon_df)[1] != 0) {
        pos_suscept_vector <- rep(0, length(taxa))
        neg_suscept_vector <- rep(0, length(taxa))
        pos_level <- taxon_df %>%
          filter(boo)
        neg_level <- taxon_df %>%
          filter(!boo)

        if ((dim(pos_level)[1] != 0)) {
          pos_suscept_vector <- grepl(paste0(pos_level$name, collapse = "|"), taxa)
          pos_suscept_vector <- pos_suscept_vector*1
        }

        if ((dim(neg_level)[1] != 0)) {
          neg_suscept_vector <- grepl(paste0(neg_level$name, collapse = "|"), taxa)
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

  log10(sum_suscept_taxa/(1-sum_suscept_taxa))
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
  show_name <- TRUE

  plotting_vector <- gsub(-Inf, -10, abx_vector)
  plotting_vector <- as.numeric(gsub(Inf, 10, plotting_vector))

  if(order) {
    plotting_vector <- sort(plotting_vector, decreasing = TRUE)
  }

  vector_cols <- ifelse(plotting_vector > 0, "green", "red")

  if(length(plotting_vector) > 50) {
    show_name <- FALSE
  }

  barplot(unlist(plotting_vector), xlab = "Samples", ylab = "Antibiotics index", names.arg = show_name, col = vector_cols, border = NA, space = 0)
}
