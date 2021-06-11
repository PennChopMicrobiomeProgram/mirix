#' Function to calculate antibiotics index for Vancomycin
#'
#' @param abundance A vector of relative abundances of bacterial taxons for a single sample
#' @param lineage Name of taxonomy lineage for each relative abundance in a sample (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia etc.)
#'
#' @return The calculated antibiotics index for the sample
#' @export
#'
#' @examples
#' apply(abx_test_df, 2, vancomycin_index, row.names(abx_test_df))
#'
vancomycin_index <- function(abundance,
                             lineage,
                             antibiotic_db = taxon_susceptibility,
                             phenotype_db = taxon_phenotypes) {
  susceptibility <- vancomycin_susceptibility(
    lineage, antibiotic_db, phenotype_db)
  antibiotic_index(abundance, susceptibility)
}

#' Function to return taxon susceptible or resistant to Vancomycin
#'
#' @param abundance A vector of relative abundances of bacterial taxons for a single sample
#' @param lineage Name of taxonomy lineage for each relative abundance in a sample (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia etc.)
#'
#' @return A dataframe of taxa abundances and whether they are susceptible or resistant to vancomycin
#' @export
#'
#' @examples
#' do.call(rbind, apply(abx_test_df, 2, vancomycin_list, row.names(abx_test_df)))
#'
vancomycin_list <- function(abundance, lineage) {
  idx <- c("gram_positive", "vancomycin")
  suscept_vector <- is_susceptible(lineage, idx)

  sorted_abundance_suscept <- sort(abundance[suscept_vector > 0.5], index.return=TRUE, decreasing = TRUE)
  susceptibles <- data.frame(lineage = character(), abundance = numeric(), phenotype = character())
  if(length(sorted_abundance_suscept$x) != 0) {
    susceptibles <- data.frame(lineage = lineage[suscept_vector > 0.5][sorted_abundance_suscept$ix],
                               abundance = sorted_abundance_suscept$x,
                               phenotype = "susceptible")
  }

  sorted_abundance_resist <- sort(abundance[suscept_vector < 0.5], index.return=TRUE, decreasing = TRUE)
  resistances <- data.frame(lineage = character(), abundance = numeric(), phenotype = character())
  if(length(sorted_abundance_resist$x) != 0) {
    resistances <- data.frame(lineage = lineage[suscept_vector < 0.5][sorted_abundance_resist$ix],
                              abundance = sorted_abundance_resist$x,
                              phenotype = "resistant")
  }
  rbind(susceptibles, resistances, make.row.names = FALSE)
}

#' Function to calculate antibiotics index for Tetracycline
#'
#' @param abundance A vector of relative abundances of bacterial taxons for a single sample
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

#' Function to return taxon susceptible or resistant to Tetracycline
#'
#' @param abundance A vector of relative abundances of bacterial taxons for a single sample
#' @param lineage Name of taxonomy lineage for each relative abundance in a sample (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia etc.)
#'
#' @return A dataframe of taxa abundances and whether they are susceptible or resistant to tetracycline
#' @export
#'
#' @examples
#' do.call(rbind, apply(abx_test_df, 2, tetracycline_list, row.names(abx_test_df)))
#'
tetracycline_list <- function(abundance, lineage) {
  idx <- "tetracycline"
  suscept_vector <- is_susceptible(lineage, idx)

  sorted_abundance_suscept <- sort(abundance[suscept_vector > 0.5], index.return=TRUE, decreasing = TRUE)
  susceptibles <- data.frame(lineage = character(), abundance = numeric(), phenotype = character())
  if(length(sorted_abundance_suscept$x) != 0) {
    susceptibles <- data.frame(lineage = lineage[suscept_vector > 0.5][sorted_abundance_suscept$ix],
                               abundance = sorted_abundance_suscept$x,
                               phenotype = "susceptible")
  }

  sorted_abundance_resist <- sort(abundance[suscept_vector < 0.5], index.return=TRUE, decreasing = TRUE)
  resistances <- data.frame(lineage = character(), abundance = numeric(), phenotype = character())
  if(length(sorted_abundance_resist$x) != 0) {
    resistances <- data.frame(lineage = lineage[suscept_vector < 0.5][sorted_abundance_resist$ix],
                              abundance = sorted_abundance_resist$x,
                              phenotype = "resistant")
  }
  rbind(susceptibles, resistances, make.row.names = FALSE)
}

#' Function to calculate antibiotics index for Penicillin-like antibiotics
#'
#' @param abundance A vector of relative abundances of bacterial taxons for a single sample
#' @param lineage Name of taxonomy lineage for each relative abundance in a sample (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia etc.)
#'
#' @return The calculated antibiotics index for the sample
#' @export
#'
#' @examples
#' apply(abx_test_df, 2, penicillin_index, row.names(abx_test_df))
#'
penicillin_index <- function(abundance, lineage) {
  idx <- "penicillin"
  suscept_vector <- is_susceptible(lineage, idx)
  calc_index(abundance, suscept_vector)
}

#' Function to return taxon susceptible or resistant to Penicillin-like antibiotics
#'
#' @param abundance A vector of relative abundances of bacterial taxons for a single sample
#' @param lineage Name of taxonomy lineage for each relative abundance in a sample (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia etc.)
#'
#' @return A dataframe of taxa abundances and whether they are susceptible or resistant to penicillin-like antibiotics
#' @export
#'
#' @examples
#' do.call(rbind, apply(abx_test_df, 2, penicillin_list, row.names(abx_test_df)))
#'
penicillin_list <- function(abundance, lineage) {
  idx <- "penicillin"
  suscept_vector <- is_susceptible(lineage, idx)

  sorted_abundance_suscept <- sort(abundance[suscept_vector > 0.5], index.return=TRUE, decreasing = TRUE)
  susceptibles <- data.frame(lineage = character(), abundance = numeric(), phenotype = character())
  if(length(sorted_abundance_suscept$x) != 0) {
    susceptibles <- data.frame(lineage = lineage[suscept_vector > 0.5][sorted_abundance_suscept$ix],
                               abundance = sorted_abundance_suscept$x,
                               phenotype = "susceptible")
  }

  sorted_abundance_resist <- sort(abundance[suscept_vector < 0.5], index.return=TRUE, decreasing = TRUE)
  resistances <- data.frame(lineage = character(), abundance = numeric(), phenotype = character())
  if(length(sorted_abundance_resist$x) != 0) {
    resistances <- data.frame(lineage = lineage[suscept_vector < 0.5][sorted_abundance_resist$ix],
                              abundance = sorted_abundance_resist$x,
                              phenotype = "resistant")
  }
  rbind(susceptibles, resistances, make.row.names = FALSE)
}

#' Function to calculate antibiotics index targeting gram positive bacteria such as Glycopeptides, Macrolides, Oxazolidinones, Lincosamides, and Lipopeptides aside from Vancomycin (see \code{vancomycin_index})
#'
#' @param abundance A vector of relative abundances of bacterial taxons for a single sample
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

#' Function to return susceptibility or resistance to antibiotics targeting gram positive bacteria
#'
#' @param abundance A vector of relative abundances of bacterial taxons for a single sample
#' @param lineage Name of taxonomy lineage for each relative abundance in a sample (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia etc.)
#'
#' @return A dataframe of taxa abundances and whether they are susceptible or resistant to anti-gram positive antibiotics
#' @export
#'
#' @examples
#' do.call(rbind, apply(abx_test_df, 2, gram_pos_list, row.names(abx_test_df)))
#'
gram_pos_list <- function(abundance, lineage) {
  idx <- "gram_positive"
  suscept_vector <- is_susceptible(lineage, idx)

  sorted_abundance_suscept <- sort(abundance[suscept_vector > 0.5], index.return=TRUE, decreasing = TRUE)
  susceptibles <- data.frame(lineage = character(), abundance = numeric(), phenotype = character())
  if(length(sorted_abundance_suscept$x) != 0) {
    susceptibles <- data.frame(lineage = lineage[suscept_vector > 0.5][sorted_abundance_suscept$ix],
                               abundance = sorted_abundance_suscept$x,
                               phenotype = "susceptible")
  }

  sorted_abundance_resist <- sort(abundance[suscept_vector < 0.5], index.return=TRUE, decreasing = TRUE)
  resistances <- data.frame(lineage = character(), abundance = numeric(), phenotype = character())
  if(length(sorted_abundance_resist$x) != 0) {
    resistances <- data.frame(lineage = lineage[suscept_vector < 0.5][sorted_abundance_resist$ix],
                              abundance = sorted_abundance_resist$x,
                              phenotype = "resistant")
  }
  rbind(susceptibles, resistances, make.row.names = FALSE)
}

#' Function to calculate antibiotics index targeting gram negatives such as Polymyxin and Aztreonam
#'
#' @param abundance A vector of relative abundances of bacterial taxons for a single sample
#' @param lineage Name of taxonomy lineage for each relative abundance in a sample (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia etc.)
#'
#' @return The calculated antibiotics index for the sample
#' @export
#'
#' @examples
#' apply(abx_test_df, 2, gram_neg_index, row.names(abx_test_df))
#'
gram_neg_index <- function(abundance, lineage) {
  idx <- "gram_positive"
  suscept_vector <- is_susceptible(lineage, idx)
  calc_index(abundance, !suscept_vector)
}

#' Function to return susceptibility or resistance to antibiotics targeting gram negative bacteria
#'
#' @param abundance A vector of relative abundances of bacterial taxons for a single sample
#' @param lineage Name of taxonomy lineage for each relative abundance in a sample (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia etc.)
#'
#' @return A dataframe of taxa abundances and whether they are susceptible or resistant to anti-gram negative antibiotics
#' @export
#'
#' @examples
#' do.call(rbind, apply(abx_test_df, 2, gram_neg_list, row.names(abx_test_df)))
#'
gram_neg_list <- function(abundance, lineage) {
  idx <- "gram_positive"
  suscept_vector <- !is_susceptible(lineage, idx)

  sorted_abundance_suscept <- sort(abundance[suscept_vector > 0.5], index.return=TRUE, decreasing = TRUE)
  susceptibles <- data.frame(lineage = character(), abundance = numeric(), phenotype = character())
  if(length(sorted_abundance_suscept$x) != 0) {
    susceptibles <- data.frame(lineage = lineage[suscept_vector > 0.5][sorted_abundance_suscept$ix],
                               abundance = sorted_abundance_suscept$x,
                               phenotype = "susceptible")
  }

  sorted_abundance_resist <- sort(abundance[suscept_vector < 0.5], index.return=TRUE, decreasing = TRUE)
  resistances <- data.frame(lineage = character(), abundance = numeric(), phenotype = character())
  if(length(sorted_abundance_resist$x) != 0) {
    resistances <- data.frame(lineage = lineage[suscept_vector < 0.5][sorted_abundance_resist$ix],
                              abundance = sorted_abundance_resist$x,
                              phenotype = "resistant")
  }
  rbind(susceptibles, resistances, make.row.names = FALSE)
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

#' Function to return susceptibility or resistance to antibiotics targeting anaerobes
#'
#' @param abundance A vector of relative abundances of bacterial taxons for a single sample
#' @param lineage Name of taxonomy lineage for each relative abundance in a sample (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia etc.)
#'
#' @return A dataframe of taxa abundances and whether they are susceptible or resistant to anti-anaerobic antibiotics
#' @export
#'
#' @examples
#' do.call(rbind, apply(abx_test_df, 2, anaerobe_list, row.names(abx_test_df)))
#'
anaerobe_list <- function(abundance, lineage) {
  idx <- "anaerobe"
  suscept_vector <- is_susceptible(lineage, idx)

  sorted_abundance_suscept <- sort(abundance[suscept_vector > 0.5], index.return=TRUE, decreasing = TRUE)
  susceptibles <- data.frame(lineage = character(), abundance = numeric(), phenotype = character())
  if(length(sorted_abundance_suscept$x) != 0) {
    susceptibles <- data.frame(lineage = lineage[suscept_vector > 0.5][sorted_abundance_suscept$ix],
                               abundance = sorted_abundance_suscept$x,
                               phenotype = "susceptible")
  }

  sorted_abundance_resist <- sort(abundance[suscept_vector < 0.5], index.return=TRUE, decreasing = TRUE)
  resistances <- data.frame(lineage = character(), abundance = numeric(), phenotype = character())
  if(length(sorted_abundance_resist$x) != 0) {
    resistances <- data.frame(lineage = lineage[suscept_vector < 0.5][sorted_abundance_resist$ix],
                              abundance = sorted_abundance_resist$x,
                              phenotype = "resistant")
  }
  rbind(susceptibles, resistances, make.row.names = FALSE)
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

#' Function to return susceptibility or resistance to antibiotics targeting aerobes
#'
#' @param abundance A vector of relative abundances of bacterial taxons for a single sample
#' @param lineage Name of taxonomy lineage for each relative abundance in a sample (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia etc.)
#'
#' @return A dataframe of taxa abundances and whether they are susceptible or resistant to anti-aerobic antibiotics
#' @export
#'
#' @examples
#' do.call(rbind, apply(abx_test_df, 2, aerobe_list, row.names(abx_test_df)))
#'
aerobe_list <- function(abundance, lineage) {
  idx <- "aerobe"
  suscept_vector <- is_susceptible(lineage, idx)

  sorted_abundance_suscept <- sort(abundance[suscept_vector > 0.5], index.return=TRUE, decreasing = TRUE)
  susceptibles <- data.frame(lineage = character(), abundance = numeric(), phenotype = character())
  if(length(sorted_abundance_suscept$x) != 0) {
    susceptibles <- data.frame(lineage = lineage[suscept_vector > 0.5][sorted_abundance_suscept$ix],
                               abundance = sorted_abundance_suscept$x,
                               phenotype = "susceptible")
  }

  sorted_abundance_resist <- sort(abundance[suscept_vector < 0.5], index.return=TRUE, decreasing = TRUE)
  resistances <- data.frame(lineage = character(), abundance = numeric(), phenotype = character())
  if(length(sorted_abundance_resist$x) != 0) {
    resistances <- data.frame(lineage = lineage[suscept_vector < 0.5][sorted_abundance_resist$ix],
                              abundance = sorted_abundance_resist$x,
                              phenotype = "resistant")
  }
  rbind(susceptibles, resistances, make.row.names = FALSE)
}

#this has not been tested yet
#' Function to calculate antibiotics index targeting gram-negative aerobes such as aminoglycoside
#'
#' @param abundance A vector of relative abundances of bacterial taxons for a single sample
#' @param lineage Name of taxonomy lineage for each relative abundance in a sample (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia etc.)
#'
#' @return The calculated antibiotics index for the sample
#' @export
#'
#' @examples
#' apply(abx_test_df, 2, aminoglycoside_index, row.names(abx_test_df))
#'
aminoglycoside_index <- function(abundance, lineage) {
  gram_neg_idx <- "gram_positive"
  aerobe_idx <- "aerobe"
  aminoglycoside_idx <- "aminoglycoside"
  gram_neg_vector <- !is_susceptible(lineage, gram_neg_idx)
  aerobe_vector <- is_susceptible(lineage, aerobe_idx)
  aminoglycoside_vector <- is_susceptible(lineage, aminoglycoside_idx)
  suscept_vector <- gram_neg_vector&aerobe_vector|aminoglycoside_vector
  calc_index(abundance, suscept_vector)
}

#' Function to return taxon susceptible or resistant to antibiotics targeting gram-negative aerobes
#'
#' @param abundance A vector of relative abundances of bacterial taxons for a single sample
#' @param lineage Name of taxonomy lineage for each relative abundance in a sample (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia etc.)
#'
#' @return A dataframe of taxa abundances and whether they are susceptible or resistant to antibiotics targeting gram-negative aerobes
#' @export
#'
#' @examples
#' do.call(rbind, apply(abx_test_df, 2, aminoglycoside_list, row.names(abx_test_df)))
#'
aminoglycoside_list <- function(abundance, lineage) {
  idx <- "penicillin"
  suscept_vector <- is_susceptible(lineage, idx)

  sorted_abundance_suscept <- sort(abundance[suscept_vector > 0.5], index.return=TRUE, decreasing = TRUE)
  susceptibles <- data.frame(lineage = character(), abundance = numeric(), phenotype = character())
  if(length(sorted_abundance_suscept$x) != 0) {
    susceptibles <- data.frame(lineage = lineage[suscept_vector > 0.5][sorted_abundance_suscept$ix],
                               abundance = sorted_abundance_suscept$x,
                               phenotype = "susceptible")
  }

  sorted_abundance_resist <- sort(abundance[suscept_vector < 0.5], index.return=TRUE, decreasing = TRUE)
  resistances <- data.frame(lineage = character(), abundance = numeric(), phenotype = character())
  if(length(sorted_abundance_resist$x) != 0) {
    resistances <- data.frame(lineage = lineage[suscept_vector < 0.5][sorted_abundance_resist$ix],
                              abundance = sorted_abundance_resist$x,
                              phenotype = "resistant")
  }
  rbind(susceptibles, resistances, make.row.names = FALSE)
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
