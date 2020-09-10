#' A dataframe of bacterial taxonomy ranks and their cell wall, oxygen tolerance, and antibiotics characteristics
#'
#' @description A dataframe of taxonomic ranks, their phenotypes, and their antibiotics indices with references
#'
#' @format A \code{data.frame} with columns which are:
#'
#' \describe{
#'   \item{attribute}{A characteristic of the bacteria taxonomic rank}
#'   \item{boo}{A boolean determinant about the attribute}
#'   \item{name}{Name of the taxon}
#'   \item{rank}{Rank of the taxon}
#'   \item{doi}{Reference of the attribute}
#' }
"abx_idx_df"

#' A test dataframe
#'
#' @description A dataframe with samples as columns and relative abundances of bacteria as rows
#'
"abx_test_df"
