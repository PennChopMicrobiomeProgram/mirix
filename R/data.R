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

#' Gram stain and aerobic status of bacterial taxa
#' @format A data frame with the following columns:
#' \describe{
#'   \item{taxon}{The name of the taxon}
#'   \item{rank}{The rank of the taxon}
#'   \item{aerobic_status}{
#'     The aerobic status. One of "aerobe", "facultative anaerobe", or
#'     "obligate anaerobe".}
#'   \item{gram_stain}{
#'     How the taxon appears when Gram-stained. One of "Gram-positive" or
#'     "Gram-negative".}
#'   \item{doi}{DOI of the publication from which the information was obtained.}
#' }
"taxon_phenotypes"
