#' Antibiotics index table from based on taxonomic levels from the Living Tree Project
#'
#' @description A dataframe of taxonomic levels, their phenotypes, and their antibiotics indices with references
#'
#' @format A \code{data.frame} with columns which are:
#' \describe{
#' \item{Kingdom}
#' \item{Phylum}
#' \item{Class}
#' \item{Order}
#' \item{Family}
#' \item{Genus}
#' \item{Species}
#' \item{gram_positive}{Index for gram positive bacteria}
#' \item{anaerobe}{Index for anaerobes}
#' \item{Phenotype_ref}{References for \code{gram_positive} and \code{anaerobe}}
#' \item{vancomycin}{Index for vancomycin susceptibility}
#' \item{vancomycin_ref}{References for \code{vancomycin}}
#' }
"LTP"

#' A test dataframe
#'
#' @description A dataframe with samples as columns and relative abundances of bacteria as rows
#'
"test_df"
