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

#' Antibiotic susceptibility of bacterial taxa
#' @format A data frame with the following columns:
#' \describe{
#'   \item{taxon}{The name of the taxon}
#'   \item{rank}{The rank of the taxon}
#'   \item{antibiotic}{The antibiotic or antibiotic class}
#'   \item{value}{
#'     The susceptibility of the taxon to the antibiotic, one of "susceptible"
#'     or "resistant".}
#'   \item{doi}{DOI of the publication from which the information was obtained.}
#' }
"taxon_susceptibility"

#' Example data from Weiss et al.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{sample_id}{The name of the sample}
#'   \item{lineage}{
#'      The taxonomic assignments from QIIME2. We removed any species-level
#'      taxa and blank taxa ("g__", etc.) from the assignments. Otherwise,
#'      the assignments are in the same format as they were when generated
#'      by the software.}
#'   \item{proportion}{The proportion of the lineage in the sample}
#'   \item{subject_id}{The subject identifier}
#'   \item{study_group}{The group, either "Healthy" or "Sepsis"}
#'   \item{study_window}{
#'     The time window in which the sample was collected. The study window is
#'     \code{NA} for healthy subjects.}
#' }
#' @references
#' Weiss SL, Bittinger K, Lee JJ, Friedman ES, Mattei LM, Graham K, Zhang D,
#' Bush J, Balamuth F, McGowan FX Jr, Bushman FD, Baldassano RN, Wu GD,
#' Wallace DC, Collman RG. Decreased Intestinal Microbiome Diversity in
#' Pediatric Sepsis: A Conceptual Framework for Intestinal Dysbiosis to
#' Influence Immunometabolic Function. Crit Care Explor. 2021 Mar
#' 17;3(3):e0360. doi: 10.1097/CCE.0000000000000360. PMID: 33786436; PMCID:
#' PMC7994045.
"weiss2021_data"
