#' Function to calculate antibiotics index
#'
#' @param df A dataframe of samples as columns and bacterial taxons as rows. Labeling of taxons in rows should follow Greengenes annotations and be separated by double underscores (e.g. k__Bacteria;p__Bacteroidetes)
#' @param abx Antibiotics of interest to calculate index
#'
#' @return The calculated antibiotics index for each sample
#' @export
#'
#' @examples
#' \dontrun{
#' abx_index(test_df)
#' }
#'
abx_index <- function(df, abx = NULL) {
  taxa <- row.names(df)
  suscept_vector <- get_susceptibles(taxa, abx)

  return(lapply(df, calc_index, suscept_vector))
}

calc_index <- function(sample_vector, suscept_vector) {
  sum_suscept_taxa <- sum(sample_vector[suscept_vector])
  return(log10(sum_suscept_taxa/(1-sum_suscept_taxa)))
}

#write the antibiotis part as an option for user, not allowing user to pick an abx

get_susceptibles <- function(taxa, abx) {
  taxa <- gsub(";", " ", taxa)
  if(is.null(abx)) {
    stop("No antibiotics specified to test")
  }
  if(abx %in% c("Glycopeptides", "Macrolides", "Oxazolidinones", "Lincosamides", "Lipopeptides", "Amoxicillin")) {
    susceptible_vector <- gram_positives(taxa)
  }
  else if(abx == "Vancomycin") {
    susceptible_vector <- vanco(taxa)
  }
  else if(abx %in% c("Polymyxins", "Aztreonam")) {
    susceptible_vector <- !gram_positives(taxa)
  }
  else if(abx =="Nitroimidazoles") {
    susceptible_vector <- anaerobes(taxa)
  }
  else if(abx =="Fluoroquinolones") {
    susceptible_vector <- aerobes(taxa)
  }
  else {
    stop("Not a listed abx")
  }
  return(susceptible_vector)
}

#implement this as phylogenetic binary tree using greengenes nomenclature

gram_pos <- c("p__Actinobacteria", "p__Firmicutes")
firm_exceptions <- c("c__Mollicutes", "c__Negativicutes")
##in LTP nomenclature, need to define the Negativicutes (Veillonella and Lactobacillus)

oblig_anaero_gn <- c("o__Bacteroidales", "g__Bacteroides", "g__Fusobacterium", "g__Porphyromonas",
                     "g__Prevotella", "g__Veillonella", "g__Parabacteroides")

oblig_anaero_gp <- c("g__Actinomyces", "g__Clostridium", "g__Peptostreptococcus", "g__Propionibacterium",
                     "Peptoniphilus", "g__Finegoldia", "g__Anaerococcus", "g__Parvimonas", "g__Lactobacillus",
                     "g__Bifidobacterium", "g__Eubacterium", "g__Gemella")

aerotol_anaero <- c("g__Staphylococcus", "g__Streptococcus", "g__Salmonella", "g__Listeria",
                    "g__Escherichia s__coli", "g__Shewanella s__oneidensis", "g__Yersinia s__pestis")

aero <- c("g__Mycobacterium s__tuberculosis", "g__Nocardia s__asteroides")

vanco <- function(taxa) {
  ##Bacteria that are known to be susceptible to vancomycin are gram positive Firmicutes and Actinobacteria, with some exceptions
  ##This resistant list is ignoring intrinsic efflux pumps
  susceptibles <- c(gram_pos, "c__Bacteroidia") ##https://doi.org/10.1126/sciadv.aax2358; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5161046/
  exceptions <- c(firm_exceptions,
                  "g__Lactobacillus s__rhamnosus", "g__Lactobacillus s__paracasei", "g__Lactobacillus s__plantarum",
                  "g__Lactobacillus s__reuteri", "g__Lactobacillus s__fermentum", ##https://doi.org/10.1093/jac/dkm035
                  "g__Enterococcus s__gallinarum", "g__Enterococcus s__casseliflavus", "g__Enterococcus s__flavescens") ##https://doi.org/10.1016/j.jiac.2018.01.001

  vanco_suscept <- grepl(paste0(susceptibles, collapse = "|"), taxa)
  vanco_resist <- !grepl(paste0(exceptions, collapse = "|"), taxa)

  return(vanco_suscept & vanco_resist)
}

##tetracycline isn't implemented yet
tetracycline <- function(taxa) { ##excluding resistance through acquired mobile elements, only trying to get intrinsic resistance
  ##All bacteria are theoretically susceptible to tetracycline except these bacteria for which some resistant clinical strains were isolated
  resistance <- c("g__Acinetobacter s__baumannii",
                      "g__Bacteroides s__fragilis",
                      "g__Escherichia s__coli",
                      "g__Enterobacter",
                      "g__Enterococcus s__faecalis",
                      "g__Klebsiella pneumoniae", ##Klebsiella generally susceptible
                      "g__Pseudomonas s__aeruginosa",
                      "g__Proteus s__mirabilis",
                      "g__Staphylococcus s__aureus",
                      "g__Stenotrophomonas s__maltophilia",
                      "g__Serratia s__marcescens", ##Intrinsic bacterial multidrug efflux pumps: doi:10.1101/cshperspect.a025387, Table 2
                      "g__Salmonella s__typhimurium", "g__Campylobacter s__jejuni", ##DOI: 10.1038/nrmicro1464
                      "g__Bacteroides") #80% of Bacteroides are resistant to tetracyclines due to mobile elements: DOI: 10.1128/mBio.00569-13

  tetra_resist <- !grepl(paste0(resistance, collapse = "|"), taxa)

  return(tetra_resist)

  ##Include distribution of tet protection genes on mobile elements? http://faculty.washington.edu/marilynr/, https://doi.org/10.1016/j.femsle.2005.02.034
  ##Can use the Clinical and Laboratory Standards Institute (CLSI) guideline for assessing resistance?
}

gram_positives <- function(taxa) {
  susceptibles <- grepl(paste0(gram_pos, collapse = "|"), taxa)
  exceptions <- !grepl(paste0(firm_exceptions, collapse = "|"), taxa)
  return(susceptibles & exceptions)
}

anaerobes <- function(taxa) {
  susceptibles <- grepl(paste0(c(oblig_anaero_gn, oblig_anaero_gp), collapse = "|"), taxa)
  return(susceptibles)
}

aerobes <- function(taxa) {
  susceptibles <- grepl(paste0(c(aerotol_anaero, aero), collapse = "|"), taxa)
  return(susceptibles)
}

##review paper basing drug class: https://doi.org/10.1016/j.bcp.2017.01.003

#TODO: Import library that can extend the full taxonomic levels in row names
