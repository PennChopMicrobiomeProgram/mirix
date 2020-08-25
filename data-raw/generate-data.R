library(dplyr)
library(tidyr)
library(taxonomizr)

selecting_LTP_cols <- function(file) {
  file %>% select(V5, V11, V1) %>%
    mutate(taxon = paste0(V11, ";", V5)) %>%
    mutate(semi = nchar(gsub("[^\\;]", "", taxon)))
}

##can rewrite this to be just given one column, do later
mutate_n_select <- function(df) {
  df %>%
    mutate(anaerobe = ifelse((is.na(anaerobe.y) | aerobic_status == "not indicated"), anaerobe.x, anaerobe.y)) %>%
    mutate(gram_positive = ifelse((is.na(gram_positive.y) | gram_stain == "not indicated"), gram_positive.x, gram_positive.y)) %>%
    select(-c(anaerobe.x, anaerobe.y, gram_positive.x, gram_positive.y, aerobic_status, gram_stain))
}

LTP_file <- utils::read.delim(here::here("data-raw", "LTPs128_SSU.tsv"), sep = "\t", header = FALSE)

tax_fp <- "/Users/tuv/Google Drive/db/accessionTaxa.sql"
blastAccessions <- LTP_file %>%
  mutate(V1 = as.character(paste0(V1, ".1"))) %>%
  pull(V1)
ids <- accessionToTaxa(blastAccessions, tax_fp)

LTP_missing <- LTP_file %>%
  selecting_LTP_cols() %>%
  filter((semi < 6)|(grepl("Unclassified|Unnamed", taxon))) %>%
  distinct() %>%
  select(taxon, V1)

LTP_5 <- as.data.frame(getTaxonomy(ids, tax_fp)) %>%
  cbind(blastAccessions) %>%
  filter(!is.na(superkingdom)) %>%
  mutate(blastAccessions = gsub("\\.1", "", blastAccessions)) %>%
  merge(LTP_missing, by.x = "blastAccessions", by.y = "V1", all.x = TRUE) %>%
  filter(!is.na(taxon)) %>%
  select(superkingdom, phylum, class, order, family, genus, species) %>%
  rename(kingdom = superkingdom) %>%
  distinct()

colnames(LTP_5) <- paste(toupper(substring(colnames(LTP_5), 1,1)), substring(colnames(LTP_5), 2), sep="")

##Collapsing LTP subclass and suborder taxonomic levels
LTP_8 <- LTP_file %>%
  selecting_LTP_cols() %>%
  filter(!grepl("Unclassified|Unnamed", taxon)) %>%
  filter(semi == 8) %>%
  distinct() %>%
  select(taxon) %>%
  separate(., col = taxon, sep = ";", into = c("Kingdom", "Phylum", "Class", "Subclass", "Order", "Suborder", "Family", "Genus", "Species")) %>%
  select(-c(Subclass, Suborder))

LTP_7_subclass <- LTP_file %>%
  selecting_LTP_cols() %>%
  filter(!grepl("Unclassified|Unnamed", taxon)) %>%
  filter(semi == 7) %>%
  distinct() %>%
  select(taxon) %>%
  separate(., col = taxon, sep = ";", into = c("Kingdom", "Phylum", "Class", "Subclass", "Order", "Family", "Genus", "Species")) %>%
  mutate(subclass = grepl("dae$", Subclass)) %>%
  filter(subclass) %>%
  select(-c(Subclass, subclass))

LTP_7_suborder <- LTP_file %>%
  selecting_LTP_cols() %>%
  filter(!grepl("Unclassified|Unnamed", taxon)) %>%
  filter(semi == 7) %>%
  distinct() %>%
  select(taxon) %>%
  separate(., col = taxon, sep = ";", into = c("Kingdom", "Phylum", "Class", "Order", "Suborder", "Family", "Genus", "Species")) %>%
  mutate(subord = grepl("neae$", Suborder)) %>%
  filter(subord) %>%
  select(-c(Suborder, subord))

##manually curated genus phenotypes from Ceylan
phenotypes_genus <- utils::read.delim(here::here("data-raw", "genera_0831.txt"), sep = "\t", header = TRUE) %>%
  select(name, aerobic_status, gram_stain, doi) %>%
  #mutate(level = "Genus") %>%
  mutate(anaerobe = ifelse(grepl("anaerobe", aerobic_status), ifelse(grepl("aerobe, microaerobe, or anaerobe|facultative anaerobe", aerobic_status), FALSE, TRUE), FALSE)) %>%
  mutate(gram_positive = ifelse(gram_stain == "Gram-positive", TRUE, FALSE)) %>%
  mutate(doi = as.character(doi)) %>%
  rename(Genus = name) %>%
  rename(Phenotype_ref = doi)

##manually curated species phenotypes from Ceylan
phenotypes_species <- utils::read.delim(here::here("data-raw", "species_0831.txt"), sep = "\t", header = TRUE) %>%
  select(name, aerobic_status, gram_stain, doi) %>%
  #mutate(level = "Species") %>%
  mutate(anaerobe = ifelse(grepl("anaerobe", aerobic_status), ifelse(grepl("aerobe, microaerobe, or anaerobe|facultative anaerobe", aerobic_status), FALSE, TRUE), FALSE)) %>%
  mutate(gram_positive = ifelse(gram_stain == "Gram-positive", TRUE, FALSE)) %>%
  mutate(doi = as.character(doi)) %>%
  rename(Species = name) %>%
  rename(Phenotype_ref = doi)

##Generate a df with taxonomic levels from the LTP
LTP <- LTP_file %>%
  selecting_LTP_cols() %>%
  filter(!grepl("Unclassified|Unnamed", taxon)) %>%
  filter(semi == 6) %>%
  distinct() %>%
  select(taxon) %>%
  separate(., col = taxon, sep = ";", into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  rbind(LTP_5, LTP_8, LTP_7_subclass, LTP_7_suborder) %>%
  distinct() %>%
  filter(!grepl("^Archaea$", Kingdom)) %>%
  arrange_all() %>%
  #mutate(uniq_id = row_number()) %>%
  ##Populate df with bacteria characteristics
  mutate(gram_positive = grepl("^Actinobacteria$|^Firmicutes$", Phylum)) %>%
  mutate(gram_positive = ifelse(grepl("^Negativicutes$", Class), FALSE, gram_positive)) %>%
  #anareobe classification reference: https://doi.org/10.1128/microbiolspec.dmih2-0015-2015
  mutate(anaerobe = (grepl("^Bacteroidales$", Order)|
                    (grepl(paste("^Fusobacterium",
                                 "Veillonella",
                                 "Actinomyces",
                                 "Clostridium",
                                 "Peptostreptococcus",
                                 "Propionibacterium",
                                 "Peptoniphilus",
                                 "Finegoldia",
                                 "Anaerococcus",
                                 "Parvimonas",
                                 "Lactobacillus",
                                 "Bifidobacterium",
                                 "Eubacterium",
                                 "Gemella$", sep = "$|^"), Genus)))) %>%
  merge(select(phenotypes_genus, Genus, anaerobe, aerobic_status, gram_positive, gram_stain, Phenotype_ref), by = "Genus", all.x = TRUE) %>%
  mutate_n_select() %>%
  merge(select(phenotypes_species, Species, anaerobe, aerobic_status, gram_positive, gram_stain, Phenotype_ref), by = "Species", all.x = TRUE) %>%
  mutate_n_select() %>%
  mutate(Phenotype_ref = ifelse(is.na(Phenotype_ref.y), Phenotype_ref.x, Phenotype_ref.y)) %>%
  select(-c(Phenotype_ref.x, Phenotype_ref.y)) %>%
  mutate(vancomycin = ((gram_positive) |
                      (grepl("^Bacteroidia$", Class)))& ##https://doi.org/10.1126/sciadv.aax2358; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5161046/
                      (!grepl(paste("Lactobacillus rhamnosus",
                                    "Lactobacillus paracasei",
                                    "Lactobacillus plantarum",
                                    "Lactobacillus reuteri",
                                    "Lactobacillus fermentum", ##https://doi.org/10.1093/jac/dkm035
                                    "Enterococcus gallinarum",
                                    "Enterococcus casseliflavus",
                                    "Enterococcus flavescens", sep = "|"), Species))) %>% ##https://doi.org/10.1016/j.jiac.2018.01.001
  mutate(vancomycin_ref = ifelse(grepl("^Bacteroidia$", Class), "https://doi.org/10.1126/sciadv.aax2358; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5161046/", NA)) %>%
  mutate(vancomycin_ref = ifelse(grepl(paste("Lactobacillus rhamnosus",
                                             "Lactobacillus paracasei",
                                             "Lactobacillus plantarum",
                                             "Lactobacillus reuteri",
                                             "Lactobacillus fermentum", sep = "|"), Species), "https://doi.org/10.1093/jac/dkm035", vancomycin_ref)) %>%
  mutate(vancomycin_ref = ifelse(grepl(paste("Enterococcus gallinarum",
                                             "Enterococcus casseliflavus",
                                             "Enterococcus flavescens", sep = "|"), Species), "https://doi.org/10.1016/j.jiac.2018.01.001", vancomycin_ref)) %>%
  mutate(nitroimidazole = anaerobe) %>%
  mutate(nitroimidazole_ref = NA) %>%
  mutate(fluoroquinolone = !anaerobe) %>%
  mutate(fluoroquinolone_ref = NA) %>%
  mutate(polymyxin_n_aztreonam = !gram_positive) %>%
  mutate(polymyxin_n_aztreonam_ref = NA) %>%
  mutate(glycopeptides_macrolides_oxazolidinones_lincosamides_lipopeptides_amoxicillin = gram_positive) %>%
  mutate(glycopeptides_macrolides_oxazolidinones_lincosamides_lipopeptides_amoxicillin_ref = NA)

#save this df as data to use in package
usethis::use_data(LTP, overwrite = TRUE)

test_df <- data.frame(a = c(0.1, 0.4, 0.1, 0.3, 0.1), b = c(0, 0.25, 0.25, 0.25, 0.25), c = c(1, 0, 0, 0, 0), d = c(0, 0, 1, 0, 0))
row.names(test_df) <- c("k__Bacteria; p__Actinobacteria; c__; o__; f__; g__; s__",
                        "k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Clostridiaceae; g__Clostridium; s__",
                        "k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Enterococcaceae; g__Enterococcus; s__gallinarum",
                        "k__Bacteria; p__Firmicutes; c__Negativicutes; o__Veillonellales; f__Veillonellaceae; g__Veillonella; s__parvula",
                        "k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales; f__Pseudomonadaceae; g__; s__"
)

usethis::use_data(test_df, overwrite = TRUE)

aerotol_anaero <- c("g__Staphylococcus", "g__Streptococcus", "g__Salmonella", "g__Listeria",
                    "g__Escherichia s__coli", "g__Shewanella s__oneidensis", "g__Yersinia s__pestis")

aero <- c("g__Mycobacterium s__tuberculosis", "g__Nocardia s__asteroides")
