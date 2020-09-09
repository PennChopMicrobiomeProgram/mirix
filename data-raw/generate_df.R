library(dplyr)
library(tidyr)

##manually curated genus phenotypes from Ceylan
phenotypes_genus <- utils::read.delim(here::here("data-raw", "genera_0831.txt"), sep = "\t", header = TRUE) %>%
  select(name, aerobic_status, gram_stain, doi) %>%
  unique() %>%
  mutate(level = "Genus") %>%
  mutate(anaerobe = ifelse(grepl("anaerobe", aerobic_status), ifelse(grepl("microaerobe or anaerobe|facultative anaerobe", aerobic_status), FALSE, TRUE), FALSE)) %>%
  mutate(aerobe = ifelse(grepl("anaerobe|not indicated|variable", aerobic_status), FALSE, TRUE)) %>%
  mutate(gram_positive = ifelse(gram_stain == "Gram-positive", TRUE, FALSE)) %>%
  mutate(gram_negative = ifelse(gram_stain == "Gram-negative", TRUE, FALSE)) %>%
  mutate(tetracycline = TRUE) %>%
  mutate(doi = as.character(doi), name = as.character(name)) %>%
  select(-c(aerobic_status, gram_stain)) %>%
  pivot_longer(cols = c("anaerobe", "aerobe", "gram_positive", "gram_negative", "tetracycline"), values_to = "boo", names_to = "attribute")
phenotypes_genus[phenotypes_genus$name=="Gemella"&phenotypes_genus$attribute=="gram_positive", "boo"] <- TRUE
phenotypes_genus[phenotypes_genus$name=="Gemella"&phenotypes_genus$attribute=="gram_negative", "boo"] <- FALSE
phenotypes_genus <- phenotypes_genus[-which(phenotypes_genus$name=="Bifidobacterium"&phenotypes_genus$attribute=="anaerobe"), ]

##manually curated species phenotypes from Ceylan
phenotypes_species <- utils::read.delim(here::here("data-raw", "species_0831.txt"), sep = "\t", header = TRUE) %>%
  select(name, aerobic_status, gram_stain, doi) %>%
  unique() %>%
  mutate(level = "Species") %>%
  mutate(anaerobe = ifelse(grepl("anaerobe", aerobic_status), ifelse(grepl("aerobe, microaerobe, or anaerobe|microaerobe or anaerobe|facultative anaerobe", aerobic_status), FALSE, TRUE), FALSE)) %>%
  mutate(aerobe = ifelse(grepl("anaerobe|not indicated|aerotolerant|variable", aerobic_status), FALSE, TRUE)) %>%
  mutate(gram_positive = ifelse(gram_stain == "Gram-positive", TRUE, FALSE)) %>%
  mutate(gram_negative = ifelse(gram_stain == "Gram-negative", TRUE, FALSE)) %>%
  mutate(tetracycline = TRUE) %>%
  mutate(doi = as.character(doi), name = as.character(name)) %>%
  select(-c(aerobic_status, gram_stain)) %>%
  pivot_longer(cols = c("anaerobe", "aerobe", "gram_positive", "gram_negative", "tetracycline"), values_to = "boo", names_to = "attribute")

grep_df <- rbind(phenotypes_genus, phenotypes_species)
grep_df <- grep_df[c("attribute", "boo", "name", "level", "doi")]

##gram_positive
grep_df[nrow(grep_df)+1,] <- list("gram_positive", TRUE, "Actinobacteria", "Phylum", "10.1128/MMBR.00019-15")
grep_df[nrow(grep_df)+1,] <- list("gram_positive", TRUE, "Firmicutes", "Phylum", "10.1099/00207713-28-1-1")
grep_df[nrow(grep_df)+1,] <- list("gram_positive", FALSE, "Negativicutes", "Class", "10.4056/sigs.2981345")
grep_df[nrow(grep_df)+1,] <- list("gram_positive", FALSE, "Salmonella", "Genus", "10.1002/9781118960608.gbm01166")

##gram_negative
grep_df[nrow(grep_df)+1,] <- list("gram_negative", FALSE, "Actinobacteria", "Phylum", "10.1128/MMBR.00019-15")
grep_df[nrow(grep_df)+1,] <- list("gram_negative", FALSE, "Firmicutes", "Phylum", "10.1099/00207713-28-1-1")
grep_df[nrow(grep_df)+1,] <- list("gram_negative", TRUE, "Negativicutes", "Class", "10.4056/sigs.2981345")
grep_df[nrow(grep_df)+1,] <- list("gram_negative", TRUE, "Salmonella", "Genus", "10.1002/9781118960608.gbm01166")

##vancomycin
grep_df[nrow(grep_df)+1,] <- list("vancomycin", TRUE, "Bacteroidia", "Class", "10.1126/sciadv.aax2358; 10.1093/jac/dkw383") ##susceptible

vanco_lacto_except <- c("Lactobacillus rhamnosus",
                        "Lactobacillus paracasei",
                        "Lactobacillus plantarum",
                        "Lactobacillus reuteri",
                        "Lactobacillus fermentum")
for(lacto in vanco_lacto_except){
  grep_df[nrow(grep_df)+1,] <- list("vancomycin", FALSE, lacto, "Species", "10.1093/jac/dkm035") ##resistant
}

vanco_entero_except <- c("Enterococcus gallinarum",
                         "Enterococcus casseliflavus",
                         "Enterococcus flavescens")
for(entero in vanco_entero_except){
  grep_df[nrow(grep_df)+1,] <- list("vancomycin", FALSE, entero, "Species", "10.1016/j.jiac.2018.01.001") ##resistant
}


##anaerobes
grep_df[nrow(grep_df)+1,] <- list("anaerobe", TRUE, "Bacteroidales", "Order", "10.1128/microbiolspec.dmih2-0015-2015")
grep_df[nrow(grep_df)+1,] <- list("anaerobe", FALSE, "Salmonella", "Genus", "10.1002/9781118960608.gbm01166")
grep_df[nrow(grep_df)+1,] <- list("anaerobe", FALSE, "Listeria", "Genus", "10.1002/9781118960608.gbm00547")
anaerobe <- c("Actinomyces",
              "Peptostreptococcus",
              "Finegoldia",
              "Parvimonas",
              "Bifidobacterium")
for(ana in anaerobe){
  grep_df[nrow(grep_df)+1,] <- list("anaerobe", TRUE, ana, "Genus", "10.1128/microbiolspec.dmih2-0015-2015")
}

##aerobes
grep_df[nrow(grep_df)+1,] <- list("aerobe", FALSE, "Bacteroidales", "Order", "10.1128/microbiolspec.dmih2-0015-2015")
grep_df[nrow(grep_df)+1,] <- list("aerobe", FALSE, "Salmonella", "Genus", "10.1002/9781118960608.gbm01166")
grep_df[nrow(grep_df)+1,] <- list("aerobe", FALSE, "Listeria", "Genus", "10.1002/9781118960608.gbm00547")
for(ana in anaerobe){
  grep_df[nrow(grep_df)+1,] <- list("aerobe", FALSE, ana, "Genus", "10.1128/microbiolspec.dmih2-0015-2015")
}

##tetracycline
tet_intrins_resist_replace <- c("Bacteroides fragilis",
                                "Escherichia coli",
                                "Enterobacter",
                                "Enterococcus faecalis",
                                "Pseudomonas aeruginosa",
                                "Staphylococcus aureus",
                                "Stenotrophomonas maltophilia")
for(tet_replace in tet_intrins_resist_replace){
  grep_df[grep_df$name==tet_replace&grep_df$attribute=="tetracycline", "boo"] <- FALSE
  grep_df[grep_df$name==tet_replace&grep_df$attribute=="tetracycline", "doi"] <- "10.1101/cshperspect.a025387"
}

tet_intrins_resist_species <- c("Acinetobacter baumannii",
                                "Klebsiella pneumoniae", ##Klebsiella generally susceptible
                                "Proteus mirabilis",
                                "Serratia marcescens")
for(tet_species in tet_intrins_resist_species){
  grep_df[nrow(grep_df)+1,] <- list("tetracycline", FALSE, tet_species, "Species", "10.1101/cshperspect.a025387")
}

grep_df[nrow(grep_df)+1,] <- list("tetracycline", FALSE, "Salmonella typhimurium", "Species", "10.1038/nrmicro1464")
grep_df[nrow(grep_df)+1,] <- list("tetracycline", FALSE, "Campylobacter jejuni", "Species", "10.1038/nrmicro1464")
grep_df[grep_df$name=="Bacteroides"&grep_df$attribute=="tetracycline", "boo"] <- FALSE
grep_df[grep_df$name=="Bacteroides"&grep_df$attribute=="tetracycline", "doi"] <- "10.1128/mBio.00569-13"
grep_df[nrow(grep_df)+1,] <- list("tetracycline", TRUE, "Klebsiella", "Genus", "10.1101/cshperspect.a025387")

#save this df as data to use in package
usethis::use_data(grep_df, overwrite = TRUE)


##review paper on drug class: https://doi.org/10.1016/j.bcp.2017.01.003
##Tet protection genes on mobile elements: http://faculty.washington.edu/marilynr/, https://doi.org/10.1016/j.femsle.2005.02.034
##the Clinical and Laboratory Standards Institute (CLSI) guidelines assesses abx resistance


abx_test_df <- data.frame(a = c(0.1, 0.4, 0.1, 0.3, 0.1),
                          b = c(0, 0.25, 0.25, 0.25, 0.25),
                          c = c(1, 0, 0, 0, 0),
                          d = c(0, 0, 1, 0, 0))
row.names(abx_test_df) <- c("k__Bacteria; p__Actinobacteria; c__; o__; f__; g__; s__",
                        "k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Clostridiaceae; g__Clostridium; s__Clostridium perfringens",
                        "k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Enterococcaceae; g__Enterococcus; s__Enterococcus gallinarum",
                        "k__Bacteria; p__Firmicutes; c__Negativicutes; o__Veillonellales; f__Veillonellaceae; g__Veillonella",
                        "k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales; f__Pseudomonadaceae; g__; s__"
)

usethis::use_data(abx_test_df, overwrite = TRUE)
