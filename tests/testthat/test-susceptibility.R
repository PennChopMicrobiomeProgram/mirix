gram_stain_db <- data.frame(
  taxon = c(
    "Actinobacteria", "Bacteroidetes", "Firmicutes",
    "Negativicutes", "Proteobacteria"),
  rank = c("Phylum", "Phylum", "Phylum", "Class", "Phylum"),
  gram_stain = c(
    "Gram-positive", "Gram-negative", "Gram-positive",
    "Gram-negative", "Gram-negative"),
  stringsAsFactors = FALSE)

gram_stain_aerobic_status_db <- data.frame(
  taxon = c(
    "Rothia aerobat", "Enterococcus", "Blautia", "Lactobacillus acidophilus",
    "Bradyrhizobium", "Rothia", "Actinobacteria"),
  rank = c(
    "Species", "Genus", "Genus", "Species", "Genus", "Genus", "Phylum"),
  aerobic_status = c(
    "aerobe", "facultative anaerobe", "obligate anaerobe", NA, "aerobe",
    "aerobe", NA),
  gram_stain = c(
    "Gram-positive", "Gram-positive", "Gram-positive", "Gram-positive",
    "Gram-negative", NA, "Gram-positive"),
  stringsAsFactors = FALSE)

vancomycin_db <- data.frame(
  taxon = c(
    "Lactobacillus", "Enterococcus flavescens", "Lactobacillus delbrueckii",
    "Lactobacillus lactis", "Lactobacillus casei"),
  rank = c("Genus", "Species", "Species", "Species", "Species"),
  value = c(
    "resistant", "resistant", "susceptible", "susceptible", "resistant"),
  antibiotic = "vancomycin",
  stringsAsFactors = FALSE)

tetracycline_db <- data.frame(
  taxon = c("Escherichia coli", "Enterococcus faecalis", "Klebsiella"),
  rank = c("Species", "Species", "Genus"),
  value = c("resistant", "resistant", "susceptible"),
  antibiotic = "tetracycline",
  stringsAsFactors = FALSE)

penicillin_db <-data.frame(
  taxon = c(
    "Bacteroides", "Stenotrophomonas maltophilia",
    "Enterococcus faecalis", "Enterococcus faecium",
    "Escherichia coli", "Mycoplasma"),
  rank = c("Genus", "Species", "Species", "Species", "Species", "Genus"),
  value = c(
    "resistant", "resistant", "resistant", "resistant", "resistant",
    "resistant"),
  antibiotic = "penicillin",
  stringsAsFactors = FALSE)

aminoglycoside_db <- data.frame(
  taxon = "Staphylococcus",
  rank = "Genus",
  value = "susceptible",
  antibiotic = "aminoglycoside",
  stringsAsFactors = FALSE)

test_antibiotic_db <- rbind(
  vancomycin_db, tetracycline_db, penicillin_db, aminoglycoside_db)

test_that("vancomycin_susceptibility works for selected taxa", {
  lineage <- c(
    "Firmicutes; Enterococcus faecalis",
    "Firmicutes; Lactobacillus",
    "Firmicutes; Lactobacillus; Lactobacillus delbrueckii")
  expect_equal(
    vancomycin_susceptibility(
      lineage = lineage,
      antibiotic_db = test_antibiotic_db,
      phenotype_db = gram_stain_db),
    c("susceptible", "resistant", "susceptible"))
})

test_that("tetracycline_susceptibility works for selected taxa", {
  lineage <- c(
    "Firmicutes; Enterococcus faecalis",
    "Firmicutes; Lactobacillus",
    "Proteobacteria; Klebsiella; Klebsiella pneumoniae")
  expect_equal(
    tetracycline_susceptibility(
      lineage = lineage,
      antibiotic_db = test_antibiotic_db),
    c("resistant", "susceptible", "susceptible"))
})

test_that("penicillin_susceptibility works for selected taxa", {
  lineage <- c(
    "Firmicutes; Enterococcus faecalis",
    "Firmicutes; Lactobacillus",
    "Proteobacteria; Klebsiella; Klebsiella pneumoniae")
  expect_equal(
    penicillin_susceptibility(
      lineage = lineage,
      antibiotic_db = test_antibiotic_db),
    c("resistant", "susceptible", "susceptible"))
})

test_that("aminoglycoside_susceptibility works for selected taxa", {
  lineage <- c(
    "Staphylococcus epidermidis",
    "Rothia", # Not enough info to determine both phenotypes, should be NA
    "Actinobacteria; Rothia",
    "Bacteroides vulgatus")
  expect_equal(
    aminoglycoside_susceptibility(
      lineage = lineage,
      antibiotic_db = test_antibiotic_db,
      phenotype_db = gram_stain_aerobic_status_db),
    c("susceptible", NA, "resistant", "resistant"))
})
