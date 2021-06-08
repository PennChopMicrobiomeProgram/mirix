gram_stain_db <- data.frame(
  taxon = c(
    "Actinobacteria", "Bacteroidetes", "Firmicutes",
    "Negativicutes", "Proteobacteria"),
  rank = c("Phylum", "Phylum", "Phylum", "Class", "Phylum"),
  value = c(
    "Gram-positive", "Gram-negative", "Gram-positive",
    "Gram-negative", "Gram-negative"),
  stringsAsFactors = FALSE)

vanco_db <- data.frame(
  taxon = c(
    "Lactobacillus", "Enterococcus flavescens", "Lactobacillus delbrueckii",
    "Lactobacillus lactis", "Lactobacillus casei"),
  rank = c("Genus", "Species", "Species", "Species", "Species"),
  value = c(
    "resistant", "resistant", "susceptible", "susceptible", "resistant"),
  stringsAsFactors = FALSE)

test_that("match_annotation works for antibiotic susceptibility", {
  expect_equal(
    match_annotation(
      c("Enterococcus", "Lactobacillus",
        "Lactobacillus delbrueckii subsp. lactis"),
      vanco_db),
    c(NA, "resistant", "susceptible"))
})

test_that("match_annotation works for normal input", {
  # For the lineage "Firmicutes Negativicutes", the phenotype value for
  # Negativicutes (Gram-negative) should override the value for Firmicutes
  # (Gram-positive)
  expect_equal(
    match_annotation(
      c("Bacteroidetes", "Firmicutes", "Firmicutes Negativicutes"),
      gram_stain_db),
    c("Gram-negative", "Gram-positive", "Gram-negative"))
})

test_that("match_annotation works for lineage of length 1", {
  expect_equal(
    match_annotation("Bacteroidetes", gram_stain_db),
    "Gram-negative")
})

test_that("match_annotation works for empty lineage vector", {
  expect_equal(match_annotation(character(), gram_stain_db), character())
})

test_that("match_annotation works for database of length 1", {
  expect_equal(
    match_annotation(
      c("Bacteroidetes", "Firmicutes", "Firmicutes Negativicutes"),
      gram_stain_db[2,]),
    c("Gram-negative", NA, NA))
})

test_that("match_annotation works for empty database", {
  expect_equal(
    match_annotation(
      c("Bacteroidetes", "Firmicutes", "Firmicutes Negativicutes"),
      gram_stain_db[integer(),]),
    as.character(c(NA, NA, NA)))
})

test_that("first_non_na_value works for multiple non-NA values", {
  expect_equal(first_non_na_value(c(NA, NA, "a", NA, "b", NA)), "a")
})

test_that("first_non_na_value returns NA if all values are NA", {
  na_chars <- as.character(c(NA, NA, NA, NA))
  expect_equal(first_non_na_value(na_chars), NA_character_)
})

test_that("match_taxa works for normal input", {
  expect_equal(
    match_taxa(
      c("Bacteroidetes; hgj", "Bacteria - Firmicutes", "kgj", "jasfdljh"),
      c("Firmicutes", "Bacteroidetes", "Proteobacteria")),
    c(2, 1, NA, NA))
})

test_that("match_taxa works for taxa with prefixes", {
  prefix_lineage <- c(
    "p__Bacteroidetes; hgj",
    "k__Bacteria - p__Firmicutes",
    "g__kgj",
    "s__jasfdljh")
  taxa <- c("Firmicutes", "Bacteroidetes", "Proteobacteria")
  expect_equal(
    match_taxa(prefix_lineage, taxa),
    c(2, 1, NA, NA))
})

test_that("match_taxa works for species names", {
  species_lineage <- c(
    "Bacteroidetes - Bacteroides vulgatus",
    "Enterobacteriaceae; Escherichia; Escherichia coli",
    "kgj vulgatus",
    "Bacteroidetes - Bacteroides")
  species_taxa <- c(
    "Escherichia coli",
    "Bacteroides vulgatus")
  expect_equal(
    match_taxa(species_lineage, species_taxa),
    c(2, 1, NA, NA))
})

test_that("match_taxa works for lineage vectors of length 1", {
  expect_equal(
    match_taxa(
      "Bacteroidetes; hgj",
      c("Firmicutes", "Bacteroidetes", "Proteobacteria")),
    2)
})

test_that("match_taxa works for taxa vectors of length 1", {
  expect_equal(
    match_taxa(
      c("Bacteroidetes; hgj", "Bacteria - Firmicutes", "kgj", "jasfdljh"),
      "Firmicutes"),
    c(NA, 1, NA, NA))
})

test_that("match_taxa works for empty lineage vectors", {
  expect_equal(
    match_taxa(
      character(),
      c("Firmicutes", "Bacteroidetes", "Proteobacteria")),
    integer())
})

test_that("match_taxa works for empty taxa vectors", {
  expect_equal(
    match_taxa(
      c("Bacteroidetes; hgj", "Bacteria - Firmicutes", "kgj", "jasfdljh"),
      character()),
    as.character(c(NA, NA, NA, NA)))
})

test_that("match_taxa gives a warning if more than one taxon is matched", {
  expect_warning(match_taxa(
    c("Bacteroidetes hfg", "Firmicutes"),
    c("Bacteroidetes", "Bacteroidetes hfg")))
})

test_that("first_true_index works for multiple TRUE values", {
  expect_equal(first_true_idx(c(FALSE, TRUE, TRUE)), 2)
})

test_that("first_true_index returns NA if there are no TRUE values", {
  expect_equal(first_true_idx(c(FALSE, FALSE, FALSE)), NA_integer_)
})
