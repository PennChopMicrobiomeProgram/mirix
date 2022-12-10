weiss_sepsis <- subset(weiss2021_data, sample_id %in% "Sepsis.22441.C")
weiss_sepsis$proportion <- round(weiss_sepsis$proportion, 3)

weiss_healthy <- subset(weiss2021_data, sample_id %in% "Healthy.148")
weiss_healthy$proportion <- round(weiss_healthy$proportion, 3)

test_that("antibiotic_susceptibility_vancomycin works on Weiss examples", {
  expect_equal(
    antibiotic_susceptibility_vancomycin(weiss_sepsis$lineage),
    c("resistant", "susceptible", "resistant", "resistant", "resistant"))

  expect_equal(
    antibiotic_susceptibility_vancomycin(weiss_healthy$lineage),
    c("resistant", "susceptible", "resistant", "resistant", "susceptible",
      "susceptible", "resistant", "susceptible", "susceptible", "susceptible",
      "susceptible", "susceptible", "susceptible", "resistant", "susceptible",
      "susceptible", "resistant", "susceptible", "resistant", "susceptible",
      "susceptible", "resistant", "susceptible", "susceptible", "susceptible",
      "resistant"))

})

test_that("mirix_vancomycin works on Weiss examples", {
  expect_equal(
    mirix_vancomycin(weiss_sepsis$proportion, weiss_sepsis$lineage),
    -0.5670429, tolerance = 1e-5)

  expect_equal(
    mirix_vancomycin(weiss_healthy$proportion, weiss_healthy$lineage),
    -0.3609116, tolerance = 1e-5)
})

test_that("antibiotic_susceptibility_tetracycline works on Weiss examples", {
  expect_equal(
    antibiotic_susceptibility_tetracycline(weiss_sepsis$lineage),
    c("susceptible", "susceptible", "resistant", "susceptible", "susceptible"))

  expect_equal(
    antibiotic_susceptibility_tetracycline(weiss_healthy$lineage),
    c("susceptible", "susceptible", "resistant", "susceptible", "susceptible",
      "susceptible", "susceptible", "susceptible", "susceptible", "susceptible",
      "susceptible", "susceptible", "susceptible", "susceptible", "susceptible",
      "susceptible", "susceptible", "susceptible", "susceptible", "susceptible",
      "susceptible", "susceptible", "susceptible", "susceptible", "susceptible",
      "susceptible"))

})

test_that("mirix_doxycycline works on Weiss examples", {
  expect_equal(
    mirix_doxycycline(weiss_sepsis$proportion, weiss_sepsis$lineage),
    1.879934, tolerance = 1e-5)

  expect_equal(
    mirix_doxycycline(weiss_healthy$proportion, weiss_healthy$lineage),
    1.045228, tolerance = 1e-5)
})

test_that("antibiotic_susceptibility_penicillin works on Weiss examples", {
  expect_equal(
    antibiotic_susceptibility_penicillin(weiss_sepsis$lineage),
    c("susceptible", "susceptible", "resistant", "susceptible", "susceptible"))

  expect_equal(
    antibiotic_susceptibility_penicillin(weiss_healthy$lineage),
    c("susceptible", "susceptible", "resistant", "susceptible", "susceptible",
      "susceptible", "susceptible", "susceptible", "susceptible", "susceptible",
      "susceptible", "susceptible", "susceptible", "susceptible", "susceptible",
      "susceptible", "susceptible", "susceptible", "susceptible", "susceptible",
      "susceptible", "susceptible", "susceptible", "susceptible", "susceptible",
      "susceptible"))
})

test_that("mirix_amoxicillin works on Weiss examples", {
  expect_equal(
    mirix_amoxicillin(weiss_sepsis$proportion, weiss_sepsis$lineage),
    1.879934, tolerance = 1e-5)

  expect_equal(
    mirix_amoxicillin(weiss_healthy$proportion, weiss_healthy$lineage),
    1.045228, tolerance = 1e-5)
})

test_that("antibiotic_susceptibility_aminoglycoside works on Weiss examples", {
  expect_equal(
    antibiotic_susceptibility_aminoglycoside(weiss_sepsis$lineage),
    c("susceptible", "resistant", "resistant", "resistant", "resistant"))

  expect_equal(
    antibiotic_susceptibility_aminoglycoside(weiss_healthy$lineage),
    c("resistant", "resistant", "resistant", "resistant", "resistant",
      "resistant", "resistant", "resistant", "resistant", "resistant",
      "resistant", "resistant", "resistant", "resistant", "resistant",
      "resistant", "resistant", "resistant", "resistant", "resistant",
      NA, "susceptible", "resistant", "resistant", "resistant", "resistant"))
})

test_that("mirix_gentamicin works on Weiss examples", {
  expect_equal(
    mirix_gentamicin(weiss_sepsis$proportion, weiss_sepsis$lineage),
    0.52420, tolerance = 1e-5)
  expect_equal(
    mirix_gentamicin(weiss_healthy$proportion, weiss_healthy$lineage),
    -2.69373, tolerance = 1e-5)
})

test_that("mirix_metronidazole works on Weiss examples", {
  expect_equal(
    mirix_metronidazole(weiss_sepsis$proportion, weiss_sepsis$lineage),
    -1.788434, tolerance = 1e-5)
  expect_equal(
    mirix_metronidazole(weiss_healthy$proportion, weiss_healthy$lineage),
    2.089021, tolerance = 1e-5)
})

test_that("mirix_metronidazole and mirix_ciprofloxacin are inverses", {
  expect_equal(
    mirix_metronidazole(weiss_sepsis$proportion, weiss_sepsis$lineage),
    -mirix_ciprofloxacin(weiss_sepsis$proportion, weiss_sepsis$lineage))
  expect_equal(
    mirix_metronidazole(weiss_healthy$proportion, weiss_healthy$lineage),
    -mirix_ciprofloxacin(weiss_healthy$proportion, weiss_healthy$lineage))
})
