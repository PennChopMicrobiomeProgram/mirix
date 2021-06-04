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
