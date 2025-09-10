# https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-9852/sdrf

library(tidyverse)
library(rvest)

covsets2 <- data.frame(ID = c("R235", "R236", "R237", "R239","R241","R242",
                          "R243", "R245","R247","R248","R249","R250"),
                 Treatment = factor(c("normal", "low", "mid", "high", "normal", "low", "mid", "high", "normal", "low", "mid", "high"),
                                    levels = c("normal", "low", "mid", "high")),
                 Batch = c(rep(c("T6", "T8", "T14_3"), each = 4)),
                 RIN = c(10, 10, 10, 9.8, 9.7, 10, 9.8, 9.6, 9.7, 9.6, 9.7, 9.3),
                 ExtractionDate = c(rep(c("2019-02-11", "2019-02-18"), times = c(8, 4))))
saveRDS(covsets2, file = "./extra-rna-seq-data/covsets2.rds")

counts <- read.table(file = "./extra-rna-seq-data/Abamectin_CountMatrix.txt", sep = "\t", header = T, dec = ".")
filter_id <- apply(counts, 1, mean) >= 1
counts2 <- counts[filter_id, ]
saveRDS(counts2, file = "./extra-rna-seq-data/counts2.rds")
