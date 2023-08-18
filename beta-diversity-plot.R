## Scatterplot of beta diversity between and within groups
## Sarah Teichman 

# We will start by loading packages 

library(speedyseq) # or library(phyloseq) if this doesn't work
library(tidyverse)

# We will also use DivNet to get beta diversity estimates
if (!("remotes" %in% row.names(installed.packages()))) {
  install.packages("remotes")
}
if (!("DivNet" %in% row.names(installed.packages()))) {
  remotes::install_github("adw96/DivNet")
}
library(DivNet)

# We are going to use data distributed from `phyloseq` called GlobalPatterns

data("GlobalPatterns")
GlobalPatterns

# We're particularly interested in looking at only the water samples so we are
# going to need to subset our samples using the variable SampleType.

water <- GlobalPatterns %>%
  subset_samples(SampleType %in% c("Freshwater",
                                   "Freshwater (creek)",
                                   "Ocean",
                                   "Sediment (estuary)"))

