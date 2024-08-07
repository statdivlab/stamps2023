## happi 16s lab - STAMPS 2022 
## Sarah Teichman, Pauline Trinh, David Clausen, Amy Willis 
## Thursday 28 July 2022 

# --------------------- Introduction -------------------------

# In the `happi_lab`, we learned about using `happi` for modeling gene presence while
# accounting for differences in genome quality. 

# We can also use `happi` to model taxon presence while accounting for
# differences in genome quality! This lets us apply `happi` to
# amplicon data, looking at presence/absence of ASVs or OTUs.
# Happi news for those of you who use amplicon data! 

# The hierarchical modeling approach to incorporate genome quality information is flexible,
# so while it was built for pangenome analysis, it can also be used in this context. 
# For example, if you were interested in understanding whether *E. coli* genomes are more
# present (in any abundance) in the microbiomes of sea otters 
# compared to narwhals, then `happi` is well-suited
# to help you answer this question! 

# `happi` models the association between 
# a covariate (e.g. narwhal microbiome vs. sea otter microbiome) and taxon presence (*E. coli*) where the 
# **covariate** is the **primary predictor** of interest and 
# **gene presence** is the **outcome** of interest. 


# --------------------- Installing happi -------------------------

# Because we're working on the RStudio Server, we already have `happi` installed. However,
# if you're working at home, uncomment the following lines and run them in your console (not in
# this script).

# To install: 
if (!("remotes" %in% row.names(installed.packages()))) {
  install.packages("remotes")
}
# Install happi using remotes and build vignettes: 
if (!("happi" %in% row.names(installed.packages()))) {
  remotes::install_github("statdivlab/happi", build_vignettes = TRUE, dependencies = TRUE) 
}

# if you want to check that `happi` is installed on your system you can run: 
"happi" %in% rownames(installed.packages()) # This should return [1] TRUE if happi is installed


# --------------------- Load required packages for this lab -------------------------

# Let's make sure we have all the packages and data that we need to follow along in this lab 
library(happi)
library(tidyverse)
library(ggplot2)

# We are going to use data from the package `corncob` lab. To do this, we will install and
# load the `corncob` package, along with the `phyloseq` package to work with the soil data.

if (!("corncob" %in% row.names(installed.packages()))) {
  remotes::install_github("bryandmartin/corncob") 
}
library(corncob)
library(phyloseq)
data(soil_phylo_sample)
data(soil_phylo_otu)
data(soil_phylo_taxa)
soil_phylo <- phyloseq::phyloseq(phyloseq::sample_data(soil_phylo_sample),
                                 phyloseq::otu_table(soil_phylo_otu, taxa_are_rows = TRUE),
                                 phyloseq::tax_table(soil_phylo_taxa))

# --------------------------- Manipulate data -------------------------------------

# We can start by looking at a description with the data with 

soil_phylo

# We want to consider our data at the phylum level

# First, let's collapse the samples to the phylum level.

soil <- soil_phylo %>% 
  tax_glom("Phylum") 

# Let's examine the phyloseq object.

soil

# We now see that we have an OTU abundance table with 39 OTUs and 119 samples. We 
# can extract using `otu_table()`. Let's examine a small subset of our data in more 
# detail.

otu_table(soil)[1:3, 1:3]

# We can also see that we have 5 sample variables. We can extract this using 
# `sample_data()`. Let's again examine a small subset in more detail.

sample_data(soil)[1:3, ]

# We're going to be looking closely at the `Amdmt` variable. 

#  - `Amdmt`: Categorical variable representing one of three soil amendments; none (0), 
# biochar (1), and fresh biomass (2), respectively.

# Finally, we have a taxonomy table with 7 taxonomic ranks. 

tax_table(soil)[1:3, ]

# We want to consider the relationship between the presence or absence of each taxon
# and `Amdmt` (soil amendments), while accounting for sequencing depth (this is the 
# variable that will give us information about the quality of each sample).

# To start, we want to include presence and absence information for each taxon. We can
# do this by taking the otu table and checking whether each element in the table is 
# greater than 0 (if greater than 0 the taxon is considered as present in the sample, 
# if equal to 0 the taxon is considered absent).

pres_mat <- as.data.frame(otu_table(soil)) > 0

# We'll transpose our matrix so that the samples represent rows and turn it into a data frame
pres_mat <- as.data.frame(t(pres_mat))

# We'll replace our otu labels with phylum names.

names(pres_mat) <- tax_table(soil)[,2]

# We will add the sample names as a separate column 

pres_mat <- pres_mat %>%
  mutate(sample_name = row.names(pres_mat))

# Next, we'll add our covariate of interest, `Amdmt`. We are going to collapse the levels
# of soil amendments down to just two, 0 for no amendments and 1 for amendments. 

sample_data_lab <- row.names(sample_data(soil))

covar_df <- data.frame(sample_data(soil)) %>%
  mutate(sample_name = sample_data_lab,
         soil_add = as.factor(ifelse(sample_data(soil)$Amdmt == 0, 0, 1))) %>%
  select(sample_name, soil_add)

soil_df <- inner_join(pres_mat, covar_df, by = "sample_name")

# Finally, we need information about sampling depth. For each sample, we'll 
# record the number of total reads across all taxa in that sample. 

soil_df <- soil_df %>%
  mutate(seq_depth = colSums(as.data.frame(otu_table(soil))))

# --------------------------- Manipulate data -------------------------------------

# Let's consider phylum Tenericutes. We can check how many samples it appears in

sum(soil_df[, "Tenericutes"])

# It appears in 64 of the 119 samples (a little more than half).
# let's subset our data to only consider Tenericutes

tnct_df <- soil_df %>%
  select(Tenericutes, soil_add, seq_depth) %>%
  mutate("Type" = ifelse(soil_add == 1, "Unamended", "Amended"))

# Let's plot the relationship between these variables. 

ggplot(tnct_df, aes(x = seq_depth, y = Tenericutes, col = Type)) + 
  geom_jitter(height = 0.08, width = 0.00) + 
  theme_bw() + 
  labs(x = "Sequencing Depth",
       col = "Amendment",
       title = "Tenericutes presence by sequencing depth") + 
  theme(plot.title = element_text(hjust = 0.5))

# Here we can see that generally the samples in which Tenericutes is not present have 
# lower sequencing depth. We can also see more samples with no amendments (a value of 0) 
# in which Tenericutes is present. Let's model this relationship! 

# --------------------------- Modeling our data --------------------------------

# We're going to start by running a method to model presence of Tenericutes based on soil
# amendments without account for sequencing depth. We'll do this with a generalized linear
# model (GLM) and test for a difference with a Rao score test. 

ha <- glm(Tenericutes ~ soil_add, family="binomial", data = tnct_df)
h0 <- glm(Tenericutes ~ 1, family="binomial", data = tnct_df)
anova(ha, h0, test = "Rao")[2, "Pr(>Chi)"]

# We get a p-value of 0.011. This means that if our alpha level is at 0.05, we can reject
# the null hypothesis that the odds that Tenericutes is present is 
# equal when comparing amended and unamended soils. 

# However, we know that this approach does not account for sequencing depth. Let's test 
# again with happi, adding in our sequencing depth information. 

# To start, we need to build a design matrix that holds our covariate of interest. 

x_matrix <- model.matrix(~soil_add, data = tnct_df)

# the main function we'll be using is happi() which has various options that you can specify 
# We'll describe a few here but to see the full list of available options you can run: 
# ?happi

# The main options we'll use to run happi() include: 
#  - outcome: Your gene presence/absence variable that needs to be coded as 0 or 1
#  - covariate: the main predictor/covariate of interest formatted as a design matrix (as we 
# did up above with x_matrix)
#  - quality_var: the quality variable of interest. Currently happi only accepts 1 quality 
#     variable but we will expand this capability to multiple quality variables 
#  - max_iterations: the maximum number of E-M steps the algorithm will run for. A larger 
#     number of steps may mean a longer run time if the algorithm doesn't reach the change_threshold.  
#  - change_threshold: the maximum percent change of the likelihood for 5 iterations in a row 
#    for both the alternative and null models that would signal termination of the algorithm 
#  - epsilon: the fixed probability of observing a gene given that it shouldn't be present 
#    (aka contamination of the genome from other genomes)
#  - method: method for estimating f. Defaults to "splines" which fits a monotone spline with
#    df determined by argument spline_df; "isotone" for isotonic regression fit
#  - firth: uses Firth penalty (default is TRUE)
#  - spline_df: degrees of freedom (in addition to intercept) to use in monotone spline fit 
#    (default is 3)
#  - verbose: do you want all the intermediate values used by happi (TRUE) or just the final output (FALSE)

happi_results <- happi(outcome = tnct_df$Tenericutes, 
                       covariate = x_matrix, 
                       quality_var = tnct_df$seq_depth,
                       max_iterations = 100, 
                       change_threshold = 0.001, 
                       epsilon = 0, 
                       method = "splines", 
                       firth = T, 
                       spline_df = 3,
                       verbose = FALSE)

# Let's look  at our results! 

happi_results$pvalue_LRT

# Our p-value here is 0.034. In this case, at an alpha level of 0.05, we still able to 
# reject the null hypothesis. However, the strength of our evidence against the null hypothesis
# is slightly lower (a higher p-value) because we've now accounted for different sequencing depth 
# in our model. 

# Recall from our picture that this taxon was less detected in lower sequencing depth samples.
# happi gives a larger p-value (less evidence for a difference in odds of presence across groups) 
# because the pattern of detection or non-detection could be attributable to genome quality. 
# So we think happi is doing a good thing here -- by saving you from getting excited by a 
# signature attributable to sequencing depth, not biology... though in this case, if we
# set out to explore Tenericutes, we'd probably still be excited. 

# If we want to get the coefficient estimates from happi we can run the following: 

happi_results$beta

# and we see that our estimates are 1.46 for our intercept beta_0 and our estimate for beta_1 
# which corresponds to the presence of soil amendments is -1.08. Based on our results, we see 
# that Tenericutes is less likely to be present in soil that is amended, and this difference is 
# not significant at the 5% significance level (p = 0.036). 

# ----------------- But what if I want to consider all of my taxa? -------------------------

# You can parallelize your analyses! We'll do this with the package `parallel`. This is already installed
# on the RStudio Server. To install on your local RStudio, uncomment the following lines.

#if (!require("parallel", quietly = TRUE))
#  install.packages("parallel") # check that parallel is installed. If not then install. 

# load parallel package 
library(parallel) 

# we need to make a new design matrix to account for our full dataset 
x_matrix_all <- model.matrix(~soil_add, data = soil_df)

# this function will take in a taxon (denoted by the column of the dataframe soil_df) and run 
# happi on that taxon

run_happi_all <- function(colnum) {
  happi(outcome = unlist(soil_df[,colnum]), 
        covariate = x_matrix_all, 
        quality_var = soil_df$seq_depth,
        method = "splines", 
        firth = T, 
        spline_df = 3,
        max_iterations = 100, 
        change_threshold = 0.01, 
        epsilon = 0,
        verbose = FALSE)
}

# Additionally for quicker run-time we're going to keep the max_iterations 100 and a change 
# threshold of 0.01 

# To run happi in parallel: 
all_taxa_results <- mclapply(1:39, run_happi_all, mc.cores=6) 
# this should only take about a minute or two to finish running... 
# 1:39 denotes the columns of `soil_df` that include taxa presence information
# run_happi_all is the function we created above 
# and mc.cores allows us to specify how many cores we want to allocate to this computational 
# task 

# we were only able to optimize our model for some taxa
error_taxa <- unlist(lapply(all_taxa_results, function(x) str_detect(x[[1]], "Error")))
sum(error_taxa)
# we were not able to optimize the model for 6 of the 39 taxa

# make a vector of p-values for taxa for which we were able to optimize the model
pvalues <- lapply(all_taxa_results[!error_taxa], function(x) x$pvalue_LRT) %>% unlist
# make a matrix of coefficients for COG functions for which we were able to optimize the model
betas <- lapply(all_taxa_results[!error_taxa], function(x) x$beta) %>% do.call("rbind",.)

hyp_results <- tibble("taxon" = (colnames(soil_df)[1:39])[!error_taxa], 
                      pvalues, 
                      betas[,1],
                      betas[,2]) %>% 
  arrange(pvalues)
head(hyp_results)

# And here we see the variables
#  - taxon: contains taxon name (here the Phylum name)
#  - pvalues: the happi p-values 
#  - betas[, 1]: the estimate for beta0 that corresponds to our intercept 
#  - betas[, 2]: the estimate for beta1 that corresponds to our primary covariate (soil amendment)

# In this case, we see that only the taxon we considered earlier (Tenericutes) is signficiant
# at the 0.05 alpha level.

# From here, you can 
# 1. choose your favorite FDR or FWER method to address multiple comparisons! 
#     e.g., with the qvalue::qvalue (you can get this with install.packages("BiocManager") and 
#           BiocManager::install("qvalue")), then as follows
# all_taxa_hyp_results %>%
#   mutate(qvalue = qvalue::qvalue(pvalues, 0.05))
#     Note that this won't work well with only 39 phyla (qvalue needs more hypotheses
#    to well-estimate the proportion of true null hypotheses)... but it should work 
#    great at a finer taxonomic level when you have thousands of taxa. 
# 2. Look at taxa that are significant. Even if nothing is highly significant, you 
#    can look at what taxa were most differentially present across your environmental
#    types. 
# 3. (Optional) Get excited about your biology!

# --------------------- Sensitivity analyses using epsilon  -------------------------

# happi has a hyperparameter epsilon, which is the probability of observing a taxon 
# given that it shouldn't be present. This can usually be thought of as 
# "contamination" of our sample from other genetic information. 
# We are interested in trying different values of epsilon to assess the robustness of our results 

# For this set of 39 taxa let's try setting epsilon = 0.05 and compare with our results when 
# epsilon = 0

run_happi_all_e05 <- function(colnum) {
  happi(outcome = unlist(soil_df[,colnum]), 
        covariate = x_matrix_all, 
        quality_var = soil_df$seq_depth,
        method = "splines", 
        firth = T, 
        spline_df = 3,
        max_iterations = 100, 
        change_threshold = 0.01, 
        epsilon = 0.05,
        verbose = FALSE)
}
all_taxa_results_e05 <- mclapply(1:39, run_happi_all_e05, mc.cores=6) 

# we were able to optimize our model for all taxa. Great!

# make a vector of p-values for taxa for which we were able to optimize the model
pvalues_e05 <- lapply(all_taxa_results_e05, function(x) x$pvalue_LRT) %>% unlist
# make a matrix of coefficients for COG functions for which we were able to optimize the model
betas_e05 <- lapply(all_taxa_results_e05, function(x) x$beta) %>% do.call("rbind",.)

hyp_results_e05 <- tibble("taxon" = colnames(soil_df)[1:39], 
                          pvalues_e05, 
                          betas_e05[,1],
                          betas_e05[,2]) %>% 
  arrange(pvalues_e05)
head(hyp_results_e05)

# To compare results let's merge both sets of results and look at the differences between the 
# p-values, beta0, and beta1 when using epsilon = 0 and 0.05 
hyp_results_comparison <- full_join(hyp_results, hyp_results_e05, by = "taxon") %>% # merge with previous results when epsilon = 0 
  mutate(pvalue_diff = round(pvalues_e05-pvalues, 3), # round differences to 3 decimal places 
         beta0_diff = round(as.numeric(`betas_e05[, 1]`) - as.numeric(`betas[, 1]`), 3), 
         beta1_diff = round(as.numeric(`betas_e05[, 2]`) - as.numeric(`betas[, 2]`), 3)) %>%
  arrange(pvalues)

# We can visualize the differences in p-values, beta0, and beta1 between happi results when 
# using epsilon = 0 and epsilon = 0.05 

hyp_results_comparison %>%
  ggplot(aes(x = pvalues_e05, y = pvalues)) +
  geom_point() +
  geom_abline() +
  xlab("happi pvalues with "~epsilon~"=0.05") +
  ylab("happi pvalues with "~epsilon~"=0") + 
  xlim(c(0,1)) + 
  ylim(c(0,1))

# Here we can see that many of our p-values are very similar when we change epsilon. However
# there is one taxon where the p-value changes by more than 0.3. Let's find which taxon this 
# change in p-value is coming from. 

hyp_results_comparison %>% 
  filter(abs(pvalue_diff) > 0.3) %>%
  select(taxon, contains("p"))

# These results are coming from OD1. Let's look closer at this Phylum! 

hyp_results_comparison %>% 
  filter(taxon == "OD1") %>%
  select(taxon, contains("p"))

# We can see that the p-value when epsilon = 0 is 1 and when epsilon = 0.05 the p-value is
# 0.133. Let's look closer at this taxon.

ggplot(soil_df, aes(x = seq_depth, y = `OD1`, col = soil_add)) + 
  geom_jitter(height = 0.08, width = 0.00) + 
  theme_bw() + 
  labs(x = "Sequencing Depth",
       col = "Amendment",
       title = "OD1 presence by sequencing depth") + 
  theme(plot.title = element_text(hjust = 0.5))
  
# This is a case where doing a sensitivity analysis is a great idea! It shows us that while
# the results for most taxa are robust to a change in epsilon, that is not the case for 
# OD1. 

# If we wanted to comprehensively understand the robustness of our results to varying degrees 
# of contamination we could continue to dial up or dial down values of epsilon (probability 
# between 0 and 1 of observing a gene given that it should be absent) and see what happens to 
# our results as epsilon increases! 

# --------------------- Final Notes -------------------------
# - `happi` can be used for testing for gene presence in pangenomes (this is its primary 
#    objective).
# - `happi` can be used when testing for gene presence between metagenomes of varying quality 
#    (e.g., sequencing depth)
# - `happi` can be used with multiple quality variables - let us know if you want this and 
#    we'll talk to you about it! 
# - If you have additional questions or issues with using `happi` please open an issue on our 
#    GitHub https://github.com/statdivlab/happi

# If you find happi useful for your work, please consider citing
# Pauline Trinh, David S. Clausen, and Amy D. Willis. happi: a hierarchical approach to 
# pangenomics inference. bioRxiv e-prints, April 2022.

