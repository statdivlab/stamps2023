### In this lab we'll explore a dataset published by Wirbel et al. (2019).
### (https://www.nature.com/articles/s41591-019-0406-6)
### This is a meta-analysis of case-control studies, meaning that Wirbel
### et al. collected raw sequencing data from studies other researchers
### conducted and re-analyzed it (in this case, they also collected some
### new data of their own).

### Wirbel et al. published two pieces of data we'll focus on today:

# metadata giving demographics and other information about participants

# an mOTU table

### We'll look at a subset of all 849 mOTUs Wirbel et al. published

### We're most interested to compare microbial abundance in cases diagnosed
### with colorectal cancer to abundances in controls (without this diagnosis)

### First let's load libraries we'll need
library(tidyverse)
library(Matrix)
library(remotes)
# install radEmu (only if you don't already have it installed)
if (!("radEmu" %in% row.names(installed.packages()))) {
  install_github("https://github.com/statdivlab/radEmu")
}
library(radEmu)

metadata <-
  read_csv("https://raw.githubusercontent.com/statdivlab/stamps2023/main/labs/radEmu-lab/data/wirbel_et_al_metadata.csv")
head(metadata)

### Let's see how many observations we have among cases ("CRC") and
### controls ("CTR")
metadata %>%
  group_by(Group) %>%
  summarize(n = length(Group))

### metadata$Group tells us which participants are cases and which are controls

### We have data from studies in 5 different countries
### How much from each study, you ask? Let's find out!
metadata %>%
  group_by(Country) %>%
  summarize(n = length(Group))

### Let's see how many cases and controls were enrolled in each study as well
metadata %>%
  with(table(Country, Group))

### Now let's load the mOTU table
mOTU_table <-
  read_csv("https://raw.githubusercontent.com/statdivlab/stamps2023/main/labs/radEmu-lab/data/wirbel_et_al_mOTUs.csv")

# let's take a peek at the mOTU table
head(mOTU_table)
# how are the columns in this object named?
# (Looking at the metadata again may be informative)

### column "X1" of our mOTU table contains names of our mOTUs
### let's store them in a separate object
mOTU_names <- mOTU_table$mOTU


### now we'll clean up the mOTU table a bit
mOTU_table <- mOTU_table %>%
  dplyr::select(-1) %>% #remove first column (mOTU names)
  as.matrix() %>% # make this whole deal a matrix
  t() %>% # transpose so that rows are samples and columns are mOTU names
  (Matrix::Matrix) #now make our transposed matrix into a ~fancy~ Matrix
# using the Matrix package

### we need to keep track of which columns of our transposed mOTU_table
### correspond to which mOTUs, so let's name them!
colnames(mOTU_table) <- mOTU_names

### we're not *quite* done with data cleaning and organization
### first, if we compare the dimensions of mOTU_table and metadata...
dim(mOTU_table)
dim(metadata)

### yike -- we have more measurements than metadata!
### in general, this would require further investigation...
### but you can take my word that we don't need the extra observations
### let's get rid of them!

rows_in_metadata <- sapply(rownames(mOTU_table),
                           function(x) x %in% metadata$Sample_ID)


mOTU_table <- mOTU_table[rows_in_metadata,]
### let's check our dimensions again now
dim(mOTU_table)
dim(metadata)

### success!
### let's check that rows of metadata and mOTU_table are in the same order
order_for_metadata <- sapply(rownames(mOTU_table),
                             function(x) which(x == metadata$Sample_ID))

# if rows are in same order, this plot should look like y = x
plot(order_for_metadata)

# ok nope - not in same order, so let's reorder:
metadata <- metadata[order_for_metadata,]

# now let's check our work and look at the ordering again
order_for_metadata <- sapply(rownames(mOTU_table),
                             function(x) which(x == metadata$Sample_ID))

plot(order_for_metadata)

# yay

### if we look at the first few rownames of mOTU_table and sample_IDs
### in metadata, we see the order also matches

metadata$Sample_ID %>% head
rownames(mOTU_table) %>% head

### one last data organizing / cleaning / beautifying step...
### Wirbel et al. published their mOTU table...
### after dividing each row (i.e., count data for each sample) by its total
### and then rounding anything below 1e-6 to zero
### ... please don't do this! Publish your count data as is -- other
### researchers can easily transform counts, but un-transforming is
### harder!

### in any case, Wirbel et al. also published library size (i.e., count totals)
### by sample, so we can recover most of the count data by multiplying each
### row by the corresponding library size (we still don't get back counts that
### were rounded to zero)

for(i in 1:nrow(mOTU_table)){
  mOTU_table[i,] <- mOTU_table[i,]*metadata$Library_Size[i]
}

### now we'll pull out some taxa we might be interested in to fit models on
### (we can also fit a model to all mOTUs, but this takes longer)
fuso <- sapply(mOTU_names,function(x) grepl("Fusobac",x,fixed = TRUE))
prevo <- sapply(mOTU_names,function(x) grepl("Prevotella ",x,fixed = TRUE))
porph <-  sapply(mOTU_names,function(x) grepl("Porphyromonas ",x,fixed = TRUE))
clostr <- sapply(mOTU_names, function(x) grepl("Clostridium",x, fixed = TRUE))
faecali <- sapply(mOTU_names, function(x) grepl("Faecalibact",x, fixed = TRUE))
eubact <- sapply(mOTU_names, function(x) grepl("Eubact",x, fixed = TRUE))

# take a look at "eubact" -- what does this look like? What information
# does it contain?

head(eubact)

### "Group" in metadata indicates whether each participant is case (diagnosed
### with colorectal cancer) or a control (no such diagnosis)
### Let's make this a factor
unique(metadata$Group)

### make "CTR" the reference category (i.e., the first category)
metadata$Group <- factor(metadata$Group,
                         levels = c("CTR","CRC"))

### we'll stick to the genera Eubacterium, Porphyromonas, Faecalibacteria,
### and Fusobacterium for now
# store names of the mOTUs in these genera
restricted_mOTU_names <- mOTU_names[eubact|porph|fuso|faecali]

# Figure out which columns of mOTU_table contain observations in
# Eubacterium, Porphyromonas, Faecalibacteria, or Fusobacterium
which_mOTU_names <- which(eubact|porph|fuso|faecali)

# We'll begin with data from a Chinese study Wirbel et al. analyzed
ch_study_obs <- which(metadata$Country %in% c("CHI"))

# let's fit a model!
# ... emuFit will talk at you a bit, but don't worry about it
ch_fit <-
  emuFit(formula = ~ Group, # this is a formula telling radEmu what predictors to
         # use in fitting a model
         # we are using Group -- i.e., an indicator for which
         # participants have CRC diagnoses and which do not
         # as our predictor
         data = metadata[ch_study_obs, #ch_study obs = we're
                         # only looking at rows
                         # containing observations
                         # from the chinese study
         ], # data
         # contains our predictor
         # data
         Y = as.matrix(mOTU_table[ch_study_obs,which_mOTU_names]),
         run_score_tests = FALSE) #which_mOTU_names = we're limiting the taxa
# we're looking at to Eubacterium, Faecalibacterium,
# Porphyromonas, and Fusobacterium
# (which_mOTU_names was constructed above)
### Ok we have estimates and confidence intervals for the group effect
### Let's take a look:

ch_df <- ch_fit %>%
  mutate(Genus = sapply(category,
                        function(x) ifelse(strsplit(x," ",fixed = TRUE)[[1]][1] %in% c("unknown","uncultured"),
                                           strsplit(x," ",fixed = TRUE)[[1]][2],
                                           strsplit(x," ",fixed = TRUE)[[1]][1]))) %>%
  mutate(category = factor(category,levels = category[order(Genus)])) 
ggplot(ch_df) + geom_point(aes(x = category, y = estimate,color = Genus),
                        size = .5) +
  geom_errorbar(aes(x = category, ymin = lower, ymax = upper,color = Genus),
                width = .25)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(-5,15))

## Interestingly, we estimate a meta-mOTU ""unknown Eubacterium [meta_mOTU_v2_7116]" 
## assigned to Eubacteria
## to have a much higher ratio of concentrations (comparing CRC group to control) 
## than is typical across the mOTUs we included in this analysis.

## The confidence interval for this effect does not include zero -- but(!!!) 
## the kind of confidence interval that is returned by default by emuFit 
## is not extremely reliable when counts are very skewed or sample size is 
## small-to-moderate. 

## To investigate further, let's run a robust score test, 
## which is more reliable in these settings (but also takes more time because 
## apparently we can't have nice things). For comparison, we'll also test 
## Fusobacterium Nucleatum, which we also estimate to have a much larger
## ratio of concentrations across groups than is typical among the taxa we 
## included in this model fit.
robust_score_tests_eubacterium_mOTU_and_f_nucleatum <- 
  emuFit(formula = ~ Group,
         data = metadata[ch_study_obs, ],
         B = ch_fit,
         test_kj = data.frame(k = c(2,2), j = c(3,36)), # j = 3 is F. nucleatum; 
         # j = 36 is the Eubacterium meta mOTU
         # (you can see this in output from the model we already fit)
         
         Y = as.matrix(mOTU_table[ch_study_obs,which_mOTU_names]))
## Let's take a look at the test output:
robust_score_tests_eubacterium_mOTU_and_f_nucleatum

## p = 0.31!!! 

## Does this make sense? Let's look at the Eubacterium mOTU counts alongside CRC group:
cbind(mOTU_table[ch_study_obs,"unknown Eubacterium [meta_mOTU_v2_7116]"],
      as.character(metadata$Group[ch_study_obs]))

## We only detect this meta-mOTU in a single sample in the Chinese study cohort!
## So, yes -- it makes sense that our test returns a relatively large p-value

## Now let's look at F. nucleatum:
cbind(mOTU_table[ch_study_obs,"Fusobacterium nucleatum s. nucleatum [ref_mOTU_v2_0777]"],
      as.character(metadata$Group[ch_study_obs]))

# this also makes sense given what we found --  F. nucleatum shows up in a sizeable minority of 
# CRC cases in relatively high counts, whereas we detect it (by Wirbel et al's standards)
# in only one control participant

# We could run robust score tests for every taxon, (there is code below), but we 
# will not now because it would too much lab time to run. 

# How do these results relate to differential abundance testing in other packages?
# Note that each software will be testing different hypotheses with different
# models and different estimators. 

# Run ALDEx2

# ALDEx2 is built to analyze RNA-seq data but can be used for 16S or shotgun data
# as well. It uses a Dirichlet-multinomial model to infer abundance from counts.
# It implements several tests, including a t test on abundance data transformed 
# using a centered log ratio transformation. 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!("ALDEx2" %in% row.names(installed.packages()))) {
  BiocManager::install("ALDEx2")
}
library(ALDEx2)

# More information can be found in vignettes and in the papers cited here: 
citation("ALDEx2")

# vector of covariate data 
covariate <- metadata$Group[ch_study_obs]

# round our data because ALDEx2 gets mad when the abundance inputs are not integers
# swap rows and columns (ALDEx2 requires samples as columns)
reads <- round(t(as.matrix(mOTU_table[ch_study_obs, which_mOTU_names])))

# run centered log ratio transformation on reads
aldex_clr <- aldex.clr(reads, covariate, mc.samples=500, denom="all", verbose=FALSE)

# calculate effect sizes and differences between conditions
# specifically, gets at the median difference between groups
aldex_effect <- aldex.effect(aldex_clr, CI=T, verbose=FALSE)

# lets check the dimension of our effect size data frame
dim(aldex_effect)
length(which_mOTU_names)

# notice that while we have 47 mOTUs we are interested in, we only have effect size
# estimates for 46. Let's look at why.
tail(rowSums(reads))

# This is because we have a single taxon (unknown Porphyromonas [meta_mOTU_v2_777])
# that is not observed in any of the samples we are considering. 
# ALDEx2 is unable to estimate parameters for a taxon that does not appear in 
# any samples, so this taxon is not included in analyses. 

# plot effect sizes and their confidence intervals
aldex_effect %>%
  # clean data (add taxon information, add Genus)
  mutate(category = restricted_mOTU_names[1:46],
         Genus = ch_df$Genus[1:46]) %>% 
  mutate(category = factor(category,levels = category[order(Genus)])) %>% 
  ggplot() + geom_point(aes(x = category, y = effect, color = Genus),
                      size = .5) +
  geom_errorbar(aes(x = category, ymin = effect.low, ymax = effect.high,
                    color = Genus),
                width = .25) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

# No specific taxa jump out here. However, let's continue by running a t-test to
# compare the two groups. 
aldex_res <- aldex(reads, covariate, mc.samples=500, test="t", effect=TRUE,
               include.sample.summary=FALSE, denom="all", verbose=FALSE)

# let's check out the p-values from ALDEx2
aldex_pvals <- aldex_res$we.ep

# therefore, append a NA to the end of p-value vector for this taxa
aldex_pvals <- c(aldex_pvals, NA)

# recall that we ran two robust score tests with radEmu and got the following:
robust_score_tests_eubacterium_mOTU_and_f_nucleatum

# check out the results for these two taxa from ALDEx2
data.frame(pval = aldex_pvals,
                       taxon = names(which_mOTU_names)) %>%
  filter(taxon %in% c("Fusobacterium nucleatum s. nucleatum [ref_mOTU_v2_0777]",
                      "unknown Eubacterium [meta_mOTU_v2_7116]"))

# Similar to radEmu, the p-value for differential abundance with respect to Group
# for the Fusobacterium is very small (< 0.05), and the p-value for the Eubacterium
# is larger. 

# Run ANCOM-BC

# ANCOM-BC tests for differential abundance by estimating unknown sampling 
# fractions (ratio of the expected absolute abundance of a taxon in a random 
# sample to its absolute abundance in a unit volume of the ecosystem), correcting 
# the bias induced by their differences through a log linear 
# regression model including the estimated sampling fraction as an offset term, 
# and identifying taxa that are differentially abundant with respect to the variable 
# of interest.

if (!("ANCOMBC" %in% row.names(installed.packages()))) {
  BiocManager::install("ANCOMBC")
}
library(ANCOMBC)
# we will also use phyloseq
library(phyloseq)

# More information can be found in vignettes and in the papers cited here: 
citation("ANCOMBC")

# ANCOM-BC requires certain data formats. One option is a phyloseq object
# turn our data into phyloseq object
sam_dat <- phyloseq::sample_data(metadata[ch_study_obs, ])
otus <- as.matrix(mOTU_table[ch_study_obs, which_mOTU_names])
row.names(otus) <- row.names(sam_dat)
otu_tab <- phyloseq::otu_table(otus, taxa_are_rows = FALSE)
phy_obj <- phyloseq::phyloseq(sam_dat, otu_tab)

# run ANCOM-BC
ancom_res <- ancombc(data = phy_obj, 
                     formula = "Group")

# which taxa do we have results for? 
res_names <- ancom_res$res$p_val$taxon
length(res_names)
# We only have 35 p-values. This means there are 12 taxa that are not included.

# Let's make a vector to tell us which taxa were included in this analysis
full_names <- colnames(otus)
taxa_used <- full_names %in% res_names

# Now let's look at some confidence intervals for the coefficient for the Group 
# variable
data.frame(category = restricted_mOTU_names[taxa_used],
                       lfc = ancom_res$res$lfc$GroupCRC,
                       lfc.low = ancom_res$res$lfc$GroupCRC - 1.96*ancom_res$res$se$GroupCRC,
                       lfc.high = ancom_res$res$lfc$GroupCRC + 1.96*ancom_res$res$se$GroupCRC,
                       Genus = ch_df$Genus[taxa_used]) %>%
  mutate(category = factor(category,levels = category[order(Genus)])) %>% 
  ggplot() + geom_point(aes(x = category, y = lfc, color = Genus),
                        size = .5) +
  geom_errorbar(aes(x = category, ymin = lfc.low, ymax = lfc.high,
                    color = Genus),
                width = .25) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

# Here we see quite a few taxa that appear to be differentially abundant! 

# Finally lets look at our p-values.
ancom_pvals <- rep(NA, ncol(otus))
ancom_pvals[taxa_used] <- ancom_res$res$p_val$GroupCRC
ancom_pvals

# again, compare to robust score tests with radEmu 
robust_score_tests_eubacterium_mOTU_and_f_nucleatum

# check out the results for these two taxa from ALDEx2
data.frame(pval = ancom_pvals,
           taxon = names(which_mOTU_names)) %>%
  filter(taxon %in% c("Fusobacterium nucleatum s. nucleatum [ref_mOTU_v2_0777]",
                      "unknown Eubacterium [meta_mOTU_v2_7116]"))

# Similar to radEmu, the p-value for differential abundance with respect to Group
# for the Fusobacterium is very small (< 0.05). Unfortunately we were not able
# to fit parameters and run a test for the Eubacterium, so we cannot compare to 
# radEmu. 

# Run DESeq2

# DESeq2 estimates the variance-mean dependence in count data from high-throughput 
# sequencing data (it was originally built for RNA-seq data) and tests for 
# differential expression based on a model using the negative binomial distribution.

if (!("DESeq2" %in% row.names(installed.packages()))) {
  BiocManager::install("DESeq2")
}
library(DESeq2)

# More information can be found in vignettes and in the papers cited here: 
citation("DESeq2")

# DESeq2 cannot run with abundances of 0
reads_pseudo <- reads
# add a pseudocount of 1 to each 0 in the abundance dataset 
reads_pseudo[reads == 0] <- 1
# create object to be used by DESeq2
dds <- DESeqDataSetFromMatrix(countData = reads_pseudo, 
                              colData = metadata[ch_study_obs, ], 
                              design = ~Group)
# run DESeq2
deseq_res <- DESeq(dds)
# obtain results
deseq_test_res <- results(deseq_res)

# Again, let's look at some confidence intervals for the coefficient for the Group 
# variable
data.frame(category = restricted_mOTU_names,
           lfc = deseq_test_res$log2FoldChange,
           lfc.low = deseq_test_res$log2FoldChange - 
             1.96*deseq_test_res$lfcSE,
           lfc.high = deseq_test_res$log2FoldChange + 
             1.96*deseq_test_res$lfcSE,
           Genus = ch_df$Genus) %>%
  mutate(category = factor(category,levels = category[order(Genus)])) %>% 
  ggplot() + geom_point(aes(x = category, y = lfc, color = Genus),
                        size = .5) +
  geom_errorbar(aes(x = category, ymin = lfc.low, ymax = lfc.high,
                    color = Genus),
                width = .25) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

# Similar to in ANCOM-BC, we see more extreme log fold differences for several 
# Fusobacterium taxa. Additionally here we see more extreme log fold differences
# for several Porphyromonas taxa. 

# Let's now check out our p-values.
deseq_pvals <- deseq_test_res$pvalue

# again, compare to robust score tests with radEmu 
robust_score_tests_eubacterium_mOTU_and_f_nucleatum

# check out the results for these two taxa from ALDEx2
data.frame(pval = deseq_pvals,
           taxon = names(which_mOTU_names)) %>%
  filter(taxon %in% c("Fusobacterium nucleatum s. nucleatum [ref_mOTU_v2_0777]",
                      "unknown Eubacterium [meta_mOTU_v2_7116]"))

# Once again, the p-value for differential abundance with respect to Group
# for the Fusobacterium is very small (< 0.05), and the p-value for the Eubacterium
# is larger and not significant at most alpha levels. 

# We can see that although each of these methods are estimating different parameters
# and using different models, when we zoom in to these two specific taxa, they 
# tend to give us similar results. 

# ------------------------------------------------------------------------------

# This is the end of this lab for today! Feel free to continue below to fit different models,
# and play around with radEmu, but several of these computations will take 
# significantly more time to run. 

# ------------------------------------------------------------------------------

## We can run score tests for all taxa by setting run_score_tests = TRUE
## This is the kind of thing you might want to let run overnight, though -- 
## it will probably take a bit. 

## For the time being, you would be justified in moving on to the next section 
## without running this. Uncomment each of the lines below to run this code. 

# test_all <- 
#   emuFit(formula = ~ Group, # this is a formula telling radEmu what predictors to
#          # use in fitting a model
#          # we are using Group -- i.e., an indicator for which
#          # participants have CRC diagnoses and which do not
#          # as our predictor
#          data = metadata[ch_study_obs, #ch_study obs = we're
#                          # only looking at rows
#                          # containing observations
#                          # from the chinese study
#          ],
#          B = ch_fit,# covariate_data
#          # contains our predictor
#          # data
#          Y = as.matrix(mOTU_table[ch_study_obs,which_mOTU_names]),
#          run_score_tests = TRUE)


### Let's look at a French study Wirbel et al. analyized
# -- do we see similar patterns?
fr_study_obs <- which(metadata$Country %in% c("FRA"))

fr_fit <-
  emuFit(~ Group,
         data = metadata[fr_study_obs, 
         ], 
         Y = as.matrix(mOTU_table[fr_study_obs,which_mOTU_names]),
         run_score_tests = FALSE)


## Let's compare to the results on the Chinese cohort:
pd <- position_dodge2(width = 0.5)
ch_fit$country <- "China"
fr_fit$country <- "France"
rbind(ch_fit,fr_fit) %>%
  mutate(Genus = sapply(category,
                        function(x) ifelse(strsplit(x," ",fixed = TRUE)[[1]][1] %in% c("unknown","uncultured"),
                                           strsplit(x," ",fixed = TRUE)[[1]][2],
                                           strsplit(x," ",fixed = TRUE)[[1]][1]))) %>%
  mutate(category = factor(category,levels = unique(category[order(Genus)]))) %>%
  ggplot() + 
  geom_point(aes(x = category,y = estimate,color = Genus),position = pd,size = 0.5) + 
  geom_errorbar(aes(x = category,ymin = lower, ymax = upper,color = Genus,linetype = country),position = pd) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, size = 6,hjust = 1))



### In the Chinese study, samples were taken from patients undergoing
### colonoscopy, with some patients providing samples before
### colonoscopy and some after -- maybe we should try to adjust for this
### by including timing relative to colonoscopy as a predictor in our model
### (Why?)

ch_fit_timing <- 
  emuFit(formula = ~ Group + Sampling_rel_to_colonoscopy, 
         data = metadata[ch_study_obs, 
         ], 
         Y = as.matrix(mOTU_table[ch_study_obs,which_mOTU_names]),
         run_score_tests = FALSE)

### Error :(
### "Error in soln_mat %*% Y_start[, j, drop = FALSE] : non-conformable arguments" -- seems like
### we might not be giving emuFit inputs that have the correct dimensions

metadata$Sampling_rel_to_colonoscopy[ch_study_obs]
timing_available <- !is.na(metadata$Sampling_rel_to_colonoscopy[ch_study_obs])
# a-ha: one participant is missing sampling time data. We'll exclude them for now.

ch_fit_timing <- 
  emuFit(formula = ~ Group + Sampling_rel_to_colonoscopy, 
         data = metadata[ch_study_obs, 
         ][timing_available,], 
         Y = as.matrix(mOTU_table[ch_study_obs,which_mOTU_names][timing_available,]),
         run_score_tests = FALSE)

### We can also fit a single model to data from multiple studies

fr_ch_study_obs <- which(metadata$Country %in% c("CHI","FRA"))

fr_ch_fit <-
  emuFit(~ Group + Country,
         data = metadata[fr_ch_study_obs,],
         Y = as.matrix(mOTU_table[fr_ch_study_obs,which_mOTU_names]),
         run_score_tests = FALSE)# tell us what's happening during optimization


fr_ch_fit %>%
  mutate(Genus = sapply(category,
                        function(x) ifelse(strsplit(x," ",fixed = TRUE)[[1]][1] %in% c("unknown","uncultured"),
                                           strsplit(x," ",fixed = TRUE)[[1]][2],
                                           strsplit(x," ",fixed = TRUE)[[1]][1]))) %>%
  filter(covariate == "GroupCRC") %>% # we only want CRC coefficients, not country coefficients
  mutate(category = factor(category,levels = unique(category[order(Genus)]))) %>%
  ggplot() + geom_point(aes(x = category, y = estimate,color = Genus),
                        size = .5) +
  geom_errorbar(aes(x = category, ymin = lower, ymax = upper,color = Genus),
                width = .25)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(-5,15))


### Notice that we didn't include Sampling_rel_to_colonoscopy in this regression
### (Optional exercise: fit the model above with Sampling_rel_to_colonoscopy
### added to the model)

## You can also fit the model to all data from all countries
all_country_fit <-
  emuFit(~ Group + Country,
         data = metadata,
         Y = as.matrix(mOTU_table[,which_mOTU_names]),
         run_score_tests = FALSE)


all_country_fit %>%
  mutate(Genus = sapply(category,
                        function(x) ifelse(strsplit(x," ",fixed = TRUE)[[1]][1] %in% c("unknown","uncultured"),
                                           strsplit(x," ",fixed = TRUE)[[1]][2],
                                           strsplit(x," ",fixed = TRUE)[[1]][1]))) %>%
  filter(covariate == "GroupCRC") %>% # we only want CRC coefficients, not country coefficients
  mutate(category = factor(category,levels = unique(category[order(Genus)]))) %>%
  ggplot() + geom_point(aes(x = category, y = estimate,color = Genus),
                        size = .5) +
  geom_errorbar(aes(x = category, ymin = lower, ymax = upper,color = Genus),
                width = .25)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_cartesian(ylim = c(-5,15))


