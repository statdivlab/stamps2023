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

# Let's look at diversity on the order level 
water_order <- water %>% speedyseq::tax_glom("Order")

# Now we can run DivNet
dv_water_order <- divnet(water_order, ncores = 4)

# Extract the Bray-Curtis dissimilarity matrix 
# (you can use whatever beta-diversity metric you care about)
# but ideally you are not using naive measures of beta-diversity! 
water_bc <- dv_water_order$`bray-curtis`

# Now we want to order the dissimilarity matrix by sample type
water_data <- sample_data(water_order)
water_data$SampleType
# Let's consider a binary covariate, freshwater or not 
freshwater <- ifelse(
  stringr::str_detect(water_data$SampleType, "Freshwater"),
  "freshwater",
  "non-freshwater")
freshwater

# It turns out that our samples are already ordered by freshwater
# If they weren't, we could order them with the following
ind_reord <- order(freshwater)
freshwater_reord <- freshwater[ind_reord]
water_bc_reord <- water_bc[ind_reord, ind_reord]

# Next we will extract the dissimilarity matrices that are within and 
# between groups
fresh_bc <- water_bc[freshwater == "freshwater", freshwater == "freshwater"]
non_fresh_bc <- water_bc[freshwater == "non-freshwater", 
                         freshwater == "non-freshwater"]
between_bc <- water_bc[freshwater == "freshwater", 
                       freshwater == "non-freshwater"]
# Let's check out the dimensions of each matrix
dim(fresh_bc)
dim(non_fresh_bc)
dim(between_bc)
# ok this seems right 

# Next we will pull out the upper triangular portion of each dissimilarity
# matrix
fresh_vec <- as.vector(fresh_bc[upper.tri(fresh_bc)])
non_fresh_vec <- as.vector(non_fresh_bc[upper.tri(non_fresh_bc)])  
between_vec <- as.vector(between_bc)  

# Now we will make a data frame using these dissimilarities
bc_df <- data.frame(
  dissim = c(fresh_vec, non_fresh_vec, between_vec),
  type = c(rep("Freshwater", length(fresh_vec)),
           rep("Non-freshwater", length(non_fresh_vec)),
           rep("Between type", length(between_vec))))

# Plot our dissimilarities
ggplot(data = bc_df, aes(x = type, y = dissim)) + 
  geom_jitter() + 
  ylim(c(0, 1)) + 
  labs(x = "", y = "Bray-Curtis dissimilarity")

# Here is looks like the dissimilarities between types are all 
# large and within types some are small and some are large
# Maybe this is because we collapsed sample types! 
# Let's make a similar plot but for the original levels of SampleType
water_data$SampleType
# freshwater
fresh_bc <- water_bc[water_data$SampleType == "Freshwater",
                     water_data$SampleType == "Freshwater"]
fresh_vec <- as.vector(fresh_bc[upper.tri(fresh_bc)])
# freshwater (creek)
creek_bc <- water_bc[water_data$SampleType == "Freshwater (creek)",
                     water_data$SampleType == "Freshwater (creek)"]
creek_vec <- as.vector(creek_bc[upper.tri(creek_bc)])
# ocean
ocean_bc <- water_bc[water_data$SampleType == "Ocean",
                     water_data$SampleType == "Ocean"]
ocean_vec <- as.vector(ocean_bc[upper.tri(ocean_bc)])
# sediment
sediment_bc <- water_bc[water_data$SampleType == "Sediment (estuary)",
                        water_data$SampleType == "Sediment (estuary)"]
sediment_vec <- as.vector(sediment_bc[upper.tri(sediment_bc)])
# between groups
between_mat <- sapply(water_data$SampleType, function(x) x != water_data$SampleType)
between_vec <- water_bc[upper.tri(water_bc) * between_mat == 1]
# data frame
sample_type_df <- data.frame(
  dissim = c(fresh_vec, creek_vec, ocean_vec, sediment_vec, between_vec),
  type = c(rep("Freshwater", length(fresh_vec)),
           rep("Creek", length(creek_vec)),
           rep("Ocean", length(ocean_vec)),
           rep("Sediment", length(sediment_vec)),
           rep("Between", length(between_vec))))

# Plot these results 
ggplot(sample_type_df, aes(x = type, y = dissim)) + 
  geom_jitter() + 
  ylim(c(0, 1)) + 
  labs(x = "", y = "Bray-Curtis dissimilarity")

# Finally, we can run a t-test to compare the mean of within versus
# between group distances
t.test(sample_type_df %>% filter(type != "Between") %>% pull(dissim),
       sample_type_df %>% filter(type == "Between") %>% pull(dissim))

# We have a 95% confidence interval for difference in means of dissimarilities
# between the two groups of (-0.625, -0.389) and a p-value < 0.05, so at
# an alpha level of 0.05 we can reject the hypothesis that the mean 
# dissimilarity within water type is equal to the mean dissimilarity
# between water type. 
