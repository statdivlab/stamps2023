### radEmu lab
### David Clausen
### July 24th, 2022


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
library(monotone)
install_github("https://github.com/statdivlab/radEmu")
library(radEmu)



metadata <-
  read_csv("https://raw.githubusercontent.com/statdivlab/stamps2022/main/Sunday-afternoon/labs/rademu_lab/wirbel_et_al_metadata.csv")
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
  read_csv("https://raw.githubusercontent.com/statdivlab/stamps2022/main/Sunday-afternoon/labs/rademu_lab/wirbel_et_al_mOTUs.csv")


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

# fun fact: I (David) use sapply a lot
# apparently very few other people do
# in any case, let me talk you through what this sapply is doing:
rows_in_metadata <- sapply( #sapply basically says "take X and do Y to each
  #element of X
  rownames(mOTU_table), # <- this is X, so
  #so far we've told sapply that we want it to do
  #something using the row names of mOTU_table
  #(the row names here tell us what sample each row
  #of the mOTU table contains observations on)
  #this next bit specifies what to do with the row names
  #namely, we're asking it to tell us if each row name
  #exists inside metadata$Sample_ID -- in other words,
  #to tell us which samples we have mOTU data for we
  #also have metadata for
  function(x) x %in% metadata$Sample_ID # <- this is Y
)


#so now that we've used everyone's favorite R function, sapply,
#we have an object "rows_in_metadata" that tells us which rows
#of mOTU_table correspond to rows in our metadata
#we'll use this to subset mOTU_table so we have data for the same
#participants in metadata and mOTU_table
mOTU_table <- mOTU_table[rows_in_metadata,]

#note: I'm using base R to do this -- you could also use tidyverse
#... which some of my colleagues might prefer


### in any case, let's check our dimensions again now
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
# again -- I'm using base R here. You could also use inner_join() to do the
# same thing (I am told)

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

mOTU_table[1:5, 1:5]

### now we'll pull out some taxa we might be interested in to fit models on
### (we can also fit a model to all mOTUs, but this takes longer)
#yay another sapply
fuso <- sapply(mOTU_names, #go through the names of our mOTUs
               #and for each name:
               
               function(x) grepl("Fusobac", #see if the name contains the string
                                 #"Fusobac"
                                 #if so, return TRUE; otherwise FALSE
                                 x,
                                 fixed = TRUE))
#now do the same thing for some other genera:
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

### we'll stick to the genera Eubacterium, Porphyromonas,
### and Fusobacterium for now
# store names of the mOTUs in these genera
restricted_mOTU_names <- mOTU_names[eubact|porph|fuso]

# Figure out which columns of mOTU_table contain observations in
# Eubacterium, Porphyromonas, or Fusobacterium
which_mOTU_names <- which(eubact|porph|fuso)

# Among just the mOTUs in Eubacterium, Porphyromonas, or Fusobacterium,
# flag those that are in Eubacterium
eubact_restr <- sapply(restricted_mOTU_names,
                       function(x) grepl("Eubact",x, fixed = TRUE))

### Ooh ok now we're going to define the constraint we'll use with radEmu
### For now, we'll make this a median over Eubacterium
### (so we require the median of effects we estimate for Eubacterium to be
### zero

# fun question time: how does using this constraint affect the interpretation of
# our estimates?

constraint_fn <- function(x){median(x[eubact_restr])}

### Great! We're ready to start fitting models!

# We'll begin with data from a Chinese study Wirbel et al. analyzed
ch_study_obs <- which(metadata$Country %in% c("CHI"))

# let's fit a model!
# ... emuFit will talk at you a bit, but don't worry about it
# ch_fit <-
#   emuFit(Y = as.matrix(mOTU_table[ch_study_obs,
#                                   which_mOTU_names]), 
#          formula = ~ Group, 
#          data = metadata[ch_study_obs, ]
#   )

# instead of fitting radEmu, let's try some other methods!
# start with ALDEx2
library(ALDEx2)

# run t test on data
conds <- metadata$Group[ch_study_obs]
# round our data because aldex gets mad when it doesn't get integers
reads <- round(t(as.matrix(mOTU_table[ch_study_obs, which_mOTU_names])))
x.all <- aldex(reads, conds, mc.samples=500, test="t", effect=TRUE,
               include.sample.summary=FALSE, denom="all", verbose=FALSE)

par(mfrow=c(1,2))
aldex.plot(x.all, type="MA", test="welch", xlab="Log-ratio abundance",
           ylab="Difference")
aldex.plot(x.all, type="MW", test="welch", xlab="Dispersion",
           ylab="Difference")

# next try ANCOM-BC
library(ANCOMBC)

# need to turn data into a phyloseq object

# covariate data for chinese study 
sam_dat <- phyloseq::sample_data(metadata[ch_study_obs, ])
otus <- as.matrix(mOTU_table[ch_study_obs, which_mOTU_names])
row.names(otus) <- row.names(sam_dat)
otu_tab <- phyloseq::otu_table(otus, taxa_are_rows = FALSE)
phy_obj <- phyloseq::phyloseq(sam_dat, otu_tab)

ancom_res <- ancombc(data = phy_obj, 
                     formula = "Group")
res_names <- ancom_res$res$p_val$taxon
full_names <- colnames(otus)
taxa_used <- which(full_names %in% res_names)
full_p_val <- rep(NA, ncol(otus))
full_p_val[taxa_used] <- ancom_res$res$p_val$GroupCRC

# compare results
res_df <- data.frame(aldex = c(x.all$we.ep, NA),
                     ancom = full_p_val)
ggplot(res_df, aes(x = aldex, y = ancom)) + 
  geom_point() + 
  xlim(c(0, 1)) + 
  ylim(c(0, 1)) + 
  geom_abline(intercept = 0, slope = 1, color = "red")
ggsave("")
