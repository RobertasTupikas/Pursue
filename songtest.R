
devtools::install_github("zhouhj1994/BMDD")
devtools::install_github("RobertasTupikas/songbird-for-R")

library(songbirdR)
library(BMDD)
library(tensorflow)
library(phyloseq)

# Load the conda environment with:
#            tensorflow 2.14.0
#            tensorflow-probability 0.23.0
#            numpy 1.26.4

use_condaenv(" ", required = TRUE)

# Loading phyloseq object with OTU table and metadata

readRDS("physeq.rds") -> my_phyloseq

otu_mat <- as.data.frame(otu_table(my_phyloseq))
otu_table <- t(as.matrix(otu_mat)) # OTU table has to have samples as rows and OTUs as columns
metadata <- as(sample_data(my_phyloseq), "data.frame")

# Running songbird on a raw OTU table and metadata

 betas <- run_songbird(
   otu_table = otu_table,        # OTU table with samples as rows and OTUs as columns
   metadata = metadata,          # metadata with samples as rows
formula = "~ brushing_event",# formula for metadata
differential_prior = 1.0,    # stddev of beta prior
learning_rate = 0.001,       # learning rate for optimizer
clipnorm = 10.0,             # gradient clipping threshold
batch_size = 5,              # size of training batches
seed = 0,                    # random seed for R and TF
num_test_samples = 5,        # number of samples held out
epochs = 9000               # number of training epochs
)

# Remove the first row (reference OTU) from betas

betasnoref <- betas[-1, ]

#Renaming OTUs by taxon dictionary from songbird data (not necessary if your OTUs are named properly)
name_map <- c(
  "51121722488d0c3da1388d1b117cd239" = "Haemophilus",
  "f6eb8537456dde7792eff60c8d1e797c" = "Haemophilus",
  "2cb16f19a154ae6de38df817b7ab8523" = "Haemophilus",
  "9e15e76a487328cc4fe9a0f2a48eb146" = "Haemophilus",
  "6dee9ad8019d7a6a17aef3ffdc75a287" = "Actinomyces",
  "76355427db027e5ec9f4404f906613a1" = "Actinomyces",
  "f58ef02fd3fe4f8f1c8ba13a3e205697" = "Actinomyces",
  "1b050f719c047f60bba533c3ba98bfd8" = "Actinomyces",
  "05c2d251996b17582af0cba08f500211" = "Actinomyces",
  "13f5aba07d66b4d4f6a20ccd480a3836" = "Actinomyces",
  "1dbf2096a837bcbda19f144fb527c499" = "Actinomyces",
  "8f0ba2d0c4a91f9303ab9e1b48d9f92b" = "Actinomyces"
)

rownames(betasnoref) <- ifelse(rownames(betasnoref) %in% names(name_map),
                               name_map[rownames(betasnoref)],
                               rownames(betasnoref))

# Plot the results

plot_songbird_ranks(betasnoref,                                       # matrix of betas from run_songbird
                    "brushing_eventbefore",                           # covariate by which to rank OTUs
                    highlight_otus = c("Haemophilus", "Actinomyces"), # OTUs to highlight (from the name_map)
                    ylim_range = c(-3, 3),                            # y-axis limits
                    xlab = "OTUs",                                    # x-axis label
                    ylab = "CLR Beta Coefficient")                    # y-axis label


extract_middle_otus(betas,                  # matrix of betas from run_songbird
                    otu_table,              # OTU table with samples as rows and OTUs as columns
                    "brushing_eventbefore") # covariate by which to rank OTUs


ratios <- logratio_from_bmdd(otu_table,                          # OTU table with samples as rows and OTUs as columns
                             metadata,                           # metadata with samples as rows
                             "brushing_event",                   # condition column in metadata
                             "41f67443ce8207be0c0a956c47823417") # Reference OTU for log-ratio calculation


result <- logratio_stat_test(ratios,            # log ratio output of logratio_from_bmdd()
                             "brushing_event",  # condition column in log-ratio table
                             test = "t")        # or "wilcox" for Wilcoxon test
