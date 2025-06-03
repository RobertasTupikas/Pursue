
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

# Plot the results

plot_songbird_ranks(betasnoref,                                       # matrix of betas from run_songbird
                    "brushing_eventbefore",                           # covariate by which to rank OTUs
                    highlight_otus = c("51121722488d0c3da1388d1b117cd239",
                                       "6dee9ad8019d7a6a17aef3ffdc75a287"), # OTUs to highlight
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
