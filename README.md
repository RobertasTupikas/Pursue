# songbird-for-R

Reimplementation of the songbird algorithm by Morton et. al. (2019) in R with additional functionalities.

Original publication:  
https://doi.org/10.1038/s41467-019-10656-5

Original repository for the python version:  
https://github.com/biocore/songbird

Songbird is a differential microbial abundance analysis method that uses a regression model to calculate beta coefficients for each OTU with one or more covariates in mind. These beta coefficients are then ranked from lowes to highest, with OTUs having the highest absolute (irrespective of sign) being considered differential, and those with the lowest absolute beta coefficients are considered stable. This ranking mainly serves as additional motivation to select a particular reference and to scout potentially differential OTUs. Significance testing is then done by calculating log ratios with a particular reference as the denominator and a t-test or regression model is used to test whether the log ratios differ between conditions (discrete variable) or are associated with a continuous variable. To remedy the effect of zeros on log ratio calculation, Bimodal Dirichlet distribution (BMDD) is introduced into the workflow of Songbird in order to model proportions and impute missing zeros.

Bimodal Dirichlet distribution:  
https://doi.org/10.1101/2025.05.08.652808

Repository:  
https://github.com/zhouhj1994/BMDD

You can install this package via:  
`devtools::install_github("RobertasTupikas/songbird-for-R")`

You will need to create a conda environment via Anaconda prompt. The environment has to contain the following:  
- tensorflow 2.14.0  
- tensorflow-probability 0.23.0  
- numpy 1.26.4  

You will have to activate the environment in R via: \
`use_condaenv("environment", required = TRUE)`

An example workflow is provided in:  
`songtest.R`

The data used is from the teeth brushing examples in the original publication by Morton et. al. (2019)

Here it's:  
`physeq.RDS`

## Workflow

1. **Run the regression model** to obtain beta coefficients. You may input either a `phyloseq` object or an OTU table + metadata. The model accepts standard training parameters. It is recommended to run the regression twice while changing the first OTU in the OTU table, since the first OTU is removed from analysis due to identifiability constraints. \

```r
run_songbird(
  physeq = my_phyloseq,        # or specify otu_table and metadata instead
  formula = "~ brushing_event",
  differential_prior = 1.0,
  learning_rate = 0.001,
  clipnorm = 10.0,
  batch_size = 5,
  seed = 0,
  num_test_samples = 5,
  epochs = 9000
)
```

2. **Make a rank plot** to visualize which OTUs have the largest absolute beta coefficients and are therefore most likely differential. You can also highlight OTUs of interest to check their position in the ranking. \

```r
plot_songbird_ranks(
  beta_clr = betas,
  coef_name = "brushing_eventbefore",
  highlight_otus = c("Haemophilus", "Actinomyces"),
  ylim_range = c(-3, 3),
  xlab = "OTUs",
  ylab = "CLR Beta Coefficient"
)
```

3. **Select a stable OTU** as the reference for log-ratio testing based on its stability (low absolute coefficient) and prevalence in samples. This OTU will be used as the denominator. \

```r
extract_middle_otus(
  beta_mat = betas,
  otu_table = otu_table, 
  covariate = "brushing_eventbefore"
)
```

4. **Model proportions using BMDD** to impute zeros and calculate log-ratios of each OTU to the chosen reference. The output is a log-ratio matrix with an added condition column for group labels. \

```r
logratio_from_bmdd(
  otu_table = otu_table,
  metadata = metadata,
  condition_col = "brushing_event",
  ref_otu = "41f67443ce8207be0c0a956c47823417"
)
```

5. **Run a statistical test** (either Wilcoxon or t-test) to compare log-ratios of each OTU between the specified conditions. Output includes test statistics, p-values, and FDR-adjusted values. \

```r
logratio_stat_test(
  logratio_df = ratios,
  condition_col = "brushing_event",
  test = "t"  # or "wilcox"
)
```
