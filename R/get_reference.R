#' Extract middle-ranked OTUs based on coefficient absolute values
#'
#' This function identifies the 10 OTUs whose coefficients (for a chosen covariate)
#' are closest to zero (i.e., least differential).
#' It returns their coefficients and also shows their minimum and maximum counts,
#' mean count, and prevalence across all samples using the provided OTU table.
#'
#' @param beta_mat A coefficient matrix from run_songbird(), with taxa as rows and covariates as columns.
#' @param otu_table A count table (samples x taxa) as a matrix or data.frame.
#'        Column names must match the taxa names in beta_mat.
#' @param covariate Character string, the column name in beta_mat to assess (e.g., "brushing_eventbefore").
#' @return A data.frame with OTU ID, coefficient, min count, max count, mean count, and prevalence.
#' @export
extract_middle_otus <- function(beta_mat, otu_table, covariate) {
  # Check inputs
  if (!covariate %in% colnames(beta_mat)) {
    stop("Covariate not found in beta matrix.")
  }

  if (!all(rownames(beta_mat) %in% colnames(otu_table))) {
    stop("OTU table is missing some OTUs present in the beta matrix.")
  }

  coefs <- beta_mat[, covariate]

  #order by absolute value, find 10 closest to 0
  middle_indices <- order(abs(coefs))[1:10]
  selected_otus <- rownames(beta_mat)[middle_indices]

  #counts for selected OTUs (columns = OTUs, rows = samples)
  counts <- otu_table[, selected_otus, drop = FALSE]

  #statistics
  min_counts   <- apply(counts, 2, min)
  max_counts   <- apply(counts, 2, max)
  mean_counts  <- colMeans(counts)
  prevalence   <- colMeans(counts > 0)

  result <- data.frame(
    OTU         = selected_otus,
    Coefficient = coefs[middle_indices],
    Min_Count   = min_counts,
    Max_Count   = max_counts,
    Mean_Count  = mean_counts,
    Prevalence  = prevalence,
    row.names   = NULL
  )

  return(result)
}
