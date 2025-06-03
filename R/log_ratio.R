#' Compute log-ratios using BMDD against a reference OTU
#'
#' This function fits BMDD to the input OTU table, computes log-ratios
#' of each OTU against a specified reference OTU, and returns the log-ratio
#' matrix (excluding the reference OTU), with an added column for the condition.
#'
#' @param otu_table A matrix or data.frame with OTUs as rows and samples as columns.
#' @param metadata A data.frame with samples as rows and metadata columns.
#' @param condition_col Character string specifying the column in metadata with group labels.
#' @param ref_otu Character string specifying the name of the reference OTU.
#' @return A data.frame of log-ratios of all OTUs to a reference OTU (removed), with an added column for the condition.
#' @export
logratio_from_bmdd <- function(otu_table, metadata, condition_col, ref_otu) {
  if (!requireNamespace("BMDD", quietly = TRUE)) {
    stop("Package 'BMDD' must be installed.")
  }

  # Transpose OTU table to samples x taxa if needed
  if (nrow(otu_table) > ncol(otu_table)) {
    warning("Transposing OTU table: expected OTUs as columns and samples as rows.")
    otu_table <- t(otu_table)
  }

  otu_table <- as.data.frame(otu_table)
  if (!all(rownames(metadata) %in% rownames(otu_table))) {
    stop("Sample names in metadata must match those in OTU table.")
  }

  # Reorder to match
  metadata <- metadata[rownames(otu_table), , drop = FALSE]

  # Run BMDD
  bmdd_fit <- BMDD::bmdd(W = t(otu_table), type = 'count')  # BMDD expects OTU x sample
  beta <- bmdd_fit$beta  # OTUs x samples

  # Check reference OTU presence in posterior
  if (!ref_otu %in% rownames(beta)) {
    stop("Reference OTU not found in BMDD posterior output. It may have been dropped due to zero counts.")
  }

  # Normalize to proportions per sample
  prop_bmdd <- t(apply(beta, 2, function(x) x / sum(x)))  # samples x OTUs
  colnames(prop_bmdd) <- rownames(beta)  # OTU names
  rownames(prop_bmdd) <- colnames(beta)  # sample names

  if (any(prop_bmdd[, ref_otu] == 0)) {
    stop("Reference OTU has zero values in posterior proportions.")
  }

  # Correct log-ratio calculation (row-wise division)
  ref_vec <- prop_bmdd[, ref_otu]
  logratios <- log(sweep(prop_bmdd, 1, ref_vec, "/"))
  logratios <- logratios[, colnames(logratios) != ref_otu, drop = FALSE]

  # Add condition column
  logratios_df <- as.data.frame(logratios)
  logratios_df[[condition_col]] <- metadata[[condition_col]]

  return(logratios_df)
}
