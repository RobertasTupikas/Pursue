#' Run Selbal with BMDD Zero-Imputation
#'
#' This function runs BMDD-based zero-imputation on the OTU table,
#' then performs selbal.cv to find a balance associated with a response variable.
#'
#' @param otu_table A matrix or data frame of raw counts (samples x OTUs)
#' @param metadata A data frame of metadata (samples x variables)
#' @param variable_name A character string, the name of the response variable in metadata
#' @param n_fold Integer, number of cross-validation folds (default = 5)
#' @param n_iter Integer, number of cross-validation iterations (default = 10)
#' @param categorical Logical, whether to treat the response variable as categorical (factor). Default is TRUE.
#'
#' @return The selbal.cv result list, especially `global.balance`
#' @export
#'
#' @examples
#' selrez <- run_selbal_balance(otu_table, metadata, "brushing_event", categorical = TRUE)
#' selrez$global.balance
run_selbal_balance <- function(otu_table, metadata, variable_name,
                               n_fold = 5, n_iter = 10, categorical = TRUE) {
  if (!requireNamespace("BMDD", quietly = TRUE)) {
    stop("Please install the 'BMDD' package.")
  }
  if (!requireNamespace("selbal", quietly = TRUE)) {
    stop("Please install the 'selbal' package.")
  }

  if (!all(rownames(otu_table) %in% rownames(metadata))) {
    stop("Row names of OTU table must match sample names in metadata.")
  }

  library(BMDD)
  library(selbal)
  metadata <- metadata[rownames(otu_table), , drop = FALSE]

  bmdd.fit <- BMDD::bmdd(W = otu_table, type = 'count')
  bmmd_otu <- bmdd.fit$beta

  if (!(variable_name %in% colnames(metadata))) {
    stop("Variable not found in metadata.")
  }
  variable <- metadata[[variable_name]]
  if (categorical) variable <- as.factor(variable)

  selrez <- selbal::selbal.cv(
    x = bmmd_otu,
    y = variable,
    n.fold = n_fold,
    n.iter = n_iter
  )

  return(selrez)
}
