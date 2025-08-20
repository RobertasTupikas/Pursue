#' Compute Log‑Ratios (Geometric‑Mean, BMDD‑imputed)
#'
#' Performs Bayesian‑multiplicative zero imputation via \code{BMDD}
#' and computes per‑sample log‑ratios between the geometric means of
#' selected numerator and denominator OTUs.
#'
#' @param otu_table  Matrix or data.frame of raw counts
#'                   (samples in rows, taxa in columns).
#' @param numerator  Character vector of OTU IDs for the numerator.
#' @param denominator Character vector of OTU IDs for the denominator.
#'
#' @return data.frame with one column \code{logratio}; row names are sample IDs.
#' @export
#'
#' @examples
#' logratios <- compute_logratios_gm(bmmd_otu,
#'                                   c("OTU1", "OTU2"),
#'                                   c("OTU3"))
compute_logratios <- function(otu_table,
                                 numerator,
                                 denominator) {

  if (is.null(dim(otu_table)))
    stop("otu_table must be a matrix or data.frame")

  if (!requireNamespace("BMDD", quietly = TRUE))
    stop("Package 'BMDD' must be installed for zero imputation")

  bmdd.fit <- BMDD::bmdd(W = as.matrix(otu_table), type = "count")
  otu_imp  <- bmdd.fit$beta
  rownames(otu_imp) <- rownames(otu_table)
  colnames(otu_imp) <- colnames(otu_table)

  num_mat <- otu_imp[, colnames(otu_imp) %in% numerator,   drop = FALSE]
  den_mat <- otu_imp[, colnames(otu_imp) %in% denominator, drop = FALSE]

  if (ncol(num_mat) == 0)  stop("No valid numerator OTUs in table.")
  if (ncol(den_mat) == 0)  stop("No valid denominator OTUs in table.")

  num_log <- rowMeans(log(num_mat))
  den_log <- rowMeans(log(den_mat))
  logratio <- num_log - den_log

  data.frame(logratio = logratio,
             row.names = rownames(otu_imp))
}
