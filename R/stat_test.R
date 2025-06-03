#' Perform statistical tests on log-ratio data
#'
#' This function takes a data.frame of log-ratios (samples x OTUs + condition column)
#' and performs Wilcoxon or t-tests on each OTU's log-ratio across the groups defined
#' by the specified condition column.
#'
#' @param logratio_df A data.frame produced by logratio_from_bmdd().
#' @param condition_col Character string of the name of the condition column.
#' @param test Character string: either "wilcox" or "t" to specify test type.
#' @return A data.frame with OTU, p-value, test statistic, and FDR-adjusted p-value.
#' @export
logratio_stat_test <- function(logratio_df, condition_col, test = "wilcox") {
  if (!condition_col %in% colnames(logratio_df)) {
    stop("Condition column not found in log-ratio table.")
  }
  
  otu_cols <- setdiff(colnames(logratio_df), condition_col)
  
  if (length(unique(logratio_df[[condition_col]])) != 2) {
    stop("Condition column must have exactly two groups for this test.")
  }
  
  group1_label <- unique(logratio_df[[condition_col]])[1]
  group2_label <- unique(logratio_df[[condition_col]])[2]
  
  results <- lapply(otu_cols, function(otu) {
    group1 <- logratio_df[logratio_df[[condition_col]] == group1_label, otu]
    group2 <- logratio_df[logratio_df[[condition_col]] == group2_label, otu]
    
    test_result <- if (test == "wilcox") {
      wilcox.test(group1, group2, exact = FALSE)
    } else if (test == "t") {
      t.test(group1, group2)
    } else {
      stop("Test must be either 'wilcox' or 't'.")
    }
    
    data.frame(
      OTU = otu,
      statistic = unname(test_result$statistic),
      p_value = test_result$p.value
    )
  })
  
  result_df <- do.call(rbind, results)
  result_df$fdr <- p.adjust(result_df$p_value, method = "fdr")
  
  return(result_df)
}
