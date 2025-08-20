#' Plot Log-Ratios by Group with Statistical Test
#'
#' Plots log-ratios grouped by a categorical variable and performs appropriate statistical test.
#'
#' @param logratios Data frame with 'logratio' column and rownames as sample IDs
#' @param metadata Data frame with rownames as sample IDs and at least one grouping variable
#' @param group_var Character, the name of the grouping variable in metadata
#'
#' @return A ggplot object
#' @export
#'
#' @import ggplot2
#'

plot_logratios_by_group <- function(logratios, metadata, group_var) {
  if (!"logratio" %in% colnames(logratios)) stop("logratios must contain a column named 'logratio'")
  if (!group_var %in% colnames(metadata)) stop(paste("Grouping variable", group_var, "not found in metadata"))

  common_samples <- intersect(rownames(logratios), rownames(metadata))
  if (length(common_samples) == 0) stop("No overlapping samples between logratios and metadata")

  df <- data.frame(
    logratio = logratios[common_samples, "logratio"],
    group = as.factor(metadata[common_samples, group_var])
  )

  group_levels <- levels(df$group)
  normality <- sapply(group_levels, function(g) {
    group_data <- df$logratio[df$group == g]
    if (length(group_data) >= 3) {
      shapiro.test(group_data)$p.value > 0.05
    } else {
      FALSE
    }
  })

  n_groups <- length(group_levels)
  test_result <- NULL
  test_label <- ""

  if (all(normality)) {
    if (n_groups == 2) {
      test_result <- t.test(logratio ~ group, data = df)
      test_label <- "t-test"
    } else {
      test_result <- aov(logratio ~ group, data = df)
      p <- summary(test_result)[[1]][["Pr(>F)"]][1]
      test_result <- list(p.value = p)
      test_label <- "ANOVA"
    }
  } else {
    if (n_groups == 2) {
      test_result <- wilcox.test(logratio ~ group, data = df)
      test_label <- "Wilcoxon"
    } else {
      test_result <- kruskal.test(logratio ~ group, data = df)
      test_label <- "Kruskal-Wallis"
    }
  }

  pval <- signif(test_result$p.value, 3)

  library(ggplot2)
  ggplot(df, aes(x = group, y = logratio)) +
    geom_jitter(width = 0.2, height = 0, alpha = 0.6, size = 2) +
    geom_boxplot(alpha = 0.2, outlier.shape = NA) +
    labs(
      x = group_var,
      y = "Log-Ratio",
      subtitle = paste(test_label, "p =", pval)
    ) +
    theme_minimal()
}
