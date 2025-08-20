#' Plot Log-Ratios vs Continuous Variable with Regression Line
#'
#' @param logratios Data frame with 'logratio' column and rownames as sample IDs
#' @param metadata Data frame with rownames as sample IDs and at least one numeric variable
#' @param continuous_var Character, the name of the continuous variable in metadata
#'
#' @return A ggplot object
#' @export
#'
#' @import ggplot2
#'
plot_logratios_by_continuous <- function(logratios, metadata, continuous_var) {
  if (!"logratio" %in% colnames(logratios)) stop("logratios must contain a column named 'logratio'")
  if (!continuous_var %in% colnames(metadata)) stop(paste("Variable", continuous_var, "not found in metadata"))

  common_samples <- intersect(rownames(logratios), rownames(metadata))
  if (length(common_samples) == 0) stop("No overlapping samples between logratios and metadata")

  df <- data.frame(
    logratio = logratios[common_samples, "logratio"],
    x = metadata[common_samples, continuous_var]
  )

  if (!is.numeric(df$x)) stop("Provided variable must be numeric")

  model <- lm(logratio ~ x, data = df)
  residuals_normal <- shapiro.test(residuals(model))$p.value > 0.05

  cor_test <- if (residuals_normal) {
    cor.test(df$logratio, df$x, method = "pearson")
  } else {
    cor.test(df$logratio, df$x, method = "spearman")
  }

  method_label <- if (residuals_normal) "Pearson" else "Spearman"
  pval <- signif(cor_test$p.value, 3)
  r_val <- signif(cor_test$estimate, 3)

  ggplot(df, aes(x = x, y = logratio)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "steelblue") +
    labs(
      x = continuous_var,
      y = "Log-Ratio",
      title = paste("Log-Ratios vs", continuous_var),
      subtitle = paste(method_label, "r =", r_val, ", p =", pval)
    ) +
    theme_minimal()
}
