#' Plot ranked OTUs from Songbird output
#'
#' @param beta_clr A matrix of CLR-transformed beta coefficients (output from `run_songbird()`)
#' @param coef_name The name of the coefficient to rank by (e.g., "Depth")
#' @param highlight_otus A character vector of OTU names to highlight (optional)
#' @param ylim_range Optional y-axis limits as a numeric vector of length 2 (default = c(-0.03, 0.03))
#' @param xlab Label for the x-axis (default = "OTUs")
#' @param ylab Label for the y-axis (default = "CLR Beta Coefficient")
#'
#' @return A ggplot2 object showing ranked OTUs
#' @export
#'
#' @examples
#' result <- run_songbird(...)
#' plot_songbird_ranks(result$beta, coef_name = "Depth", highlight_otus = c("OTU1", "OTU2"))
plot_songbird_ranks <- function(beta_clr, coef_name, highlight_otus = NULL, ylim_range = c(-0.03, 0.03), xlab = "OTUs", ylab = "CLR Beta Coefficient") {
  library(ggplot2)
  library(dplyr)

  if (!(coef_name %in% colnames(beta_clr))) {
    stop("Coefficient name not found in beta matrix.")
  }

  coef_vals <- beta_clr[, coef_name]
  df <- data.frame(
    OTU = rownames(beta_clr),
    Beta = coef_vals
  ) %>%
    arrange(Beta) %>%
    mutate(Rank = row_number())

  df$Highlight <- ifelse(df$OTU %in% highlight_otus, df$OTU, NA)
  df$IsHighlighted <- !is.na(df$Highlight)

  base_plot <- ggplot(df, aes(x = Rank, y = Beta)) +
    geom_col(data = df[!df$IsHighlighted, ], fill = "gray90", width = 1) +
    geom_col(data = df[df$IsHighlighted, ], aes(fill = Highlight), width = 1.5) +
    scale_fill_manual(
      values = setNames(rainbow(length(highlight_otus)), highlight_otus),
      breaks = highlight_otus
    ) +
    coord_cartesian(ylim = ylim_range, expand = FALSE) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.title = element_blank(),
      legend.title = element_blank(),
      legend.key = element_blank(),
      legend.box.margin = margin(0, 0, 0, 20)
    ) +
    labs(
      x = xlab,
      y = ylab
    )

  return(base_plot)
}
