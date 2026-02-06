#'############################################################################
#' Two-segment generalized additive models for feature screening
#'
#' This function utilizes two-segment generalized additive models to identify
#' variables showing significant differences between young and old segments
#' based on age information.
#'
#' @import ggplot2 mgcv
#'
#' @param y A numeric vector variables of the feature exxpression.
#' @param Age A numeric vector variables of length \eqn{n} containing age
#' information for each sample.
#' @param Group A binary classification vector variables of length \eqn{n}
#' containing age segment for each sample.
#' @param k An integer specifying the basis dimension of spline smoothers
#' used in GAM fitting.
#' @param bs a two letter character string indicating the (penalized) smoothing
#' basis to use.
#' @param ci_level confidence interval.
#' @param n_grid the number of grids in the diagram.
#' @param show_lm Logical, set to \code{TRUE}, if you want to show the linear
#' line in the diagram.
#' @param main the title of the diagram.
#'
#' @examples
#' ## Not run:
#' ## First example
#' sim <- simData(N = 40,
#'                mean.pat = 2,
#'                sd.pat = 1,
#'                mean.diff = 2,
#'                sd.diff = 1,
#'                mean.bg = 0,
#'                sd.bg = 0.5,
#'                p_total = 2000,
#'                p_pattern_relevant = 40,
#'                p_diff_relevant = 40,
#'                percentage_overlap_variables = 0.5,
#'                seeds = 123)
#'
#' plot_feature_young_old_gam(sim$data[, 5], sim$Age, sim$Group,
#' k = 6, show_lm = TRUE)
#'
#' @export
plot_feature_young_old_gam <- function(y, Age, Group,
                                       k = 6,
                                       bs = "cs",
                                       ci_level = 0.95,
                                       n_grid = 200,
                                       show_lm = TRUE,
                                       main = "Feature") {


  df <- data.frame(
    y = as.numeric(y),
    Age = as.numeric(Age),
    Group = factor(as.character(Group), levels = c("Young", "Old"))
  )

  dY <- df[df$Group == "Young", , drop = FALSE]
  dO <- df[df$Group == "Old",   , drop = FALSE]

  mY <- mgcv::gam(y ~ s(Age, k = k, bs = bs), data = dY)
  mO <- mgcv::gam(y ~ s(Age, k = k, bs = bs), data = dO)

  gY <- seq(min(dY$Age), max(dY$Age), length.out = n_grid)
  gO <- seq(min(dO$Age), max(dO$Age), length.out = n_grid)

  pY <- predict(mY, newdata = data.frame(Age = gY), se.fit = TRUE)
  pO <- predict(mO, newdata = data.frame(Age = gO), se.fit = TRUE)

  z <- qnorm(1 - (1 - ci_level) / 2)

  fit_df <- rbind(
    data.frame(
      Age = gY, Group = factor("Young", levels = levels(df$Group)),
      fit = as.numeric(pY$fit),
      lwr = as.numeric(pY$fit - z * pY$se.fit),
      upr = as.numeric(pY$fit + z * pY$se.fit)
    ),
    data.frame(
      Age = gO, Group = factor("Old", levels = levels(df$Group)),
      fit = as.numeric(pO$fit),
      lwr = as.numeric(pO$fit - z * pO$se.fit),
      upr = as.numeric(pO$fit + z * pO$se.fit)
    )
  )

  if (show_lm) {
    lmY <- lm(y ~ Age, data = dY)
    lmO <- lm(y ~ Age, data = dO)
    fit_df$lm <- NA_real_
    fit_df$lm[fit_df$Group == "Young"] <- predict(lmY, newdata = data.frame(Age = fit_df$Age[fit_df$Group == "Young"]))
    fit_df$lm[fit_df$Group == "Old"]   <- predict(lmO, newdata = data.frame(Age = fit_df$Age[fit_df$Group == "Old"]))
  }

  ggplot(df, aes(x = Age, y = y, color = Group)) +
    geom_point(alpha = 0.6) +
    geom_line(data = fit_df, aes(x = Age, y = fit, color = Group),
              linewidth = 1.3, inherit.aes = FALSE) +
    geom_line(data = fit_df, aes(x = Age, y = lwr, color = Group),
              linetype = "dashed", linewidth = 0.8, alpha = 0.7, inherit.aes = FALSE) +
    geom_line(data = fit_df, aes(x = Age, y = upr, color = Group),
              linetype = "dashed", linewidth = 0.8, alpha = 0.7, inherit.aes = FALSE) +
    {if (show_lm) geom_line(data = fit_df, aes(x = Age, y = lm, color = Group),
                            linetype = "dotdash", linewidth = 1.0, alpha = 0.9,
                            inherit.aes = FALSE)} +

    theme_classic() +
    labs(title = main, y = "Expression")
}
