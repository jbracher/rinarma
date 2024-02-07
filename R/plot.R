#' Plotting functions for INARMA(1, 1) models
#'
#' Currently three plot types are available: model fit (`type = "fit"`),
#' residuals (`type = "residuals"`) and residual autocorrelation (`type = "racf"`).
#'
#' Plots are done in the base R plotting framework.
#' @export
#' @param model an object of class `"inarma"` to plot.
#' @param type the type of plot. One out of `"fit", "residuals", "racf"`.
#' @param interval_level for `type == "fit"`: the level for uncertainty intervals.
#' @param residual_type type of residuals to be used; one of `"raw"` or `"Pearson"`.
#' @param ylim the y-limit; only used for `type == "racf"`.
#' @param legend should a legend be shown? Only used for `type == "fitted`.
#' @param lag.max maximum lag to display if `type == "racf"`.

plot.inarma <- function(model, type = c("fit", "residuals", "racf"),
                        interval_level = 0.8, residual_type = "Pearson",
                        ylim = NULL, legend = TRUE, lag.max = 10){

  if(model$fitting_method != "maximum_likelihood"){
    stop("Plotting functions require that the model has been fitted using fit_inarma (rather than fit_inarma_moments).")
  }

  if(type == "fit"){
    plot_inarma_fitted(model = model, interval_level = interval_level,
                       ylim = ylim, legend = legend)
  }

  if(type == "residuals"){
    plot_inarma_residuals(model = model, type = residual_type)
  }

  if(type == "racf"){
    plot_inarma_racf(model, type = residual_type, lag.max = lag.max)
  }
}

#' Plot fitted values of INARMA model
#'
#' @param model an object of class `"inarma"` to plot.
#' @param interval_level a level for the uncertainty intervals (shaded area) displayed.
#' @param ylim a y-limit.
#' @param legend should a legend be shown?
plot_inarma_fitted <- function(model, interval_level = 0.8, ylim = NULL, legend = TRUE){
  # determine ylim:
  if(is.null(interval_level)){
    ylim <- c(0, 1.2*max(model$observed))
  }
  # empty plot:
  plot(model$observed, xlab = "t", ylab = "observed", ylim = ylim, col = "white")
  # define colors:
  col_fitted <- "red"
  col_intervals <- rgb(1, 0, 0, alpha = 0.5)

  # compute cdfs:
  cdfs <- apply(model$lik_distr, 1, cumsum)

  # compute lower and upper interval ends:
  level_lower <- (1 - interval_level)/2
  lower <- apply(cdfs, MARGIN = 2, function(vect) max(min(which(vect > level_lower)) - 1, 0))

  level_upper <- 1 - (1 - interval_level)/2
  upper <- apply(cdfs, MARGIN = 2, function(vect) max(max(which(vect < level_upper)), 0))

  # add to plot:
  t <- seq_along(model$observed)
  polygon(c(t, rev(t)), c(lower, rev(upper)), col = col_intervals, border = NA)

  # add fitted values:
  lines(model$fitted_values, col = col_fitted)

  # add observations:
  points(model$observed, pch = 20, cex = 0.75)

  # add legend if requested:
  if(legend){
    legend("top", pch = c(20, NA, 15), lty = c(NA, "solid", NA),
           col = c("black", col_fitted, col_intervals),
           legend = c("observed", "fitted", paste0("intervals (", interval_level, ")")),
           bty = "n", cex = 0.75)
  }
}

#' Plot residuals of an INARMA(1, 1) model
#'
#' @param model an object of class `"inarma"`.
#' @param type the type of residuals to be plotted; one of `"raw"` or `"Pearson"`.
plot_inarma_residuals <- function(model, type = c("raw", "Pearson"), ylim = NULL){
  if(type == "raw"){
    residuals_to_plot <- model$observed - model$fitted_values
    type <- "Raw" # for plot
  }
  if(type == "Pearson"){
    residuals_to_plot <- model$pearson_residuals
  }

  plot(residuals_to_plot,
       xlab = "t", ylab = "residual",
       main = paste(type, "residuals"),
       type = "h",
       ylim = ylim)
}

#' Plot residual autocorrelation function of an INARMA(1, 1) model
#'
#' @param model an object of class `"inarma"`
#' @param type the type of residuals to be used; one of `"raw"` or `"Pearson"`.
plot_inarma_racf <- function(model, type = "Pearson", lag.max = 10){
  if(type == "raw"){
    residuals_to_plot <- model$observed - model$fitted_values
  }
  if(type == "Pearson"){
    residuals_to_plot <- model$pearson_residuals
  }

  acf(residuals_to_plot, main = paste("ACF of", type, "residuals"), lag.max = lag.max)
}