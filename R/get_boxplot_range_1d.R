#' Create data ranges based boxplots statistics and outliers
#'
#' Take a vector of data and returns the bounds corresponding the boxplot specifications.
#' @author Issoufou Liman
#' @inheritParams stats::quantile
#' @param x numeric. A vector from which the ranges are to be derived based on quantiles.
#' NA and NaN values are not allowed unless na.rm is TRUE.
#' @param split_IQR logical. Should the inter-quartile range be split at the median
#' to form different ranges? The default is TRUE
#' @details
#' The ranges are bounded by the lower outlier (if any), extreme of the lower whisker,
#' the lower ‘hinge’, the median, the upper ‘hinge’, the extreme of the upper whisker,
#' and the upper outlier (if any).
#' @return
#' A numeric vector containing the bounds labelled as out_min (the lower outlier),
#' ymin (the extreme of the lower whisker), lower (the lower ‘hinge’), middle (the median),
#' upper (the upper ‘hinge’), ymax (the extreme of the upper whisker), out_max (the upper outlier).
#' @examples
#' sample <- rnorm(100)
#' get_boxplot_range_1d (sample)
#' @importFrom stats IQR quantile
#' @export
get_boxplot_range_1d <- function(x, split_IQR = TRUE, na.rm = TRUE , type = 5, ...) {
  # Lower hinge to upper hinge using 25% and 75% percentile
  qnt <- quantile(x, probs=c(0.25, 0.5, 0.75), na.rm = na.rm, type = type, ...)
  # interquartile range
  H <- 1.5 * IQR(x, na.rm = na.rm, type = type)
  # H <- 1.58 * IQR(x, na.rm = na.rm, type = type) / sqrt(length(x))
  # y_min: lower whisker = smallest observation greater than or equal to lower hinge - 1.5 * IQR
  min_x <- min(x, na.rm = TRUE)
  y_min <- qnt[[1]] - H
  #y_min <- ifelse(y_min >= min_x, y_min, min_x)
  y_min <- min(x[which(x >= y_min)])

  # y_max: upper whisker = largest observation less than or equal to upper hinge + 1.5 * IQR
  max_x <- max(x, na.rm = TRUE)
  y_max <- qnt[[2]] + H
  #y_max <- ifelse(y_max <= max_x, y_max, max_x)
  y_max <- max(x[which(x <= y_max)])
  if(split_IQR){
    out <-  c('out_min' = min_x, # lowest values will be the lowest outliers
              'ymin' = y_min,
              'lower' = qnt[[1]],
              'middle' = qnt[[2]],
              'upper' = qnt[[3]],
              'ymax' = y_max,
              'out_max' = max_x) # highest value will be the highest outliers
  } else{
    out <- c('out_min' = min_x, # lowest values will be the lowest outliers
             'ymin' = y_min,
             'lower' = qnt[[1]],
             'upper' = qnt[[2]],
             'ymax' = y_max,
             'out_max' = max_x) # highest value will be the highest outliers
  }
  return(out)
}
