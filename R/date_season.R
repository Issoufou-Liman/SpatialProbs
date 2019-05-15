#' Get the begin and end dates of seasons in seasonal time series
#'
#' Getting the begin and end dates of seasons in seasonal time series such as rainfall data.
#' When the rainy season is bimodal, for instance, the function detect the 2 seasons first, then
#' return the begin and end dates of each.
#' @author Issoufou Liman
#' @param x A numeric vector or an object of class seasons
#' @author Issoufou Liman
#' @param ... further argument to be passed to the class constructor function `seasons ()`.
#' Ignored when `x` is of class `seasons`.
#' @return The beginning and the end of the season.
#' @details detecting the dates of the biginning and the end of the season is more likely to
#' yield meaningful result when `x` has been derived via `get_phases (x,..., returned = 'seasonal')` which
#' rescales the data to make more consistant ordering.
#' @seealso \code{\link[SpatialProbs]{phases}}, \code{\link[SpatialProbs]{get_phases}}
#' @examples
#' nam<- seq.Date(from = as.Date('2016-01-01'), to = as.Date ('2018-12-31'), by = 16)
#' dy11 <- c(1.40, 1.00, 1.50, 2.00, 5.00, 3.00, 1.00, 0.76, 2.00, 1.00, 3.50, 3.00, 1.50)
#' dy12 <- c(1.30, 1.10, 1.40, 2.01, 5.50, 2.80, 1.01, 1, 2.03, 1.09, 3.10, 3.00, 1.50)
#' dy1 <- c(dy11, dy12)
#' names(dy1) <- nam[1:length(dy1)]
#'
#' y <- get_phases (dy1, ts_freq = 12)
#' z <- lapply(y, date_season)
#' @export
#' @rdname date_season
date_season <- function(x, ...) {
    UseMethod("date_season")
}

#' @rdname date_season
#' @export
date_season.default <- function(x, ...) {
  y <- seasons(x, ...)
  y <- y$season_dates
  y <- as.Date(y)
}

#' date_season method for class seasons
#'
#' This is the method on which date_season () dispatches when x of class seasons.
#' see ?seasons
#' @rdname date_season
#' @export
date_season.seasons <- function(x, ...) {
  y <- x$season_dates
}
