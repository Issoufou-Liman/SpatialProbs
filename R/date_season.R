#' S3 Generic for class seasons
#'
#' @param x A numeric vector or an object of class seasons
#' @author Issoufou Liman
#' @param ... further argument to be passed to the class constructor function seasons ().
#' @return The beginning and the end of the season
#' @details detecting the dates of the biginning and the end of the season is more likely to
#' yield meaningful result when x has been derived via get_phases (x,..., returned = 'seasonal') which
#' rescales the data to make more consistant ordering.
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
date_season <- function(x, ...) {
    UseMethod("date_season")
}
