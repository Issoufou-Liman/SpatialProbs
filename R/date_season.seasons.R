#' date_season method for class seasons
#'
#' This is the method on which date_season () dispatches when x of class seasons.
#' see ?seasons
#' @author Issoufou Liman
#' @export
date_season.seasons <- function(x) {
    y <- x$season_dates
}
