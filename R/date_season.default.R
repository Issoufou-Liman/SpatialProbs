#' date_season default method for class numeric.
#'
#' This is the method on which date_season () dispatches when x is simply the raw numeric vector.
#' see ?phases and ?get_phases
#' @author Issoufou Liman
#' @export
date_season.default <- function(x) {
    y <- seasons(x)
    y <- y$season_dates
    y <- as.Date(y)
}
