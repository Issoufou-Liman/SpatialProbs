#' get_phases method for class phases
#'
#' This is the method on which get_phases () dispatches when x is of class phases.
#' see ?phases
#' @author Issoufou Liman
#' @param x An object of class phases.
#' @export
get_phases.phases <- function(x) {
    x$phases
}
