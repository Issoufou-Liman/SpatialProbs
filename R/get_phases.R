#' S3 Generic for class phases
#'
#' get_phases method for getting the phasal data according to the specified type arguments
#' (peaks or v_points). By default, the phases are derived according to v_points. see ?phases,
#' the constructor function for the class phases.

#' @param x Either a numeric vector or an object of class phases.
#' @param ... further arguments passed to the class constructor function phases.
#' Ignored when `x` is of class `phases`
#' @seealso phases
#' @return A list containing the phasal data.
#' @examples
#' nam<- seq.Date(from = as.Date('2016-01-01'), to = as.Date ('2018-12-31'), by = 16)
#' dy11 <- c(1.40, 1.00, 1.50, 2.00, 5.00, 3.00, 1.00, 0.76, 2.00, 1.00, 3.50, 3.00, 1.50)
#' dy12 <- c(1.30, 1.10, 1.40, 2.01, 5.50, 2.80, 1.01, 1, 2.03, 1.09, 3.10, 3.00, 1.50)
#' dy1 <- c(dy11, dy12)
#' names(dy1) <- nam[1:length(dy1)]
#'
#' ## dispatches on the get_phases.default when x is the raw numeric vector
#' get_phases (dy1, ts_freq = 12)
#'
#' ## dispatches on the get_phases.phases when x is of class phases.
#' y <- phases(dy1, ts_freq = 12)
#' get_phases (y)
#' @export
get_phases <- function(x, ...) {
    UseMethod("get_phases")
}

#' @rdname get_phases
#' @export
get_phases.default <- function(x, ...) {
  y <- phases(x, ...)
  y$phases
}

#' @rdname get_phases
#' @export
get_phases.phases <- function(x, ...) {
  x$phases
}
