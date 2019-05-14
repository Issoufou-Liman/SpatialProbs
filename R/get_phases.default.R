#' get_phase method for class numeric.
#'
#' This is the method on which get_phases () dispatches when x is simply the raw numeric vector.
#' see ?phases and ?get_phases
#' @author Issoufou Liman
#' @export
get_phases.default <- function(x, type = c("v_points", "peaks"), n_criticals = 1, steps = 2, ts_freq = 23, returned = c("ts_seasonal",
    "original")) {
    type <- match.arg(type)
    returned <- match.arg(returned)
    y <- phases(x, type = type, n_criticals = n_criticals, steps = n_criticals, ts_freq = ts_freq, returned = returned)
    y$phases
}
