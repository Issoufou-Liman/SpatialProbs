#' Test whether a phases has a peak.
#'
#' Takes a vector and returns a list containing the indices of one or more peaks and
#' TRUE as test if it found a peak. if no peak is found, then has_peak () will returned
#' null as peak  and set the test to FALSE. Applied to phases slot of the phases object (see ?phases)
#' this function can help to determine whether a phase qualifies for a season or not.
#' Note that a season is a phases having a peak.
#' @author Issoufou Liman
#' @param x A numeric vector on which to test the existance of one or more peaks
#' @param npeaks An integer specifying the number of peaks to look for.
#' @param steps An interger specifying the number of points for a given
#' points to qualify for a peak.
#' @return a list containing of the indices of the peaks and logical test.
#' @details `has_peak ()` extends `get_falls ()` in the sense that it just internally applies
#' `get_falls ()` function to - x to locate the peaks and add the logical test.
#' @seealso \code{\link[SpatialProbs]{get_falls}}
#' @examples
#' ## data with complete cases
#' nam<- seq.Date(from = as.Date('2016-01-01'), to = as.Date ('2018-12-31'), by = 16)
#' dy11 <- c(1.40, 1.00, 1.50, 2.00, 5.00, 3.00, 1.00, 0.76, 2.00, 1.00, 3.50, 3.00, 1.50)
#' dy12 <- c(1.30, 1.10, 1.40, 2.01, 5.50, 2.80, 1.01, 1, 2.03, 1.09, 3.10, 3.00, 1.50)
#' dy1 <- c(dy11, dy12)
#' names(dy1) <- nam[1:length(dy1)]
#'
#' default_par <- par()
#' layout(rbind(c(1, 1), c(2, 3)))
#' par(mar = c(2, 2, 1, 1))
#' plot (dy1, type = 'o')
#' y <- phases (dy1, ts_freq = 12)
#' y <- y$phases
#' lapply (y, has_peak)
#' lapply (y, plot, type = 'o')
#'
#'
#' ## data with manageable number of NAs
#' dy2 <- dy1
#' dy2[23:length(dy2)] <- NA
#'
#' plot (dy2, type = 'o')
#' y <- phases (dy2, ts_freq = 12)
#' y <- y$phases
#' lapply (y, has_peak)
#' lapply (y, plot, type = 'o')
#'
#'
#' ## data with unmanageable number of NAs. This does not make sense as a season!
#' dy3 <- dy1
#' dy3[13:20] <- NA
#'
#' plot (dy3, type = 'o')
#' y <- phases (dy3, ts_freq = 12)
#' y <- y$phases
#' lapply (y, has_peak)
#' lapply (y, plot, type = 'o')
#' par(default_par)
#' @export
has_peak <- function(x, npeaks = 1, steps = 2) {
    peaks <- get_falls((-x), n_v_shape = npeaks, steps = steps)
    if (is.null(peaks) || length(peaks) == 0) {
        return(list(peak = peaks, test = FALSE))
    } else {
        return(list(peak = peaks, test = TRUE))
    }
}
