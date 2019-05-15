#' A constructor function for class phase
#'
#' Takes vector (typically seasonal data such as rainfall time series) makes phases or seasons,
#' and returns an object of class phases or seasons. phases, are here defined
#' as a cycle of some sort that can be distinguished in the time series.
#' The only important difference between phases and seasons is that seasons are phases having a peak.
#' @author Issoufou Liman
#' @param x a numeric vector which names are dates at which the data were acquired.
#' @param type character string. either `v_points` or `peaks` depending on how `x` should be
#' broken into phases. `v_points` (default) may make more sens for expression seasonal data.
#' @param n_criticals Numeric. The number of points to be considered. Default to 1.
#' @inheritParams get_falls
#' @param ts_freq The frequence of time series (see frequency argument in \code{\link[stats]{ts}}). Default to 23.
#' @param returned character string either `ts_seasonal` or `original.` if `ts_seasonal` (default),
#' the returned phases data will be the seasonal component of the decomposed time series (see \code{\link[stats]{decompose}}).
#' if 'original', the original data will be returned.
#' @return an object of class phases comprise a list of:
#' \itemize{
#'  \item{"phases"}{ list of phases}
#'  \item{"criticals"}{index of the points separating the phases}
#' }
#' @seealso \code{\link[SpatialProbs]{get_falls}}, \code{\link[SpatialProbs]{get_phases}}.
#' @examples
#' ## 2 years of a uni-modal pixel with complete without NAs
#' nam<- seq.Date(from = as.Date('2016-01-01'), to = as.Date ('2018-12-31'), by = 16)
#' dx11 <- c(1.30, 1.15,  1.50,  2.00,  2.01,  3.00, 3.20,  4.76,  3.50,  3.00,  2.40,  2.00,  1.50)
#' dx12 <-c(1.29, 1.1, 1.49, 1.99, 2, 3.1, 4.5, 4, 2.8, 2.5, 2.3, 1.6, 1.59)
#' dx1 <- c(dx11, dx12)
#' names(dx1) <- nam[1:length(dx1)]
#'
#' ## plotting the data
#' default_par <- par()
#' layout(rbind(c(1, 1), c(2, 3)))
#' par(mar = c(2, 2, 1, 1))
#' plot(dx1, type = 'o', main = 'raw data') # note 2 phases = 2 years and the starting point.
#'
#' ## returning the phases as seasonal data.
#' y1 <- phases(dx1, ts_freq = 12)
#' y1
#' lapply (X = y1$phases, FUN = plot, type = 'o',
#'  main = 'phases extracted as a seasoanal component of ts object')
#'
#' ## returning the phases as original data
#' y2 <- phases(dx1, ts_freq = 12, returned = 'original')
#' y2
#' lapply (X = y2$phases, FUN = plot, type = 'o',
#' main = 'phases extracted as raw')
#' par(default_par)
#'
#' ## 2 years of uni-modal pixel with many Nas towards the end.
#' dx2 <- dx1
#' dx2[21:length(dx2)] <- NA
#' default_par <- par()
#' layout(rbind(c(1, 1), c(2, 3)))
#' par(mar = c(2, 2, 1, 1))
#' plot(dx2, type = 'o')
#'## returning the phases as seasonal data.
#' y1 <- phases(dx2, ts_freq = 12)
#' y1
#' lapply (X = y1$phases, FUN = plot, type = 'o')
#'
#' ## returning the phases as original data
#' y2 <- phases(dx2, ts_freq = 12, returned = 'original')
#' y2
#' lapply (X = y2$phases, FUN = plot, type = 'o') # note 1 phase and the starting point
#' par(default_par)
#' @export
phases <- function(x, type = c("v_points", "peaks"), n_criticals = 1, steps = 2,
                   ts_freq = 23, returned = c("ts_seasonal", "original")) {
    type <- match.arg(type); returned <- match.arg(returned)
    # keeping a copy as x will be altered after interpolate_gaps
    x_copy <- x
    # replacing NAs (if any) and bringing back the names
    x <- chillR::interpolate_gaps(x)$interp ; names(x) <- names(x_copy)
    # make a time series and extract latter the one cyle.
    x_decomp <- decompose(ts(x, frequency = ts_freq))
    # this is to make sure the data start with the biggining of a season
    x <- as.numeric(x_decomp$seasonal)  # getting the seasonal component
    names(x) <- names(x_decomp$x)  # giving it the corresponding names
    # as in the raw data stored in the ts object
    cycle_ids <- which(x == min(x))  # range of indeces of complete cycles.
    cycle_ids <- cycle_ids[1:2]  # taking just the first 2, so just one cycle.
    x <- x[cycle_ids[1]:cycle_ids[2]]  # subesetting our interpolated data with the
    if (type == "v_points") {
        # if we want to have the data broken down as 1:v_points1, v_point_1:end in other words if we want separe the
        # seasons.
        pts <- get_falls(x = x, n_v_shape = n_criticals, steps = steps)  # this will return the index of the v_point.
    } else if (type == "peaks") {
        pts <- get_falls(x = -x, n_v_shape = n_criticals, steps = steps)
    }
    # lastly use the range of data indexes to break it down into the corresponding pieces: for each range of index,
    # subset the data.
    critics <- c(pts, length(x))  # take just from the v_point to the end,
    # the part 1:v_point will be handled in the lapply loop below if we want get back the seasonal data (not the
    # original)
    if (returned == "ts_seasonal") {
        j <- 1  # here is the index of the first value.
        critics <- lapply(critics, function(i) {
            dat <- x[j:i]  # subset the data within this range
            j <<- i  # j was equal to 1 (the index of the first value) at the beginning but change it values
            # but we want it change to v_point, then the index of the last value upon itteration.
            dat  # will get a list of the phases back.
        })
    } else if (returned == "original") {
        # if we want the original data instead.
        j <- 1
        critics <- lapply(critics, function(i) {
            dat <- x_copy[j:i]
            j <<- i
            dat
        })
    }
    # let's have some names for slots and a class holding them
    names(critics) <- c(paste("phase", rep(1:length(critics)), sep = "_"))
    critics <- list(phases = critics, criticals = pts)
    class(critics) <- "phases"
    return(critics)
}

#' @rdname phases
#' @importFrom grDevices boxplot.stats
#' @importFrom stats decompose runif ts
#' @export
seasons <- function(x) {
    # setting seed for reproductibility
    set.seed(123)
    test <- has_peak(x)$test
    if (test == TRUE) {
        # calculating the NDVI ratio
        ratio_numerator <- x - min(x, na.rm = T)
        ratio_denominator <- max(x, na.rm = T) - min(x, na.rm = T)
        x <- ratio_numerator/ratio_denominator
        # getting a sample of values between the 2 first minimal values for comparaison purpose latter. the ideas is to
        # check at each time lag whether the ndvi value has significantly changed from the previous value or not.  if it
        # has change then the loop has to breack. otherwise continue onto the next values.  4 (or say 4 +1 = 5) values
        # because it looks like a minimum of 5 values are needed to detect a difference
        minVi <- runif(4, min(x), min(x[x != min(x)]))
        # making 1 cycle x <- c(x, x[1]) place holder for the loop
        id <- vector(mode = "numeric", length = length(x))
        for (i in 1:length(x)) {
            # subset the data from the first minimal value to the ith value
            ref <- c(minVi, x[i])
            # check whether the ith value is an outlier
            ref <- ref[ref %in% boxplot.stats(ref)$out]
            # if it is not, continue onto the next one
            if (length(ref) == 0)
                next
            # if it is return the indices index and break the loop, you are done!
            if (length(ref) == 1) {
                id[i] <- i
                break
            }
        }
        for (i in length(x):1) {
            # subset the data from the first minimal value to the ith value
            ref <- c(minVi, x[i])
            # check whether the ith value is an outlier
            ref <- ref[ref %in% boxplot.stats(ref)$out]
            # if it is not, continue onto the next one
            if (length(ref) == 0)
                next
            # if it is return the indices index and break the loop, you are done!
            if (length(ref) == 1) {
                id[i] <- i  #(length(x) - (i-1))
                break
            }
        }
        id <- id[id != 0]
        seas <- list(season = x, season_dates = c(begin = as.Date(names(x)[id][1]), end = as.Date(names(x)[id][2])))
        class(seas) <- "seasons"
    } else if (test == FALSE) {
        seas <- list(season = x, season_dates = c(begin = NA, end = NA))
    }
    return(seas)
}
