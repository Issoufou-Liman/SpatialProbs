#' A constructor function for class 'phase'
#'
#' Takes vector, makes phases, and returns an s3 object of class phase. A phase may not
#' a lot of sense if it does have a peack.
#' The raison behind the class phases is better illustrated using an example. considering a bi-modal season having
#' consecutives NAs at the end (see example), then the last season (see class season) may not make
#' sense as it will be truncated or incomplete.Thus, this class will allow further analysis to be
#' conducted on real seasonal data.
#' @author Issoufou Liman
#' @param vec a numeric vector which names are dates at which the data were acquired. These are typically pixels values
#' extracted from a stack of NDVI time series.
#' @param type character string. either 'v_points'or 'peaks' depending on how the vec should be broken into phases.
#' 'v_points' (default) may make more sens for expression seasonal data.
#' @param n_criticals Numeric. The number of points to be considered. Default to 1.
#' @param steps see check_v_shapes ()
#' @param ts_freq The frequence of time series (see ?ts, argument frequency).
#' default is 1 year (23 = 2 values per months 16 days temporal resolution data).
#' @param returned character string either 'ts_seasonal' or 'original'. if 'ts_seasonal' (default),
#' the returned phases data will be seasonal component of the decomposed time series (see ?decompose).
#' if 'original', the original data will be returned.
#' @return an object of class phases comprise a list of:
#'   - phases, list of phases
#'   - criticals, index of the points separating the phases
#' @details see get_phases method for this class.
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
#' lapply (X = y1$phases, FUN = plot, type = 'o', main = 'phases extracted as a seasoanal component of ts object') # note 1 phase and the starting point
#'
#' ## returning the phases as original data
#' y2 <- phases(dx1, ts_freq = 12, returned = 'original')
#' y2
#' lapply (X = y2$phases, FUN = plot, type = 'o', main = 'phases extracted as raw') # note 1 phase and the starting point
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
phases <- function(vec, type = c("v_points", "peaks"), n_criticals = 1, steps = 2, ts_freq = 23, returned = c("ts_seasonal",
    "original")) {
    type <- match.arg(type)
    returned <- match.arg(returned)
    saved_names <- names(vec)  # making a copy of the names which are the dates we are targeting here as after
    # the interpolation the names will be lost.
    vec_copy <- vec  # also making a copy of the data
    vec <- chillR::interpolate_gaps(vec)  # replacing the NAs
    vec <- vec$interp  # taking the interpolated values
    names(vec) <- saved_names  # bringing back the names
    vec_decomp <- decompose(ts(vec, frequency = ts_freq))  # make a time series and extract latter the one cyle.
    # this is to make sure the data start with the biggining of a season
    vec <- as.numeric(vec_decomp$seasonal)  # getting the seasonal component
    names(vec) <- names(vec_decomp$x)  # giving it the corresponding names
    # as in the raw data stored in the ts object
    cycle_ids <- which(vec == min(vec))  # range of indeces of complete cycles.
    cycle_ids <- cycle_ids[1:2]  # taking just the first 2, so just one cycle.
    vec <- vec[cycle_ids[1]:cycle_ids[2]]  # subesetting our interpolated data with the
    if (type == "v_points") {
        # if we want to have the data broken down as 1:v_points1, v_point_1:end in other words if we want separe the
        # seasons.
        pts <- check_v_shapes(vektor = vec, n_v_shape = n_criticals, steps = steps)  # this will return the index of the v_point.
    } else if (type == "peaks") {
        pts <- check_v_shapes(vektor = -vec, n_v_shape = n_criticals, steps = steps)
    }
    # lastly use the range of data indexes to break it down into the corresponding pieces: for each range of index,
    # subset the data.
    critics <- c(pts, length(vec))  # take just from the v_point to the end,
    # the part 1:v_point will be handled in the lapply loop below if we want get back the seasonal data (not the
    # original)
    if (returned == "ts_seasonal") {
        j <- 1  # here is the index of the first value.
        critics <- lapply(critics, function(i) {
            dat <- vec[j:i]  # subset the data within this range
            j <<- i  # j was equal to 1 (the index of the first value) at the beginning but change it values
            # but we want it change to v_point, then the index of the last value upon itteration.
            dat  # will get a list of the phases back.
        })
    } else if (returned == "original") {
        # if we want the original data instead.
        j <- 1
        critics <- lapply(critics, function(i) {
            dat <- vec_copy[j:i]
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
