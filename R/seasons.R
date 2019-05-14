#' @author Issoufou Liman
seasons <- function(vektor) {
    # setting seed for reproductibility
    set.seed(123)
    test <- has_peak(vektor)$test
    if (test == TRUE) {
        # calculating the NDVI ratio
        ratio_numerator <- vektor - min(vektor, na.rm = T)
        ratio_denominator <- max(vektor, na.rm = T) - min(vektor, na.rm = T)
        vektor <- ratio_numerator/ratio_denominator
        # getting a sample of values between the 2 first minimal values for comparaison purpose latter. the ideas is to
        # check at each time lag whether the ndvi value has significantly changed from the previous value or not.  if it
        # has change then the loop has to breack. otherwise continue onto the next values.  4 (or say 4 +1 = 5) values
        # because it looks like a minimum of 5 values are needed to detect a difference
        minVi <- runif(4, min(vektor), min(vektor[vektor != min(vektor)]))
        # making 1 cycle vektor <- c(vektor, vektor[1]) place holder for the loop
        id <- vector(mode = "numeric", length = length(vektor))
        for (i in 1:length(vektor)) {
            # subset the data from the first minimal value to the ith value
            ref <- c(minVi, vektor[i])
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
        for (i in length(vektor):1) {
            # subset the data from the first minimal value to the ith value
            ref <- c(minVi, vektor[i])
            # check whether the ith value is an outlier
            ref <- ref[ref %in% boxplot.stats(ref)$out]
            # if it is not, continue onto the next one
            if (length(ref) == 0)
                next
            # if it is return the indices index and break the loop, you are done!
            if (length(ref) == 1) {
                id[i] <- i  #(length(vektor) - (i-1))
                break
            }
        }
        id <- id[id != 0]
        seas <- list(season = vektor, season_dates = c(begin = as.Date(names(vektor)[id][1]), end = as.Date(names(vektor)[id][2])))
        class(seas) <- "seasons"
    } else if (test == FALSE) {
        seas <- list(season = vektor, season_dates = c(begin = NA, end = NA))
    }
    return(seas)
}
