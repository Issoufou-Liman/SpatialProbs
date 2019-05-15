#' Check seasonal change point in a vector of data
#'
#' Takes a numeric vector, locates and returns the indices of a maximum number of n_v_shape points
#' given a specified number of relative points (that are higher at both side of a v shape point).
#' @author Issoufou Liman
#' @param x A numeric vector. in which to check the indices of the fall (v-like) shape points.
#' @param n_v_shape Integer. The number of points at which x takes a v-like shape. Default to 1.
#' @param steps Integer. The minimum number of points at either sides of each v-like point
#' to be considered higher for it to qualify as v-like point. Default to 3
#' @param force_steps Logical. Should steps argument be enforced on both side each v_shape point?
#' This is related to v-like points at the beiggining and/or end. Default to TRUE
#' @return A vector of `lenght (n_v_shape)`, at most, containing the indices of v-like points
#' @details The number of points at both side of a v-like point will determine the number of
#' v-like points to be return when there are many candidates since the higher
#' the steps argument, the less the number of final points to be returned. while
#' `get_falls ()` will increase step-wise the number of steps to get closer to the number
#' of specified v-like points (`n_v_shape`), it cannot handle cases where there are 2 or
#' more equal points located at a v shape. Called on `- x`, `get_falls ()` will return
#' the peaks instead of v-like points.
#' @examples
#' dx11 <- c(1.30, 1.15,  1.50,  2.00,  2.01,  3.00, 3.20,  4.76,  3.50,  3.00,  2.40,  2.00,  1.50)
#' dx12 <-c(1.29, 1.1, 1.49, 1.99, 2, 3.1, 4.5, 4, 2.8, 2.5, 2.3, 1.6, 1.59)
#' dx1 <- c(dx11, dx12)
#' ## checking existing v-like shape
#' v_pts <- get_falls (dx1, n_v_shape = 1, steps = 3)
#' ## ploting the v_pts in red
#' plot(dx1, col=ifelse(dx1==dx1[v_pts], 'red', 'black'),
#' pch=ifelse(dx1==dx1[v_pts], 19, 1), type = 'o')
#' ## checking non existing v-like shape
#' get_falls (dx11, n_v_shape = 1, steps = 3)
#' get_falls (dx12, n_v_shape = 1, steps = 3)
#' @export
get_falls <- function(x, n_v_shape = 1, steps = 3, force_steps = TRUE) {
    # points where a curve has v shape are point where the derivate of the signs of the derivates of the points
    # representing the curve is equal to 2: i.e. diff(sign(diff(x)))==2 suppose the vector starts a -Infinity to
    # handle cases of v_shape a the beggining
    shape <- diff(sign(diff(x)))
    # applying the concept above to every points and returning the the index of those points that qualify for v_shape
    # areas Note that the number of v_shapes depends on the steps (number of points aside the targeted point)
    # considered.  Therefore we start with a small number of steps and increase it if needed depending on the number
    # of desired peacks.  This is achieved using the deep assignment arrow <<- as below inside the increase_steps ()
    # function.  This helper function increase the steps (number of points to be higher than the targeted point at
    # both side) each time repeat () is called.  This way, the (main) get_falls () function will returned its
    # value depending on the number of v_shape (n_v_shape argument) specified.  <<-, instead of creating a variable
    # in the current environment, as as does <-, modifies a variable existing in the upper level (parent
    # environments).  (see ?'<<-')
    steps <- steps - 1  # substract 1 to account for value specified as argument
    increase_steps <- function() {
        steps <<- steps + 1  # remembering the previous value of steps.
        return(steps)
    }
    repeat {
        steps <- increase_steps()
        v_shape <- sapply(X = which(shape == 2), FUN = function(i) {
            # apply this function to every points where diff(sign(diff()))==2 steps is the minimum number of points to
            # consider at both side of the targeted point.  the +1 account for the differences in index due the derivative as
            # diff function uses x[(1+lag):n] - x[1:(n-lag)].  see ?diff
            before <- i - steps + 1
            # making sure to not get out of bound index
            before <- ifelse(before > 0, before, 1)
            after <- i + steps + 1
            after <- ifelse(after < length(x), after, length(x))
            # getting steps points before the targeted point and steps points after the targeted point
            test <- x[c(before:i, (i + 2):after)]
            # getting the targeted point
            target <- x[i + 1]
            # we know the point is some how located at a v shape area but to which extent? so, testing the extent to which it
            # qualifies (given steps and n_v_shape arguments).  force to have the required number of steps on each side if
            # all points are higher than the targeted point,
            if (force_steps == TRUE) {
                if (all(test > target) && length(test) >= 2 * steps) {
                  return(i + 1)  # then return the target.
                } else {
                  # otherwise, return null
                  return(numeric(0))
                }
            } else if (force_steps == FALSE) {
                if (all(test > target)) {
                  # if all points are higher than the targeted point,
                  return(i + 1)  # then return the target.
                } else {
                  # otherwise, return null
                  return(numeric(0))
                }
            }
        })
        # save to v_shape object the vector of the index of the v shape points.
        v_shape <- c(unlist(v_shape))
        if (length(v_shape) <= n_v_shape) {
            # stop at most when you reach n_v_shape as i asked!
            break
        }
    }
    return(v_shape)
}
