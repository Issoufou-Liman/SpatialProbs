#' Create reclassification matrix
#'
#' Create reclassification matrix based on boxplot statistics and outliers or user defined specifications.
#' @author Issoufou Liman
#' @param x numeric vector, data.frame, or raster object.
#' @inheritParams get_boxplot_range_1d
#' @param ... Additional arguments passed to get_boxplot_range_1d.
#' @param custum_rclmat numeric vector indicating the limits based on which to reclassify x.
#' @param size Integer. Sample size (if sample is TRUE or if the data cannot
#' be process in memory).
#' @param sample Logical. Should the data be sampled? If TRUE, the result is
#' based on a ramdon sampling of the original data.
#' @return
#' a matrix (if x is one dimensional) or list of matrix.
#' @examples
#' sample <- rnorm(100)
#' reclassifion_matrix (sample)
#' @importFrom raster canProcessInMemory getValues sampleRandom
#' @export
reclassifion_matrix <- function(x, split_IQR = FALSE, na.rm = TRUE , type = 5, ..., custum_rclmat=NULL,
                                sample = FALSE, size = 1000){
  if (!(class(x) %in% c("integer","numeric", 'matrix', 'data.frame', "RasterLayer", 'RasterBrick', 'RasterStack'))){
    stop(paste(deparse(substitute(xx)), "must be integer, numeric, matrix, data.frame, RasterLayer, RasterBrick, RasterStack"))
  } else {
    if (class(x) %in% c("RasterLayer", 'RasterBrick', 'RasterStack')){
      if((!canProcessInMemory(x)) | sample){
        x <- sampleRandom(x, size = size)
      } else {
        x <- getValues(x)
      }
    }
  }

  fun <- function(x, split_IQR, na.rm , type, ...,custum_rclmat){
    if(!is.null(custum_rclmat)){
      rclmat <- unique(custum_rclmat)
    } else {
      rclmat <- get_boxplot_range_1d(x, split_IQR = split_IQR, na.rm = na.rm , type = type, ...)
      rclmat <- unique(rclmat)
    }

    if (length(rclmat) <= 2){
      stop(paste('only', length(rclmat), "elemnts in the reclassification matrix."))
    }
    rclmat <- sort(c(rclmat, rclmat[2:(length(rclmat)-1)]), decreasing = FALSE)
    rclmat[1] <- (-Inf); rclmat[length(rclmat)] <- Inf
    rclmat <- matrix(rclmat, ncol = 2, byrow = TRUE)
    rclmat <- cbind(rclmat, 1:nrow(rclmat))
    return(rclmat)
  }
  if (class(x) %in% c('matrix', 'data.frame')){
    sapply(X=colnames(x), function(i){
      fun (x=x[, i], split_IQR = split_IQR, na.rm = na.rm , type=type, custum_rclmat = custum_rclmat)
    }, simplify = FALSE, USE.NAMES = TRUE)
  } else {
    fun (x=x, split_IQR = split_IQR, na.rm = na.rm , type=type, custum_rclmat = custum_rclmat)
  }
}
