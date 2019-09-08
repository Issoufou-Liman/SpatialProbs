#' Discretize raster
#'
#' Discretize raster based boxplots statistics and outliers
#' @author Issoufou Liman
#' @inheritParams raster::raster
#' @param split_IQR logical. Should the inter-quartile range be split at the median
#' to form different ranges? The default is TRUE
#' @details
#' The ranges are bounded by the lower outlier (if any), extreme of the lower whisker,
#' the lower ‘hinge’, the median, the upper ‘hinge’, the extreme of the upper whisker,
#' and the upper outlier (if any).
#' @return
#' A raster layer or list of raster layers
#' @importFrom raster nlayers reclassify
#' @export
discretize_raster <- function(x, split_IQR = FALSE){
  rclmat <- reclassifion_matrix (x, split_IQR = split_IQR)
  if (class(x) == 'RasterLayer'){
    reclassify(x, rclmat)
  } else if (class(x) %in% c('RasterBrick', 'RasterStack')){
    tmp <- 1:nlayers(x); names(tmp) <- names(x)
    sapply(X=tmp, function(i){
      reclassify(x[[i]], rcl=rclmat[[i]])
    }, simplify = FALSE, USE.NAMES = TRUE)
  }
}
