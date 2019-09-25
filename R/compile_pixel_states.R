#' Deriving qualitative pixel states from raster data
#'
#' make_pixel_states takes a raster stack or brick to qualify pixel
#' following scale either defined by the user or derived from the natural
#' breaks of data following the following boxplot statistics and outliers.
#' @author Issoufou Liman
#' @param x a raster stack or brick
#' @inheritParams make_pixel_states
#' @inheritParams raster::focal
#' @param as_states Logical, are the states already made (e.g. x is an output of make_pixel_states function call)?
#' If TRUE, the x is assumed to be a stack of ready made states and is directly processed.
#' @details
#' The ranges are bounded by the lower outlier (if any), extreme of the lower whisker,
#' the lower ‘hinge’, the median, the upper ‘hinge’, the extreme of the upper whisker,
#' and the upper outlier (if any).
#' @return
#' A raster brick where each layers correspond to one of the states.
#' @examples
#' sample <- rnorm(100)
#' get_boxplot_range_1d (sample)
#' @importFrom raster focal ncell extract xyFromCell predict overlay
#' @importFrom caret confusionMatrix
#' @export
compile_pixel_states <- function(x, w = matrix(1, 3, 3),
                                 split_IQR = FALSE, custum_rclmat=NULL,
                                 sample = FALSE, size = 1000, ties = 'NA', op = c("sampler", "proba"), as_states=FALSE, inparallel = NULL, ...){
  if (!as_states){
    x <- make_pixel_states (x, split_IQR = split_IQR, custum_rclmat=custum_rclmat,
                            sample = sample, size = size, ties = 'NA', op='proba', inparallel = inparallel)
  }

  # fun: a function to be used by calc to decide the node states depending on unique
  # highest probability
  f <- function(x){
    out <- which(x==max(x, na.rm = TRUE)) # getting the indices of the max values
    out <- ifelse(length(out) == 1, out, NA) # if only one recode it to the layer number else assign NA to it
    out
  }
  output <- calc(x, f)
  ## Majority rule for 1st order confusing pixels ####

  # Deciding the confusing pixels (pixels that have more than 1 max values)
  # using majority rule via a 3 by 3 moving window.

  # fun: a function to be used by focal
  g <- function(x){
    out <- modal(x, na.rm = TRUE, ties = 'NA') # getting the modal values
    out <- ifelse(length(out) == 1, out, NA) # if only one recode it to the most likely in the window else assign NA to it
    out
  }
  # # the moving window
  # window <- matrix(1, 3, 3)
  # calling focal with well tuned specifications.
  # - 1) the fun and the windows above
  # - 2) na.rm = FALSE, consider NA during focal computation which will then be ignore when computing the mode in fun.
  # - 3) NAonly = TRUE, operate only on NA which we are interested in only
  # - 4) pad = TRUE to get ride of egde effect.
  output <- focal(x=output, w=w, fun=g , na.rm=FALSE, NAonly=TRUE, pad=TRUE)
  ## Prediction for 2nd order confusing pixels ####
  # extract_many_max_cells, function for extraction the values of the confusing pixels
  extract_many_max_cells <- function(r_stack, include_coords=FALSE, include_cells=FALSE){
    cells <- sapply(seq_len(ncell(r_stack)), function (i){
      tmp <- which(r_stack[i] == max(r_stack[i], na.rm = TRUE))
      ifelse(length(tmp) > 1, i, NA)
    }, simplify = T)
    cells <- cells[!is.na(cells)]
    out <- extract(r_stack, cells)
    if (include_coords){
      coordinate <- xyFromCell(r_stack, cell=cells, spatial=FALSE)
      out <- data.frame(coordinate, out)
    }
    if(include_cells){
      cell_numbers <- cells
      out <- data.frame(cell_numbers, out)
    }
    out
  }
  # extracting the cell values only if there are 2 or more max in the cell
  tmp <- extract_many_max_cells(x, include_cells = FALSE, include_coords = FALSE)
  # matrix to dataframe
  tmp <- as.data.frame(tmp)
  # extracting the max values for each row
  tmp_1 <- lapply(1:nrow(tmp), function(i) {
    j <- tmp[i,]
    j[which(j==max(j))]
  })
  # long data formating
  tmp_1 <- reshape2::melt(tmp_1)
  # fitting a decision tree model
  cart.model = party::ctree(variable ~ value, tmp_1)
  # filling in the remaining pixels with the fitted model prediction
  predict_fun <- function(x, .cart.model, original_order){
    out <- which(x==max(x, na.rm = TRUE))
    if(length(out) > 1){
      out <- as.factor(predict(.cart.model, type = "response", simplify=T, newdata=data.frame(value=max(x))))
      predicted_order <- .cart.model@responses@levels$variable[out]
      out <- which(original_order==predicted_order)
    } else {
      out <- NA
    }
    out
  }
  tmp <- calc(x=x, fun = function(i) predict_fun(i, .cart.model = cart.model, original_order=names(tmp)))
  h <- function(x,y){
    x[is.na(x)] = y[is.na(x)]
    x
  }
  output <- overlay(output, tmp, fun=h)
  list (output=output, model=cart.model, accuracy=confusionMatrix(predict(cart.model), tmp_1$variable), data=tmp_1)

}
