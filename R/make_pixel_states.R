#' Deriving qualitative pixel states from raster data
#'
#' make_pixel_states takes a raster stack or brick to qualify pixel
#' following scale either defined by the user or derived from the natural
#' breaks of data following the following boxplot statistics and outliers.
#' @author Issoufou Liman
#' @param x a raster object
#' @inheritParams reclassifion_matrix
#' @inheritParams raster::calc
#' @inheritParams raster::modal
#' @param ... Additional arguments as for writeRaster
#' @param inparallel numeric indicating the number of processes to run in parallel
#' @details
#' The ranges are bounded by the lower outlier (if any), extreme of the lower whisker,
#' the lower ‘hinge’, the median, the upper ‘hinge’, the extreme of the upper whisker,
#' and the upper outlier (if any).
#' @return
#' A raster brick where each layers correspond to one of the states.
#' @examples
#' sample <- rnorm(100)
#' get_boxplot_range_1d (sample)
#' @importFrom raster calc reclassify canProcessInMemory writeRaster addLayer brick rasterTmpFile modal extension
#' @export
make_pixel_states <- function(x, split_IQR = FALSE, custum_rclmat=NULL,
                              sample = FALSE, size = 1000, ties = 'NA', filename = '', inparallel = NULL,...){
  tmp_filename_1 <- rasterTmpFile(prefix = 'r_tmp_1_')
  tmp_filename_2 <- rasterTmpFile(prefix = 'r_tmp_2_')
  tmp_filename_3 <- rasterTmpFile(prefix = 'r_tmp_3_')
  rclmat <- reclassifion_matrix (x = x, split_IQR = split_IQR, custum_rclmat = custum_rclmat, sample = sample, size = size)
  disaggregate_reclmat <- function(rclmat){
    tmp <- 1:nrow(rclmat); names(tmp) <- rownames(rclmat)
    sapply(tmp, function(i){
      rclmat <- rclmat[i, , drop = FALSE]
      if(i == 1){
        rclmat <- rbind(rclmat,
                        c(rclmat[2], Inf, 0))
      } else if (i == length(tmp)){
        rclmat <- rbind(c(-Inf, rclmat[1], 0),
                        rclmat)
      } else {
        rclmat <- rbind(c(-Inf, rclmat[1], 0),
                        rclmat,
                        c(rclmat[2], Inf, 0))
      }
      rclmat[,3] <- as.integer(as.logical(rclmat[, 3]))
      rownames(rclmat) <- NULL
      rclmat
    }, simplify = FALSE, USE.NAMES = TRUE)
  }
  if(class(x) == 'RasterLayer'){ # Case of single layer
    rclmat <- disaggregate_reclmat(rclmat = rclmat)
    if(canProcessInMemory(x, 4*length(rclmat))){ # Single layer can be processed in memory
      out <- sapply(1:length(rclmat), function(i){
        reclassify(x, rclmat[[i]])
      })
      out_target <- brick(out)
      if (filename != '') {
        writeRaster(out_target, filename = filename, ...) # problems with native raster format
      }
      return(out_target)
    } else { # single layer can't be processed in memory
      for(i in 1:length(rclmat)){
        tmp_filename_2 <- rasterTmpFile(prefix = paste0(names(rclmat)[[i]], "_"))
        if(i == 1){
          out_target <- reclassify(x, rclmat[[i]], filename = tmp_filename_1)
        } else {
          out_tmp <- reclassify(x, rclmat[[i]], filename = tmp_filename_2)
          out_target <- addLayer(out_target, out_tmp)
        }
        out_target
      }
      names(out_target) <- names(rclmat)
      if (filename != '') {
        out_target <- writeRaster(out_target, filename = filename, ...) # problems with native raster format
      }
      return(out_target)
    }
  } else { # Case of Multi layer
    rclmat <- sapply(X = rclmat, FUN = disaggregate_reclmat, simplify = FALSE, USE.NAMES = TRUE)
    elem_size <- max(sapply(rclmat, length))
    complet_names <- rclmat[which(lapply(rclmat, length)==elem_size)]
    complet_names <- names(complet_names[[1]])
    rclmat <- sapply(complet_names, function(i){
      sapply(rclmat, '[[', i, simplify = FALSE)
    }, simplify = FALSE)
    rclmat <- sapply (rclmat, function (i) Filter(Negate(is.null), i), simplify = FALSE)
    if(canProcessInMemory(x, 4*elem_size*nlayers(x))){ # Multi layer can be processed in memory
      out <- rclmat
      out <- sapply(1:length(rclmat), function(i){ # something like 5 for boxplot stats and out
        sapply(1:length(rclmat[[i]]), function(j){ # something like 138 for a 3 year modis 16 days ts
          reclassify(x[[j]], rclmat[[i]][[j]])
        }, simplify = FALSE)
      }, simplify = FALSE)
      out <- sapply(out, function(i){
        calc(x=brick(i), fun = function (x) modal(x, ties = ties))
      }, simplify = FALSE, USE.NAMES = TRUE)
      out <- brick(out)
      if (filename != '') {
        out <- writeRaster(out, filename = filename,  ...)
      }
      return(out)
    } else { # Multi layer can't be processed in memory
      for(i in 1:length(rclmat)){ # something like 5 for boxplot stats and out
        if (filename == '') {
          tmp_filename_1 <- rasterTmpFile(prefix = paste0(names(rclmat)[[i]], "_"))
        } else {
          tmp_filename_1 <- gsub(pattern = extension(filename), replacement = paste0('_', i, extension(filename)), filename)
        }
        for (j in 1:length(rclmat[[i]])) { # something like 138 for a 3 year modis 16 days ts
          if(j==1){
            out <- reclassify(x[[j]], rclmat[[i]][[j]], filename = tmp_filename_1, overwrite = TRUE)
          } else {
            tmp_filename_2 <- rasterTmpFile(prefix = 'tmp_file_')
            out_tmp <- reclassify(x[[j]], rclmat[[i]][[j]], filename = tmp_filename_2, overwrite = TRUE)
            out <- addLayer(out, out_tmp)
          }
        }
        names(out) <- names(rclmat[[i]])
        if (i==1){
          out_target <- calc(x=out, fun = function (x) modal(x, ties = ties), filename = tmp_filename_3)
        } else {
          out_tmp <- calc(x=out, fun = function (x) modal(x, ties = ties), filename = tmp_filename_1, overwrite = TRUE)
          out_target <- addLayer(out_target, out_tmp)
        }
      }
      names(out_target) <- names(rclmat)
      if (filename != '') {
        out_target <- writeRaster(out_target, filename = filename,  ...)
      }
      return(out_target)
    }
  }
}
