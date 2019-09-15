#' Spatially explicit Principal Component Analysis (PCA)
#'
#' Performs PCA on multi layers raster data. This function is based on rasterPCA from the RStoolbox package (Leutner et al., 2018).
#' @param x RasterBrick RasterStack
#' @param nSamples Integer or NULL. Number of pixels to sample for PCA fitting. If NULL, all pixels will be used.
#' @param nComp number of principal components to consider
#' @param n Integer, sample size for specific sampling method sampling_type (see sp::spsample)
#' @param spca Logical. If TRUE, perform standardized PCA. Corresponds to centered and scaled input image.
#' This is usually beneficial for equal weighting of all layers. (FALSE by default)
#' @param maskCheck Logical. Masks all pixels which have at least one NA.
#' @param na.action a function which indicates what should happen when the data contain NAs. The default is set by the na.action
#' setting of options, and is na.fail if that is unset. The ‘factory-fresh’ default is na.omit.
#' @param scale. A logical value indicating whether the variables should be scaled to have unit variance
#' before the analysis takes place. The default is FALSE for consistency with S, but in general scaling
#' is advisable. Alternatively, a vector of length equal the number of columns of x can be supplied. The value is passed to scale.
#' @param sampling_type Type specific sampling to be used
#' @param ... Additional arguments passed to writeRaster.
#' @references
#' Benjamin Leutner, Ned Horning and Jakob Schwalb-Willmann (2018). RStoolbox: Tools for Remote Sensing Data Analysis. R
#' package version 0.2.3. https://CRAN.R-project.org/package=RStoolbox
#' @return RasterBrick
#' @importFrom  methods as
#' @importFrom stats na.exclude na.omit prcomp var
#' @importFrom raster crs crs<- extent cellStats mask crop
#' @importFrom sp spsample
#' @export
raster_pca <- function (x, nSamples = NULL, n = ceiling(ncell(x)*0.2), nComp = nlayers(x), spca = FALSE,
                        maskCheck = TRUE,na.action=na.exclude, scale. = TRUE, sampling_type="regular", ...)
{
  set.seed(25)
  x <- scale(x, center = T, scale = T)

  if (nlayers(x) <= 1)
    stop("Need at least two layers to calculate PCA.")
  ellip <- list(...)
  if ("norm" %in% names(ellip)) {
    warning("Argument 'norm' has been deprecated. Use argument 'spca' instead.\nFormer 'norm=TRUE' corresponds to 'spca=TRUE'.",
            call. = FALSE)
    ellip[["norm"]] <- NULL
  }
  if (nComp > nlayers(x))
    nComp <- nlayers(x)
  if (!is.null(nSamples)) {
    trainData <- sampleRandom(x, size = nSamples, na.rm = TRUE)
    trainData <- t(trainData)
    trainData <- trainData[ , apply(X=trainData, MARGIN=2, FUN=var, na.rm = TRUE) != 0]
    #trainData <- trim
    colnames(trainData) <- paste0("x", col(trainData)[1, ])
    trainData <- scale(trainData, center = TRUE, scale = scale.)
    trainData <- data.frame(t(na.omit(t(trainData))))


    if (nrow(trainData) < nlayers(x))
      stop("nSamples too small or x contains a layer with NAs only")
    model <- prcomp(~., data = trainData, na.action=na.action, scale. = scale.)

  } else {
    if (maskCheck) {
      totalMask <- !sum(calc(x, is.na))
      if (cellStats(totalMask, sum) == 0)
        stop("x contains either a layer with NAs only or no single pixel with valid values across all layers")
      x <- mask(x, totalMask, maskvalue = 0)
    }
    sp_trainer <- as(extent(x), "SpatialPolygons"); crs(sp_trainer) <- crs(x); sp_trainer <- crop(sp_trainer, x)
    sp_trainer <- spsample(sp_trainer, n=n, type=sampling_type)
    trainData_sp <- extract(x, sp_trainer, method='bilinear', fun=mean, na.rm=T, sp=TRUE)
    trainData <- trainData_sp@data
    #colnames(trainData) <- paste0("x", col(trainData)[1, ])
    trainData <- as.data.frame(trainData)
    trainData <- scale(trainData, center = TRUE, scale = scale.)
    #trainData <- data.frame(t(na.omit(t(trainData))))
    trainData <- data.frame(na.omit(trainData))

    model <- prcomp(trainData, retx = TRUE)
    #model <- prcomp(~., data=trainData, na.action=na.action, scale. = scale.)
    # formule <- as.formula(paste("~", paste(names(trainData), collapse = " + ")))
    # model <- prcomp(formule, data=trainData, na.action=na.action, scale. = scale.)
    #

  }
  .paraRasterFun <- function(raster, rasterFun, args = list(), wrArgs = list()){
    if (isTRUE( getOption('rasterCluster'))) {
      do.call("clusterR", args = c(list(x = raster, fun = rasterFun, args=args), wrArgs))
    } else {
      do.call("rasterFun", args=c(raster, args, wrArgs))
    }
  }
  # x <- t(x)
  # names(x) <- colnames(trainData)

  out <- .paraRasterFun(x, rasterFun = raster::predict,
                        args = list(model = model, na.rm = TRUE, index = 1:nComp),
                        wrArgs = ellip)
  rownames(model$rotation) <- as.character(as.Date(substring(rownames(model$rotation), 2), "%Y.%m.%d"))
  structure(list(call = match.call(), model = model, map = out),
            class = "raster_pca")
  # model
}
