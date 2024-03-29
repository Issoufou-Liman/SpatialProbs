% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_boxplot_range_1d.R
\name{get_boxplot_range_1d}
\alias{get_boxplot_range_1d}
\title{Create data ranges based boxplots statistics and outliers}
\usage{
get_boxplot_range_1d(x, split_IQR = TRUE, na.rm = TRUE, type = 5,
  ...)
}
\arguments{
\item{x}{numeric. A vector from which the ranges are to be derived based on quantiles.
NA and NaN values are not allowed unless na.rm is TRUE.}

\item{split_IQR}{logical. Should the inter-quartile range be split at the median
to form different ranges? The default is TRUE}

\item{na.rm}{logical; if true, any \code{\link{NA}} and \code{NaN}'s
    are removed from \code{x} before the quantiles are computed.}

\item{type}{an integer between 1 and 9 selecting one of the
    nine quantile algorithms detailed below to be used.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
A numeric vector containing the bounds labelled as out_min (the lower outlier),
ymin (the extreme of the lower whisker), lower (the lower ‘hinge’), middle (the median),
upper (the upper ‘hinge’), ymax (the extreme of the upper whisker), out_max (the upper outlier).
}
\description{
Take a vector of data and returns the bounds corresponding the boxplot specifications.
}
\details{
The ranges are bounded by the lower outlier (if any), extreme of the lower whisker,
the lower ‘hinge’, the median, the upper ‘hinge’, the extreme of the upper whisker,
and the upper outlier (if any).
}
\examples{
sample <- rnorm(100)
get_boxplot_range_1d (sample)
}
\author{
Issoufou Liman
}
