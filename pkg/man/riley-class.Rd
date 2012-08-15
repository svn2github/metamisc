\name{riley-class}
\alias{riley-class}
\alias{plot.riley}
\alias{summary.riley}
\alias{print.riley}

\concept{random effects model}
\concept{diagnostic meta-analysis}
\concept{meta-analysis}

\title{
Methods for \code{riley} objects.
}
\description{
Objects of the class \code{\link{riley}} are output by the function with the same name. 
}
\usage{
\method{print}{riley}(x, ...)
\method{summary}{riley}(object, level = 0.95, ...)
\method{plot}{riley}(x, plotsumm = TRUE, plotnumerics = TRUE, level = 0.95, 
     main="", ylim = c(0,1), xlim = c(0,1), pch = 1, lty = 1, lwd = 1, pch = 1, 
     cex.numerics = 0.45, add = F, ...)
}
\arguments{
  \item{object}{a \code{riley} object.}
  \item{level}{numeric, the level for calculations of confidence intervals (\code{summary}) or regions (\code{plot})}
  \item{add}{logical, should the confidence region be added to the current plot?}
  \item{plotsumm}{logical, should the summary pair of sensitivity and false positive rate together with its confidence region be plotted?}
  \item{plotnumerics}{logical, should the plot contain a summary table of sensitivity and false positive rate?}
  \item{ylim}{numeric of length 2, which section of the sensitivities to plot?}
  \item{xlim}{numeric of length 2, which section of the false positive rates to plot?}
  \item{pch}{integer, symbol for the pair of mean sensitivity and false positive rate}
  \item{lty}{integer, line type of confidence curve}
  \item{lwd}{integer, line width of the confidence curve}  
  \item{cex.numerics}{numeric, text size 
  \item{...}{arguments to be passed on to other functions}
}
\details{

}
\value{

}
\author{
Thomas Debray <thomas.debray@gmail.com>
}

\seealso{
\code{\link{riley}}
}

\keyword{regression}
\keyword{classes}
