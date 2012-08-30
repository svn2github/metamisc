\name{riley-package}
\alias{riley-package}
\alias{riley}
\docType{package}
\title{Diagnostic and prognostic meta-analysis}
\description{This package provides functions for diagnostic and prognostic meta-analyses. It estimates univariate, bivariate and multivariate models, and allows the aggregation of previously published prediction models with new data.}
\details{
\tabular{ll}{
Package: \tab metamisc\cr
Type: \tab Package\cr
Version: \tab 0.0.2\cr
Date: \tab 2012-08-30\cr
License: \tab GPL-2 \cr
}
The package provides tools for the meta-analysis of diagnostic accuracy data. For this the number true positives (TP), false negatives (FN), true negatives (TN) and false postives (FP) for each study must be known. The package can fit the alternative model for bivariate random-effects meta-analysis by Riley et al. (2008). The main function to fit such model is \code{\link{riley}} and standard methods are available for the output of this function.
}
\author{
Thomas Debray <thomas.debray@gmail.com>

Maintainer: Thomas Debray <thomas.debray@gmail.com>
}
\references{
Reitsma, J., Glas, A., Rutjes, A., Scholten, R., Bossuyt, P., & Zwinderman, A. (2005).
\dQuote{Bivariate analysis of sensitivity and specificity produces informative summary
measures in diagnostic reviews.} \emph{Journal of Clinical Epidemiology}, \bold{58}, 982--990.

Riley, RD., Thompson, JR., & Abrams, KR. (2008).
\dQuote{An alternative model for bivariate random-effects meta-analysis when the within-study correlations are unknown.} \emph{Biostatistics}, \bold{9}, 172--186.
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation ~~
~~ directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
~~ simple examples of the most important functions ~~
}
