\name{logLik.riley}
\alias{logLik.riley}

\title{Print the log-likelihood}

\description{This function provides the (restricted) log-likelihood of the fitted model.}

\usage{logLik.riley(object, ...)}

\arguments{
  \item{object}{a \code{riley} object.}
  \item{\dots}{arguments to be passed on to other functions, currently ignored}
}

\value{The log-likelihood of a fitted \code{riley} model}

\author{Thomas Debray <thomas.debray@gmail.com>}

\examples{
data(Scheidler)
ds <- Scheidler[which(Scheidler$modality==1),]
fit <- riley(ds,correction.control="single")
logLik(fit)
}

\keyword{ likelihood }
\keyword{ riley }
\keyword{ bivariate }
\keyword{ meta-analysis }
