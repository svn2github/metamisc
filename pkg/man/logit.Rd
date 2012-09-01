\name{logit}
\alias{logit}

\title{
Apply logit tranformation
}
\description{
 Transforms values between 0 and 1 to values between -Inf and Inf.
}
\usage{
logit(x)
}

\arguments{
  \item{x}{
   A vector of numerics (between 0 and 1)
}
}

\value{
 A vector of numerics (between -Inf and Inf).
}

\author{
  Thomas Debray <thomas.debray@gmail.com>
}

\seealso{
  \code{\link{inv.logit}}
}

\examples{
  inv.logit(0.50) #A probability of 50% corresponds to a linear predictor of 0.
}

\keyword{logit}
\keyword{transformation}
