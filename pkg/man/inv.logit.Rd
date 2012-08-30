\name{inv.logit}
\alias{inv.logit}

\title{
Inverse of logistic tranformation
}
\description{
Transforms a linear predictor into a probability.
}
\usage{
inv.logit(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A vector of numerics between -Inf and Inf.
}
}

\value{
A vector of numerics between 0 and 1.
}
\author{
Thomas Debray <thomas.debray@gmail.com>
}

\examples{
inv.logit(0) #A linear predictor of 0 corresponds to a probability of 50%
}

\keyword{logistic}
\keyword{inverse}
