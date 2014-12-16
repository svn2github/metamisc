\name{adjust.pm}
\alias{adjust.pm}

\title{Update a prediction model}

\description{
Update a prediction model.
}
\usage{
\method{adjust}{pm}(x, data, method="is", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class "\code{\link{pm}}".}
  \item{data}{A data frame with information on predictor and outcome variables of the prediction model.}
  \item{method}{Character indicating the updating method. Choose "\code{i}" for an intercept update, or "\code{is}" for logistic calibration (update intercept + common slope).}
  \item{\ldots}{Other arguments (currently not used).}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{An object of class "\code{\link{pm}}".
}
\references{
%% ~put references to the literature/web site here ~
}

\author{Thomas Debray <thomas.debray@gmail.com>}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(DVTipd)
data(DVTmodels)
weights.Wells <- DVTmodels$weights["Wells",]
coefs <- weights.Wells[!is.na(weights.Wells)]
fmla <- as.formula(paste("dvt~0+",paste(names(coefs), collapse="+")))
model.Wells <- pm(coefs, fmla, binomial()) 

#Update common slope and add intercept term to Wells model
update(model.Wells, data=DVTipd, method="is")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
