\name{predict.riley}
\alias{predict.riley}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
predict.riley(object, level = 0.95, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{level}{
%%     ~~Describe \code{level} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, level = 0.95, ...) 
{
    alpha = (1 - level)/2
    predint = array(NA, dim = c(2, 3))
    rownames(predint) = c("Sens", "FPR")
    colnames(predint) = c("Estimate", paste((alpha * 100), "\%"), 
        paste(((1 - alpha) * 100), "\%"))
    df = object$df
    if (df > 0) {
        predint[1, ] = inv.logit(c(coefficients(object)["mu1"], 
            (qt(alpha, df = df) * sqrt((coefficients(object)["psi1"]^2) + 
                diag(vcov(object))[1]) + coefficients(object)["mu1"]), 
            (qt((1 - alpha), df = df) * sqrt((coefficients(object)["psi1"]^2) + 
                diag(vcov(object))[1]) + coefficients(object)["mu1"])))
        predint[2, ] = inv.logit(c(coefficients(object)["mu2"], 
            (qt(alpha, df = df) * sqrt((coefficients(object)["psi2"]^2) + 
                diag(vcov(object))[2]) + coefficients(object)["mu2"]), 
            (qt((1 - alpha), df = df) * sqrt((coefficients(object)["psi2"]^2) + 
                diag(vcov(object))[2]) + coefficients(object)["mu2"])))
    }
    else {
        predint[1, ] = c(inv.logit(mu1), 0, 1)
        predint[2, ] = c(inv.logit(mu2), 0, 1)
    }
    predint
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
