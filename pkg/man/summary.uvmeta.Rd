\name{summary.uvmeta}
\alias{summary.uvmeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summary.uvmeta(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
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
function (object, ...) 
{
    se <- c(sqrt(object$fixef$var), sqrt(object$ranef$var), NA, 
        NA, NA)
    lower.conf = c(object$fixef$lowerconf, object$ranef$lowerconf, 
        NA, object$H2$H2.lowerconf, object$I2$I2.lowerconf)
    upper.conf = c(object$fixef$upperconf, object$ranef$upperconf, 
        NA, object$H2$H2.upperconf, object$I2$I2.upperconf)
    TAB = cbind(Estimate = c(object$fixef$mean, object$ranef$mean, 
        object$Q$Q, object$H2$H2, object$I2$I2), StdErr = se, 
        lower.conf = lower.conf, upper.conf = upper.conf)
    rownames(TAB) = c("Fixed Effects", "Random Effects", "Cochran Q", 
        "H-square index", "I-square index")
    res = list(call = object$call, estimates = TAB)
    class(res) = "summary.macc"
    res
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
