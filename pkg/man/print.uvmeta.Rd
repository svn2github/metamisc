\name{print.uvmeta}
\alias{print.uvmeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.uvmeta(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
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
function (x, ...) 
{
    cat("Call:\n")
    print(x$call)
    cat("\nCoefficients:\n\n")
    fe = array(NA, dim = c(2, 3))
    fe[1, ] = cbind(round(x$fixef$mean, 5), round(sqrt(x$fixef$var), 
        5), x$fixef$p2t)
    fe[2, ] = cbind(round(x$ranef$mean, 5), round(sqrt(x$ranef$var), 
        5), x$ranef$p2t)
    colnames(fe) = c("Estimate", "StdErr", "p.value (2-tailed)")
    rownames(fe) = c("Fixed Effects", "Random Effects")
    print(fe, quote = F)
    cat(paste("\nTau squared: \t\t", round(x$tau_sq, 5), sep = ""))
    cat(paste("\nCochran's Q statistic: \t", round(x$Q$Q, 5), 
        " (p.value: ", round(x$Q$p.value, 5), ")", sep = ""))
    cat(paste("\nH-square index: \t", round(x$H2$H2, 5), " (95\% CI: ", 
        round(x$H2$H2.lowerconf, 5), "; ", round(x$H2$H2.upperconf, 
            5), ")", sep = ""))
    cat(paste("\nI-square index: \t", round(x$I2$I2 * 100, 5), 
        "\% (95\% CI: ", round(x$I2$I2.lowerconf * 100, 5), "\% ; ", 
        round(x$I2$I2.upperconf * 100, 5), "\%)\n", sep = ""))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
