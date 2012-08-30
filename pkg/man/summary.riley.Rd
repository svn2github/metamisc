\name{summary.riley}
\alias{summary.riley}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parameter summaries of \code{\link{riley}} models.
}
\description{
Provides the summary estimtes of the alternative model for bivariate random-effects meta-analysis by Riley et al. (2008) with their corresponding confidence intervals. The model parameters are given as mu1 (logit sensitivity), mu2 (logit FPR), psi1, psi2 and rho; the summary sensitivity and false positive rate are added for completeness. Confidence intervals are derived from the inverse Hessian. 
}
\usage{
summary.riley(object, level = 0.95, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a \code{riley} object.
}
  \item{level}{
numeric, the level for calculations of confidence intervals
}
  \item{\dots}{
arguments to be passed on to other functions
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
Riley, RD., Thompson, JR., & Abrams, KR. (2008).
\dQuote{An alternative model for bivariate random-effects meta-analysis when the within-study correlations are unknown.} \emph{Biostatistics}, \bold{9}, 172--186.
}
\author{
Thomas Debray <thomas.debray@gmail.com>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{riley}}
\code{\link{plot.riley}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, level = 0.95, ...) 
{
    confints <- cbind(object$coefficients, confint(object, level = level))
    colnames(confints)[1] <- "Estimate"
    confints <- rbind(inv.logit(confints[1:2, ]), confints)
    rownames(confints)[1:2] <- c("Sens", "FPR")
    confints[7, ] = inv.logit(confints[7, ]) * 2 - 1
    rownames(confints)[7] = "rho"
    res <- list(call = object$call, confints = confints)
    class(res) <- "summary.riley"
    res
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{riley}
\keyword{confidence}
