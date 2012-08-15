\name{plot.riley}
\alias{plot.riley}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.riley(x, plotsumm = TRUE, plotnumerics = TRUE, level = 0.95, main = "", ylim = c(0, 1), xlim = c(0, 1), pch = 1, lty = 1, lwd = 1, cex.numerics = 0.45, add = F, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{plotsumm}{
%%     ~~Describe \code{plotsumm} here~~
}
  \item{plotnumerics}{
%%     ~~Describe \code{plotnumerics} here~~
}
  \item{level}{
%%     ~~Describe \code{level} here~~
}
  \item{main}{
%%     ~~Describe \code{main} here~~
}
  \item{ylim}{
%%     ~~Describe \code{ylim} here~~
}
  \item{xlim}{
%%     ~~Describe \code{xlim} here~~
}
  \item{pch}{
%%     ~~Describe \code{pch} here~~
}
  \item{lty}{
%%     ~~Describe \code{lty} here~~
}
  \item{lwd}{
%%     ~~Describe \code{lwd} here~~
}
  \item{cex.numerics}{
%%     ~~Describe \code{cex.numerics} here~~
}
  \item{add}{
%%     ~~Describe \code{add} here~~
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
function (x, plotsumm = TRUE, plotnumerics = TRUE, level = 0.95, 
    main = "", ylim = c(0, 1), xlim = c(0, 1), pch = 1, lty = 1, 
    lwd = 1, cex.numerics = 0.45, add = F, ...) 
{
    if (!add) 
        plot(-500, -500, type = "l", xlim = xlim, ylim = ylim, 
            xlab = "1-Specificity", ylab = "Sensitivity", main = main)
    alpha = (1 - level)/2
    FP <- x$data$FP
    negatives <- FP + x$data$TN
    FPR <- FP/negatives
    mu = x$coefficients[c("mu1", "mu2")]
    Sigma = x$vcov[c(1, 2), c(1, 2)]
    mu.ellipse <- ellipse(Sigma, centre = mu, level = level)
    polygon(inv.logit(mu.ellipse[, 2]), inv.logit(mu.ellipse[, 
        1]), lty = lty, lwd = lwd)
    if (plotsumm) 
        points(inv.logit(mu[2]), inv.logit(mu[1]), pch = pch)
    if (plotnumerics) {
        ci = summary(x, level = level)[2]$confints
        text(0.8, 0.15, labels = "Estimate", pos = 2, cex = cex.numerics)
        text(0.9, 0.15, labels = paste((alpha * 100), "\% CI", 
            sep = ""), pos = 2, cex = cex.numerics)
        text(1, 0.15, labels = paste(((1 - alpha) * 100), "\% CI", 
            sep = ""), pos = 2, cex = cex.numerics)
        text(0.5, 0.1, labels = "Sensitivity", pos = 4, cex = cex.numerics)
        text(0.8, 0.1, labels = paste("", formatC(round(ci["Sens", 
            1], 2), 2, format = "f", flag = "0")), pos = 2, cex = cex.numerics)
        text(0.9, 0.1, labels = paste("", formatC(round(ci["Sens", 
            2], 2), 2, format = "f", flag = "0")), pos = 2, cex = cex.numerics)
        text(1, 0.1, labels = paste("", formatC(round(ci["Sens", 
            3], 2), 2, format = "f", flag = "0")), pos = 2, cex = cex.numerics)
        text(0.5, 0.05, labels = "1-Specificity", pos = 4, cex = cex.numerics)
        text(0.8, 0.05, labels = paste("", formatC(round(ci["FPR", 
            1], 2), 2, format = "f", flag = "0")), pos = 2, cex = cex.numerics)
        text(0.9, 0.05, labels = paste("", formatC(round(ci["FPR", 
            2], 2), 2, format = "f", flag = "0")), pos = 2, cex = cex.numerics)
        text(1, 0.05, labels = paste("", formatC(round(ci["FPR", 
            3], 2), 2, format = "f", flag = "0")), pos = 2, cex = cex.numerics)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
