\name{riley.default}
\alias{riley.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
riley.default(X = NULL, TP, FN, FP, TN, correction = 0.5, correction.control = "all", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{TP}{
%%     ~~Describe \code{TP} here~~
}
  \item{FN}{
%%     ~~Describe \code{FN} here~~
}
  \item{FP}{
%%     ~~Describe \code{FP} here~~
}
  \item{TN}{
%%     ~~Describe \code{TN} here~~
}
  \item{correction}{
%%     ~~Describe \code{correction} here~~
}
  \item{correction.control}{
%%     ~~Describe \code{correction.control} here~~
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
function (X = NULL, TP, FN, FP, TN, correction = 0.5, correction.control = "all", 
    ...) 
{
    if (!is.null(X)) {
        X <- as.data.frame(X)
        origdata <- X
        TP <- X$TP
        FN <- X$FN
        FP <- X$FP
        TN <- X$TN
    }
    if (correction.control == "all") {
        if (any(c(TP, FN, FP, TN) == 0)) {
            TP <- TP + correction
            FN <- FN + correction
            FP <- FP + correction
            TN <- TN + correction
        }
    }
    if (correction.control == "single") {
        correction = ((((TP == 0) | (FN == 0)) | (FP == 0)) | 
            (TN == 0)) * correction
        TP <- correction + TP
        FN <- correction + FN
        FP <- correction + FP
        TN <- correction + TN
    }
    numstudies = length(TP)
    df <- numstudies - 5
    if (df < 0) {
        warning("There are very few primary studies!")
    }
    number.of.pos <- TP + FN
    number.of.neg <- FP + TN
    sens <- TP/number.of.pos
    fpr <- FP/number.of.neg
    var.sens = sens * (1 - sens)/number.of.pos
    var.fpr = fpr * (1 - fpr)/number.of.neg
    logit.sens <- logit(sens)
    logit.fpr <- logit(fpr)
    var.logit.sens <- 1/(sens * (1 - sens) * number.of.pos)
    var.logit.fpr <- 1/(fpr * (1 - fpr) * number.of.neg)
    vars = cbind(var.logit.sens, var.logit.fpr)
    Y = array(NA, dim = c((length(logit.sens) * 2), 1))
    for (i in 1:length(logit.sens)) {
        Y[((i - 1) * 2 + 1)] = logit.sens[i]
        Y[((i - 1) * 2 + 2)] = logit.fpr[i]
    }
    sumlsens <- uvmeta(r = logit.sens, v = var.logit.sens, method = "MOM")
    sumlfpr <- uvmeta(r = logit.fpr, v = var.logit.fpr, method = "MOM")
    negfullloglik <- function(pars, Y, vars) {
        beta1 = pars[1]
        beta2 = pars[2]
        psisq1 = pars[3]^2
        psisq2 = pars[4]^2
        rho = inv.logit(pars[5]) * 2 - 1
        k = 2
        n = dim(Y)[1]/2
        Beta = rbind(beta1, beta2)
        X = array(NA, dim = c(n * 2, 2))
        X[, 1] = rep(c(1, 0), n)
        X[, 2] = rep(c(0, 1), n)
        Phi = array(0, dim = c((n * 2), (n * 2)))
        for (i in 1:n) {
            Phi[((i - 1) * 2 + 1), ((i - 1) * 2 + 1)] = vars[i, 
                1] + psisq1
            Phi[((i - 1) * 2 + 2), ((i - 1) * 2 + 2)] = vars[i, 
                2] + psisq2
            Phi[((i - 1) * 2 + 1), ((i - 1) * 2 + 2)] = rho * 
                sqrt((vars[i, 1] + psisq1) * (vars[i, 2] + psisq2))
            Phi[((i - 1) * 2 + 2), ((i - 1) * 2 + 1)] = rho * 
                sqrt((vars[i, 1] + psisq1) * (vars[i, 2] + psisq2))
        }
        0.5 * ((n - k) * log(2 * pi) - log(det(t(X) \%*\% X)) + 
            log(det(Phi)) + log(det(t(X) \%*\% solve(Phi) \%*\% X)) + 
            (t(Y - X \%*\% Beta) \%*\% solve(Phi) \%*\% (Y - X \%*\% 
                Beta)))
    }
    pars.start = c(sumlsens$ranef$mean, sumlfpr$ranef$mean, sqrt(sumlsens$ranef$var), 
        sqrt(sumlfpr$ranef$var), 0)
    fit = optim(pars.start, negfullloglik, Y = Y, vars = vars, 
        hessian = T)
    mu1 = fit$par[1]
    mu2 = fit$par[2]
    psi1 = abs(fit$par[3])
    psi2 = abs(fit$par[4])
    rhoT = fit$par[5]
    coefficients = c(mu1, mu2, psi1, psi2, rhoT)
    names(coefficients) = c("mu1", "mu2", "psi1", "psi2", "rhoT")
    Sigma = solve(fit$hessian)
    colnames(Sigma) = c("mu1", "mu2", "psi1", "psi2", "rhoT")
    rownames(Sigma) = c("mu1", "mu2", "psi1", "psi2", "rhoT")
    iterations <- fit$iterations
    logLik <- fit$value
    output <- list(coefficients = coefficients, vcov = Sigma, 
        df = df, nobs = 2 * numstudies, logLik = logLik, iterations = (iterations + 
            1), call = match.call(), data = origdata)
    class(output) <- "riley"
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
