\name{valmeta}
\alias{valmeta}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Meta-analysis of prediction model performance
}
\description{
This function allows to meta-analyze the performance of a prediction model, and to obtain summary estimates of the c-statistic, the total observed-expected ratio and the calibration slope. Where appropriate, data transformations are applied and missing information is derived from available quantities.
}
\usage{
valmeta(cstat, cstat.se, cstat.95CI, OE, OE.95CI, N, O, E, method="REML", 
        knha=TRUE, verbose=FALSE, method.restore.c.se="Newcombe.4", 
        scale.c = "logit", scale.oe = "log", n.chains = 4, pars, \ldots)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cstat}{vector with the estimated c-statistic for each valiation}
  \item{cstat.se}{vector with the standard error of the estimated c-statistics}
  \item{cstat.95CI}{2-dimensional array with the lower (first column) and upper (second column) boundary of the 95\% confidence interval of the estimated c-statistics}
  \item{OE}{vector with the estimated ratio of total observed versus total expected events}
  \item{OE.95CI}{2-dimensional array with the lower (first column) and upper (second column) boundary of the 95\% confidence interval of the total O:E ratio}
  \item{N}{vector with the total number of participants for each valiation}
  \item{O}{vector with the total number of observed events for each valiation}
  \item{E}{vector with the total number of expected events for each valiation}
  \item{method}{character string specifying whether a fixed- or a random-effects model should be fitted. A fixed-effects model is fitted when using \code{method="FE"}. Random-effects models are fitted by setting method equal to one of the following: \code{"REML"} (Default), \code{"DL"}, \code{"HE"}, \code{"SJ"}, \code{"ML"}, \code{"EB"}, \code{"HS"}, \code{"GENQ"} or \code{"BAYES"}. See 'Details'.}
  \item{knha}{logical specifying whether the method by Knapp and Hartung (2003) should be used for adjusting test statistics and confidence intervals (default is TRUE). Type '\code{?rma}' for more details.}
  \item{verbose}{if TRUE then messages generated during the fitting process will be displayed.}
  \item{method.restore.c.se}{character string specifying which method should be used to restore missing estimates for the standard error of the c-statistic. So far, only method 2 (\code{"Newcombe.2"}) and method 4 (\code{"Newcombe.4"}) are supported. These methods have been described by Newcombe in 2006.}
  \item{scale.c}{character string specifying whether the meta-analysis should use transformed values of the c-statistic and the corresponding standard error. Options are \code{scale="logit"} (logistic transformation) or \code{"identity"} (no transformation).}
  \item{scale.oe}{character string specifying whether the meta-analysis should use transformed values of the total observed:expected ratio and the corresponding standard error. Options are \code{scale="log"} (log transformation) or \code{"identity"} (no transformation).}
  \item{n.chains}{the number of chains to use in the Gibbs sampler (\code{method="BAYES"}). More chains will improve the sensitivity of the convergence diagnostic, but will cause the simulation to run more slowly. The default number of chains is 4.}
  \item{pars}{A list with additional arguments.  The following parameters configure the MCMC sampling procedure:  \code{hp.mu.mean} (Hyperparameter: mean of the prior distribution of the random effects model, defaults to zero), \code{hp.mu.var} (Hyperparameter: variance of the prior distribution of the random effects model, defaults to 1E6), \code{hp.tau.min} (Hyperparameter: minimum value for the between-study standard deviation, defaults to 0), \code{hp.tau.max} (Hyperparameter: maximum value for the between-study standard deviation, defaults to 2)}
  \item{\ldots}{Additional arguments that are passed to \code{rma} or \code{runjags} (if \code{method="BAYES"}).}
}
\details{
Unless specified otherwise, all meta-analysis models assume random effects and are fitted  using restricted maximum likelihood estimation with the \code{metafor} package (Viechtbauer 2010).  Further, confidence intervals for the average performance are based on the Hartung-Knapp-Sidik-Jonkman method, to better account for the uncertainty in the estimated between-study heterogeneity (Debray 2016). A Bayesian meta-analysis can be performed by specifying \code{method="BAYES"}. In that case, the R packages \code{runjags} and \code{rjags} must be installed. Details on the prior distributions and assumptions are provided below.


\bold{Meta-analysis of the concordance statistic (c-statistic)}

The c-statistic is a measure of discrimination, and indicates the ability of a prediction model to distinguish between patients developing and not developing the outcome. The c-statistic typically ranges from 0.5 (no discriminative ability) to 1 (perfect discriminative ability). A meta-analysis for the c-statistic will be performed if the c-statistics (\code{cstat}) and their respective standard errors (\code{cstat.se}) are defined. For studies where the standard error is unknown, it can be derived from the 95\% confidence interval, or from \code{cstat}, \code{O} and \code{N} (Newcombe 2006). By default, c-statistics are rescaled using the logit transformation to improve the validity of Normality assumptions.

\bold{Mta-analysis of the total observed versus expected ratio (O:E ratio)}

The total O:E ratio provides a rough indication of the overall model calibration (across the entire range of predicted risks). A meta-analysis for the O:E ratio will be performed if the total number of observed (\code{O}) and expected (\code{E}) events are defined, or if their ratio (\code{OE}) is provided. If the standard error of the O:E ratio is not provided, it can be approximated from the 95\% confidence interval. Alternatively, if the number of expected events is treated as a fixed constant, the standard error can be derived from the error variance of the binomial or Poisson (if \code{O} << \code{N}) distribution. In that case, continuitiy corrections are applied for studies where \code{E=0} (if \code{scale.oe = "identity"}) or \code{O=0} (if \code{scale.oe = "log"}) by adding 0.5 to \code{O}, \code{E} and \code{N}. Standard errors are always prioritized as source of within-study variability in case multiple sources are available (e.g. from confidence intervals).   By default, O:E ratios are rescaled using the (natural) log transformation to improve the validity of Normality assumptions.


\bold{Bayesian meta-analysis}

For Bayesian meta-analysis models that involve the Gibbs sampler (\code{method="BAYES"}), the R packages \code{runjags} and \code{rjags} must be installed.  
}
\value{
An object of class \code{valmeta} with the following elements:
\item{cstat}{object of class \code{vmasum} containing the meta-analysis results of the c-statistic.}
\item{oe}{object of class \code{vmasum} containing the meta-analysis results of the total O:E ratio.}
\item{method}{character string specifying whether a fixed- or a random-effects model was fitted.}
%\item{results}{vector with the summary c-statistic, the corresponding standard error, 95\% confidence interval and 95\% prediction interval.}
%\item{scale}{the scale used for meta-analysis.}
%\item{data}{2-dimensional array with the effect estimates (\code{theta}) and (restored) error variances (\code{theta.var}) of each validation. For \code{method="BAYES"}, the estimates represent draws from the posterior distribution.}
%\item{rma}{if applicable, fitted \code{rma} object of the meta-analysis.}

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Debray TPA, Damen JAAG, Snell KIE, Ensor J, Hooft L, Reitsma JB, et al. A guide to systematic review and meta-analysis of prediction model performance. \emph{BMJ}. 2017; 356:i6460. 

Hanley JA, McNeil BJ. The meaning and use of the area under a receiver operating characteristic (ROC) curve. \emph{Radiology}. 1982; 143(1):29--36. 

Newcombe RG. Confidence intervals for an effect size measure based on the Mann-Whitney statistic. Part 2: asymptotic methods and evaluation. \emph{Stat Med}. 2006; 25(4):559--73. 

Viechtbauer W. Conducting Meta-Analyses in R with the metafor Package. \emph{Journal of Statistical Software}. 2010; 36(3). Available from: \url{http://www.jstatsoft.org/v36/i03/}
}
\author{Thomas Debray <thomas.debray@gmail.com>}

\seealso{
\code{\link{plot.valmeta}}
}

\examples{
data(EuroSCORE)
attach(EuroSCORE)

# Meta-analysis of the c-statistic (random effects)
fit <- valmeta(cstat=c.index, cstat.se=se.c.index, 
               cstat.95CI=cbind(c.index.95CIl,c.index.95CIu), 
               N=n, O=n.events)
fit
plot(fit)

# Meta-analysis of the total O:E ratio (random effects)
valmeta(O=n.events, E=e.events, N=n)     # Option 1
valmeta(O=n.events, E=e.events)          # Option 2 
valmeta(O=(Po/100)*n, E=(Pe/100)*n, N=n) # Option 3


# Meta-analysis of the c-statistic and the total O:E ratio (random effects)
valmeta(cstat=c.index, cstat.se=se.c.index, 
        cstat.95CI=cbind(c.index.95CIl,c.index.95CIu), 
        N=n, O=n.events, E=e.events)

\dontrun{
# Bayesian meta-analysis of the c-statistic (random effects)
fit2 <- valmeta(cstat=c.index, cstat.se=se.c.index, 
                cstat.95CI=cbind(c.index.95CIl,c.index.95CIu),
                N=n, O=n.events, method="BAYES")
plot(fit2)} 

#Detach the dataset
detach(EuroSCORE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ meta-analysis }
\keyword{ discrimination }
\keyword{ calibration }
