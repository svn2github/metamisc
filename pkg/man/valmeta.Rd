\name{valmeta}
\alias{valmeta}
\alias{valmeta.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Meta-analysis of prediction model performance
}
\description{
This function allows to meta-analyze the performance of a prediction model, and to obtain summary estimates of the c-statistic, the total observed-expected ratio and the calibration slope. Where appropriate, data transformations are applied and missing information is derived from available quantities.
}
\usage{
valmeta(x, ...)
valmeta.default(cstat, cstat.se, cstat.95CI,
                N, O, method="REML", knha=TRUE, verbose=FALSE, 
                method.restore.se="Newcombe.4", scale = "logit", 
                ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{cstat}{vector with the estimated c-statistic for each valiation}
  \item{cstat.se}{vector with the standard error of the estimated c-statistics}
  \item{cstat.95CI}{2-dimensional array with the lower (first column) and upper (second column) boundary of the 95\% confidence interval of the estimated c-statistics}
  \item{N}{vector with the total number of participants for each valiation}
  \item{O}{vector with the total number of observed events for each valiation}
  \item{method}{character string specifying whether a fixed- or a random-effects model should be fitted. A fixed-effects model is fitted when using \code{method="FE"}. Random-effects models are fitted by setting method equal to one of the following: \code{"DL"}, \code{"HE"}, \code{"SJ"}, \code{"ML"}, \code{"REML"}, \code{"EB"}, \code{"HS"}, \code{"GENQ"} or \code{"BAYES"}. Default is \code{"REML"}. See 'Details'.}
  \item{knha}{logical specifying whether the method by Knapp and Hartung (2003) should be used for adjusting test statistics and confidence intervals (default is TRUE). Type '\code{?rma}' for more details.}
  \item{verbose}{if TRUE then messages generated during the fitting process will be displayed.}
  \item{method.restore.c.se}{character string specifying which method should be used to restore missing estimates for the standard error of the c-statistic. So far, only method 2 (\code{"Newcombe.2"}) and method 4 (\code{"Newcombe.4"}) are supported. These methods have been described by Newcombe in 2006.}
  \item{scale.c}{character string specifying whether the meta-analysis should use transformed values of the c-statistic and the corresponding standard error. Options are \code{scale="logit"} (logistic transformation) or \code{"identity"} (no transformation).}
}
\details{
\bold{Meta-analysis of the c-statistic}

The c-statistic is a measure of discrimination, and indicates the ability of a prediction model to distinguish between patients developing and not developing the outcome. The c-statistic typically ranges from 0.5 (no discriminative ability) to 1 (perfect discriminative ability). Because the c-statistic is highly dependent on patient spectrum (case-mix), between-study heterogeneity likely to occur in a meta-analysis. For this reason, it is generally recommended to adopt random-effects models (Debray, BMJ 2016). By default, extracted c-statistics are rescaled using the logistic transformation to improve the validity of the Normality assumptions in the meta-analysis models. 

\bold{Meta-analysis models}

Most meta-analysis models are implemented using \emph{metafor}, and can directly be specified using the parameter \code{method}. More information on these can be obtained by typing '\code{?rma}'.
}
\value{
\item{cstat}{object of class \code{vmasum} containing the meta-analysis results of the c-statistic.}
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
Debray TPA, Damen JAAG, Snell KIE, Ensor J, Hooft L, Reitsma JB, et al. A guide to systematic review and meta-analysis of prediction model performance. \emph{BMJ}. 2017;356:i6460. 

Hanley JA, McNeil BJ. The meaning and use of the area under a receiver operating characteristic (ROC) curve. \emph{Radiology}. 1982 Apr;143(1):29--36. 

Newcombe RG. Confidence intervals for an effect size measure based on the Mann-Whitney statistic. Part 2: asymptotic methods and evaluation. \emph{Stat Med}. 2006 Feb 28;25(4):559--73. 

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
data(EuroSCORE)
attach(EuroSCORE)

#Execute a random effects meta-analysis on the logit scale
fit <- valmeta(cstat=c.index, cstat.se=se.c.index, N=n, O=n.events)
fit

#Detach the dataset
detach(EuroSCORE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
