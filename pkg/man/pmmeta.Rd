\name{pmmeta}
\alias{pmmeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Meta-analysis of risk prediction models
}
\description{
This function performs performs a meta-analysis of multiple risk prediction models by applying Stacked Regressions or Model Averaging (Debray 2014).
}
\usage{
pmmeta(models=list(), names.models=NULL, outcome=NA, data=NULL,
       method="stacked.regressions")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{models}{A list of prediction models. Each element may represent a named vector of regression coefficients, or a \code{glm} object. Note that the coefficients (e.g. intercept terms or interaction terms) should be named according to the output of a \code{glm} analysis.}
  \item{names.models}{A vector with the names of the provided prediction models.}
  \item{outcome}{The name of the outcome variable of interest.}
  \item{data}{A data frame containing the data used for meta-analyzing the prediction models.}
  \item{method}{Method of meta-analysis: \code{stacked.regressions} or \code{model.averaging}}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
An object of the class \code{metamodel} for which many standard methods are available.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Debray TPA, Koffijberg H, Nieboer D, Vergouwe Y, Steyerberg EW, Moons KGM. Meta-analysis and aggregation of multiple published prediction models. \emph{Statistics in Medicine} 2014; Early View.
}
\author{
Thomas Debray <thomas.debray@gmail.com>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{DVTmodels}}
}
\examples{
## Load example
data(DVTmodels)

  
## Generate artificial validation dataset
  
  #Apply meta-analysis
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~meta-analysis }

