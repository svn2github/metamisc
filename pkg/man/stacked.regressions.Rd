\name{stacked.regressions}

\title{Meta-analysis of risk prediction models: Stacked Regression.}

\description{This function performs performs a meta-analysis of multiple risk prediction models by applying Stacked Regressions (Debray 2014). }

\usage{
stacked.regressions <- function (models=list(), names.models=NULL, outcome=NA, data=NULL) 
       }

\arguments{
  \item{models}{A list of prediction models. Each element may represent a named vector of regression coefficients, or a \code{glm} object. Note that the coefficients (e.g. intercept terms or interaction terms) should be named according to the output of a \code{glm} analysis.}
  \item{names.models}{A vector with the names of the provided prediction models.}
  \item{outcome}{The name of the outcome variable of interest.}
  \item{data}{A data frame containing the data used for meta-analyzing the prediction models.}
}


\value{ An object of the class \code{metamodel} for which many standard methods are available.}

\references{
Debray TPA, Koffijberg H, Nieboer D, Vergouwe Y, Steyerberg EW, Moons KGM. Meta-analysis and aggregation of multiple published prediction models. \emph{Statistics in Medicine} 2014; Early View.
}

\author{ Thomas Debray <thomas.debray@gmail.com> }


\examples {
  #Load example
  data(DVTmodels)
  models <- list(Wells=DVTmodels[1,], ModifiedWells=DVTmodels[2,], Gagne=DVTmodels[3,], Hamilton=DVTmodels[4,], Oudega=DVTmodels[5,])
  outcome <- "dvt"
  
  #Generate a dataset
  
  #Apply meta-analysis

}

\keyword{ meta-analysis }
