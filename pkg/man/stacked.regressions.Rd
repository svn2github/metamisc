\name{stacked.regressions}

\title{Meta-analysis of prediction models (stacked regression).}

\description{This function performs performs a meta-analysis of multiple prediction models (Debray 2014). }

\usage{
stacked.regressions <- function (models=list(), names.models=NULL, outcome=NA, data=NULL) 
       }

\arguments{
  \item{models}{A list of prediction models. Each element may represent a named list of regression coefficients, or a \code{glm} fit.}
  \item{names.models}{A vector with the names of the provided prediction models.}
  \item{outcome}{The name of the outcome variable of interest.}
  \item{data}{A data frame containing the data used for meta-analyzing the prediction models.}
}

\details{  
}

\value{ An object of the class \code{metamodel} for which many standard methods are available.}

\references{
Debray TPA, Koffijberg H, Nieboer D, Vergouwe Y, Steyerberg EW, Moons KGM. Meta-analysis and aggregation of multiple published prediction models. \emph{Statistics in Medicine} 2014; Early View.
}

\author{ Thomas Debray <thomas.debray@gmail.com> }

\seealso{

}

\examples {
  #Load example
  data(DVTmodels)
  models <- list(Wells=DVTmodels[1,], ModifiedWells=DVTmodels[2,], Gagne=DVTmodels[3,], Hamilton=DVTmodels[4,], Oudega=DVTmodels[5,])
  outcome <- "dvt"
  
  #Generate a dataset
  
  #Apply meta-analysis

}

\keyword{ meta-analysis }
