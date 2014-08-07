\name{DVTipd}
\alias{DVTipd}
\docType{data}
\title{Hypothetical dataset for diagnosis of Deep Vein Thrombosis (DVT)
}
\description{
A hypothetical dataset with 500 subjects suspected of having deep vein thrombosis (DVT).
}
\usage{data(DVTipd)}
\format{
  A data frame with 500 observations of 15 variables.
  \describe{
    \item{\code{sex}}{gender (0=female, 1=male)}
    \item{\code{malign}}{active malignancy (0=no active malignancy, 1=active malignancy)}
    \item{\code{par}}{paresis (0=no paresis, 1=paresis)}
    \item{\code{surg}}{recent surgery or bedridden}
    \item{\code{tend}}{tenderness venous system}
    \item{\code{oachst}}{oral contraceptives or hst}
    \item{\code{leg}}{entire leg swollen}
    \item{\code{notraum}}{absence of leg trauma}
    \item{\code{calfdif3}}{calf difference >= 3 cm}
    \item{\code{pit}}{pitting edema}
    \item{\code{vein}}{vein distension}
    \item{\code{altdiagn}}{alternative diagnosis present}
    \item{\code{histdvt}}{history of previous DVT}
    \item{\code{ddimdich}}{dichotimized D-dimer value}
    \item{\code{dvt}}{final diagnosis of DVT}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(DVTipd)
str(DVTipd) 
summary(apply(DVTipd,2,as.factor))

## Develop a prediction model to predict presence of DVT
model.dvt <- glm("dvt~sex+oachst+malign+surg+notraum+vein+calfdif3+ddimdich", 
                  family=binomial, data=DVTipd)
summary(model.dvt)
}
\keyword{datasets}
