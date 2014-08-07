\name{validate}
\alias{validate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Validate a prediction model
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
validate(X, ds.ipd, time.calibration=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{An object of class "\code{\link{glm}}" or "\code{pmmodel}", representing a prediction model.}
  \item{ds.ipd}{A validation dataset with information on the predictor and outcome variables of \code{X}. }
  \item{time.calibration}{Not implemented at this stage.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{roc}{Model discrimination results, presented as an object of class "\code{roc}".}
  \item{cal}{Model calibration results ("estimate") with corresponding standard error ("se") and/or 95\% confidence interval ("2.5\%" and "97.5\%"). The following statistics are included: actual and predicted number of events ("$events"), ratio observed versus expected number of events ("$O:E"), calibration slope ("$slope") and calibration-in-the-large ("$intercept")}
  \item{predictions}{A data frame with the linear predictor ("lp"), predicted outcome ("yhat") and observed outcome ("y") for each subject from the validation dataset.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Thomas Debray <thomas.debray@gmail.com>}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(DVTipd)

## Develop a prediction model to predict presence of DVT
model.dvt <- glm("dvt~sex+oachst+malign+surg+notraum+vein+calfdif3+ddimdich", 
                  family=binomial, data=DVTipd)
val.dvt <- validate(model.dvt, DVTipd)
val.dvt
plot(val.dvt, type="calibration")

## Example for a model without an intercept term
data(DVTmodels)
weights.Wells <- DVTmodels$weights["Wells",]
coefs <- weights.Wells[!is.na(weights.Wells)]
fmla <- as.formula(paste("dvt~0+",paste(names(coefs)[-1], collapse="+")))
model.Wells <- pmmodel(coefs, fmla, binomial()) 
val.Wells <- validate(model.Wells, DVTipd)
val.Wells

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{validation}
\keyword{performance}
