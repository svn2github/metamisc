\name{validate}
\alias{validate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
validate(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
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
data(DVTipd)

## Develop a prediction model to predict presence of DVT
model.dvt <- glm("dvt~sex+oachst+malign+surg+notraum+vein+calfdif3+ddimdich", family=binomial, data=DVTipd)
val.dvt <- validate(model.dvt, DVTipd)
val.dvt

## Example for a model without an intercept term
data(DVTmodels)
weights.Wells <- DVTmodels$weights["Wells",]
coefs <- weights.Wells[!is.na(weights.Wells)]
fmla <- as.formula(paste("dvt~0+",paste(names(coefs)[-1], collapse="+")))
model.Wells <- as.model(coefs, fmla, binomial()) 
val.Wells <- validate(model.Wells, DVTipd)
val.Wells

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
