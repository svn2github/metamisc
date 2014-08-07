\name{pm}
\alias{pm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create a prediction model object
%%  ~~function to do ... ~~
}
\description{
Create an object of the class \code{pm}, for which several functions are available.
}
\usage{
pm(x, formula, family=binomial())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class "\code{\link{glm}}" or a named vector of regression coefficients.}
  \item{formula}{An object of class "\code{\link{formula}}" providing a symbolic description of the fitted model.}
  \item{family}{An object of class "\code{\link{family}}" describing the error distribution and link function to be used in the model.}
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
model.dvt <- glm("dvt~sex+oachst+malign+surg+notraum+vein+calfdif3+ddimdich", 
                  family=binomial, data=DVTipd)
pm(model.dvt)                  
                  
data(DVTmodels)

# extract coefficients Oudega model
weights.Oudega <- DVTmodels$weights["Oudega",]
coefs <- weights.Oudega[!is.na(weights.Oudega)]
fmla <- as.formula(paste("dvt~",paste(names(coefs)[-1], collapse="+")))

# create pm object and validate the model in external data
model.Oudega <- pm(coefs, fmla, binomial())
val.Oudega <- validate(model.Oudega, DVTipd)
val.Oudega
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
