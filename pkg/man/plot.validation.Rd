\name{plot.validation}
\alias{plot.validation}

\title{Plot the results of a validation study}

\description{
This function plots the ROC curve or calibration curve of a validated prediction model.
}
\usage{
plot.validation(x, type="discrimination", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class "\code{pm}"}
  \item{type}{A character indicating whether to plot the ROC curve ("discrimination") or calibration curve ("calibration").}
  \item{...}{Additional arguments for \code{\link{plot}}.}
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

\author{Thomas Debray <thomas.debray@gmail.com>}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Apparent validation of a novel prediction model
data(DVTipd)
model.dvt <- glm("dvt~sex+oachst+malign+surg+notraum+vein+calfdif3+ddimdich", 
                  family=binomial, data=DVTipd)
val.dvt <- validation(pm(model.dvt), DVTipd)
plot(val.dvt, type="discrimination")
plot(val.dvt, type="calibration")

## External validation of a literature model
data(DVTmodels)
weights.Oudega <- DVTmodels$weights["Oudega",]
coefs <- weights.Oudega[!is.na(weights.Oudega)]
fmla <- as.formula(paste("dvt~",paste(names(coefs)[-1], collapse="+")))
model.Oudega <- pm(coefs, fmla, binomial()) 
val.Oudega <- validation(model.Oudega, DVTipd)
plot(val.Oudega, type="calibration")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
