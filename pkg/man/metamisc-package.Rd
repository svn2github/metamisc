\name{metamisc-package}
\alias{metamisc-package}
\docType{package}
\title{Diagnostic and prognostic meta-analysis}
\description{This package provides functions for diagnostic and prognostic research. It allows to perform meta-analysis of univariate (e.g. added value of a prognostic marker), bivariate (e.g. diagnostic test accuracy) and multivariate (e.g. previously published prediction model) statistics. Functionalities to update and externally validate prediction models are also implemented.}

\details{
%\tabular{ll}{
%Package: \tab metamisc\cr
%Type: \tab Package\cr
%Version: \tab 0.1.2\cr
%Date: \tab 2014-12-16\cr
%License: \tab GPL-2 \cr
%}
The following functionality is currently implemented: univariate meta-analysis of summary data (\code{\link{uvmeta}}), bivariate meta-analysis of correlated outcomes (\code{\link{riley}}), bivariate meta-analysis of diagnostic test accuracy data (\code{reitsma}), meta-analysis of risk prediction models (\code{\link{pmmeta}}) and external validation of a prediction model (\code{\link{validation}}).}

\author{Thomas Debray <thomas.debray@gmail.com>}

\references{
Debray TPA, Koffijberg H, Vergouwe Y, Moons KGM, Steyerberg EW. Aggregating published prediction models with individual participant data: a comparison of different approaches. \emph{Statistics in Medicine} 2012; \bold{31}: 2697--2712. DOI: 10.1002/sim.5412.

Debray TPA, Koffijberg H, Nieboer D, Vergouwe Y, Steyerberg EW, Moons KGM. Meta-analysis and aggregation of multiple published prediction models. \emph{Statistics in Medicine} 2014; \bold{33}: 2341--2362. DOI: 10.1002/sim.6080.

DerSimonian R, Laird N. Meta-analysis in clinical trials. \emph{Controlled Clinical Trials} 1986; \bold{7}: 177--188.

Reitsma J, Glas A, Rutjes A, Scholten R, Bossuyt P, Zwinderman A. Bivariate analysis of sensitivity and specificity produces informative summary
measures in diagnostic reviews. \emph{Journal of Clinical Epidemiology} 2005; \bold{58}: 982--990.

Riley RD, Thompson JR, Abrams KR. An alternative model for bivariate random-effects meta-analysis when the within-study correlations are unknown. \emph{Biostatistics} 2008; \bold{9}: 172--186.
}

\keyword{package}
\keyword{regression}
\keyword{meta-analysis}
\keyword{validation}

