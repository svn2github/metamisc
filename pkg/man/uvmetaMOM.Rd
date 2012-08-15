\name{uvmetaMOM}
\alias{uvmetaMOM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
uvmetaMOM(r, vars)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{vars}{
%%     ~~Describe \code{vars} here~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (r, vars) 
{
    w = 1/vars
    weighted_Tbar = sum(r * w)/sum(w)
    var_T = 1/sum(w)
    se_T = sqrt(var_T)
    z_T = weighted_Tbar/se_T
    Q = sum(w * (r - weighted_Tbar)^2)
    I_sq = 0
    H_sq = 0
    se_lnH = 0
    dfr = length(r) - 1
    if (Q > dfr) {
        re_C = sum(w) - sum(w^2)/sum(w)
        between_study_var = (Q - dfr)/re_C
        I_sq = (Q - dfr)/Q
        H_sq = Q/dfr
        se_lnH = (log(Q) - log(dfr))/(2 * (sqrt(2 * Q) - sqrt((2 * 
            length(r)) - 3)))
    }
    else {
        between_study_var = 0
        se_lnH = sqrt((1/(2 * (length(r) - 2))) * (1 - (1/(3 * 
            ((length(r) - 2)^2)))))
    }
    re_v = vars + between_study_var
    re_w = 1/re_v
    re_weighted_Tbar = sum(r * re_w)/sum(re_w)
    re_var_T = 1/sum(re_w)
    re_se_T = sqrt(re_var_T)
    re_z_T = re_weighted_Tbar/re_se_T
    fe_p1t = pnorm(weighted_Tbar/sqrt(var_T), lower = F)
    fe_p2t = 2 * pnorm(weighted_Tbar/sqrt(var_T), lower = F)
    fe.lowerconf = weighted_Tbar + qnorm(0.025) * sqrt(var_T)
    fe.upperconf = weighted_Tbar + qnorm(0.975) * sqrt(var_T)
    fixef.results = list(mean = weighted_Tbar, var = var_T, p2t = fe_p2t, 
        lowerconf = fe.lowerconf, upperconf = fe.upperconf)
    re_p1t = pnorm(re_weighted_Tbar/sqrt(re_var_T), lower = F)
    re_p2t = 2 * pnorm(re_weighted_Tbar/sqrt(re_var_T), lower = F)
    re.lowerconf = re_weighted_Tbar + qnorm(0.025) * sqrt(re_var_T)
    re.upperconf = re_weighted_Tbar + qnorm(0.975) * sqrt(re_var_T)
    ranef.results = list(mean = re_weighted_Tbar, var = re_var_T, 
        p2t = re_p2t, lowerconf = re.lowerconf, upperconf = re.upperconf)
    lnH = log(sqrt(H_sq))
    lnH_p2t = 2 * pnorm(lnH/se_lnH, lower = F)
    lnH.lowerconf = lnH + qnorm(0.025) * se_lnH
    lnH.upperconf = lnH + qnorm(0.975) * se_lnH
    H2.lowerconf = (exp(lnH.lowerconf))^2
    H2.upperconf = (exp(lnH.upperconf))^2
    H2.results = list(H2 = H_sq, lnH = lnH, lnH_p2t = lnH_p2t, 
        lnH.lowerconf = lnH.lowerconf, lnH.upperconf = lnH.upperconf, 
        H2.lowerconf = H2.lowerconf, H2.upperconf = H2.upperconf)
    I2.lowerconf = (H2.lowerconf - 1)/H2.lowerconf
    I2.upperconf = (H2.upperconf - 1)/H2.upperconf
    I2.results = list(I2 = I_sq, I2.lowerconf = I2.lowerconf, 
        I2.upperconf = I2.upperconf)
    Q.critical = qchisq(0.95, df = (length(r) - 1))
    Q_p = 1 - pchisq(Q, df = (length(r) - 1))
    Q.results = list(Q = Q, critical = Q.critical, p.value = Q_p)
    out <- list(fixef = fixef.results, ranef = ranef.results, 
        Q = Q.results, H2 = H2.results, I2 = I2.results, tau_sq = between_study_var)
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
