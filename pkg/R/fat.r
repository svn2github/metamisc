### To add / change:
# 
# 
#' Regression tests for detecting funnel plot asymmetry
#'
#' Formal tests to assess the presence of funnel plot asymmetry typically estimate the association between the reported effect size and their standard error, the total sample size or the inverse of the total sample size
#' 
#' @param b Vector with the effect size of each study. Examples are log odds ratio, log hazards ratio, log relative risk. 
#' @param b.se vector with the standard error of the effect size of each study
#' @param method Method for testing funnel plot asymmetry, defaults to \code{"E-UW"} (Egger's test). More info in "Details".
#'
#' @details A common method to test the presence of small-study effects (\code{method="E-UW"}) is given as the following unweighted regression model: \deqn{\hat{b}_k = \beta_0 + \beta_1\, \widehat \mathrm{SE}(\hat{b}_k) + \epsilon_k \;,\; \epsilon_k \sim \mathcal{N}\left(0, \sigma^2 \right) }{b=A+B*b.se + e;  e~N(0,s^2)}
#' Whereas \eqn{\beta_0}{A} indicates the size and direction of the treatment effect, \eqn{\beta_1}{B} provides a measure of asymmetry; the larger its deviation from zero the more pronounced the asymmetry.
#' @return a list containing the following entries:
#' \itemize{
##'  \item{"pval"} {P-value indicating statistical significance of funnel plot asymettry test}
##'  \item{"model"}{A fitted \code{glm} object, representing the estimated regression model}
##' }
#' @author Thomas Debray
#' 
#' @references Debray TPA, Moons KGM, Riley RD. Detecting small-study effects and funnel plot asymmetry in meta-analysis of survival data: a comparison of new and existing tests. \emph{Res Syn Meth}. 2017 Oct 3.
#'
#' @examples 
#' data(Fibrinogen)
#' b <- log(Fibrinogen$HR)
#' b.se <- ((log(Fibrinogen$HR.975) - log(Fibrinogen$HR.025))/(2*qnorm(0.975)))
#' fat(b=b, b.se=b.se)

#'
#' @import stats
#' @importFrom stats pt
#' @importFrom stats qnorm
#' 
#' @exportMethod fat
#' 
#'
fat <- function(b, b.se, method="E-UW") 
{
  
  if (method %in% c("E-UW", "E-FIV")) {
    if (length(b) != length(b.se))
      stop("Incompatible vector sizes for 'b' and 'b.se'!")
    
    # Identify which studies can be used
    nstudies <- sum(!is.na(b) & !is.na(b.se))
  }
  
  
  w <- NA
  
  if (method=="E-FIV") {
    suppressWarnings(m.fat <- try(glm(b~b.se, weights=w), silent=T))
  } else if (method=="E-UW")  {
    suppressWarnings(m.fat <- try(glm(b~b.se), silent=T))
  } else {
    stop("Invalid method for Egger's test")
  }
  
  if ("try-error" %in% attr(m.fat,"class")) {
    warning("Estimation of the regression model unsuccessful, results omitted.")
    z.fat <- NA
    p.fat <- NA
    m.fat <- NA
  } else {
    z.fat <- coefficients(m.fat)[2]/sqrt(diag(vcov(m.fat))[2])
    p.fat <- 2*pt(-abs(z.fat),df=(nstudies-2))
  }

  out <- list(pval = p.fat, model = m.fat)
  return(out)
}