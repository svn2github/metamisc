### To add / change:
# 
# 
#' Regression tests for detecting funnel plot asymmetry
#'
#' The presence of small-study effects is a common threat to systematic reviews and meta-analyses, especially when 
#' it is due to publication bias, which occurs when small primary studies are more likely to be reported (published) 
#' if their findings were positive. The presence of small-study effects can be verified by visual inspection of 
#' the funnel plot, where for each included study of the meta-analysis, the estimate of the reported effect size is 
#' depicted against a measure of precision or sample size. 
#' The premise is that the scatter of plots should reflect a funnel shape, if small-study 
#' effects do not exist. However, when small studies are predominately in one direction (usually the 
#' direction of larger effect sizes), asymmetry will ensue.\cr \cr
#' The \code{\link{fat}} function implements several previously proposed tests for detecting funnel plot asymmetry, 
#' which can be used when the presence of between-study heterogeneity is low or neglegible.
#' 
#' @param b Vector with the effect size of each study. Examples are log odds ratio, log hazards ratio, log relative risk. 
#' @param b.se vector with the standard error of the effect size of each study
#' @param method Method for testing funnel plot asymmetry, defaults to \code{"E-UW"} (Egger's test). Other options are \code{E-FIV}. 
#' More info in "Details"
#'
#' @details A common method to test the presence of small-study effects is given as the 
#' following unweighted regression model (\code{method="E-UW"}, Egger 1997): 
#' \deqn{\hat{b}_k = \beta_0 + \beta_1\, \widehat \mathrm{SE}(\hat{b}_k) + \epsilon_k \;,\; \epsilon_k \sim \mathcal{N}\left(0, \sigma^2 \right) }{b = B0 + B1*b.se + e;  e~N(0, s^2)}
#' Whereas \eqn{\beta_0}{B0} indicates the size and direction of the treatment effect, \eqn{\beta_1}{B1} provides 
#' a measure of asymmetry; the larger its deviation from zero the more pronounced the asymmetry. Otherwise, if 
#' \eqn{\beta_1=0}{B1=0}, there is no association between the estimated effect sizes \eqn{\hat{b}_k}{b} and their 
#' corresponding estimates for the standard error \eqn{\widehat \mathrm{SE}(\hat{b}_k)}{b.se} among the reported 
#' studies, indicating no asymmetry and thus no small-study effects. \cr \cr
#' It is possible to allow for between-study heterogeneity by adopting a multiplicative overdispersion parameter 
#' by which the variance in each study is multiplied (\code{method="E-FIV"}, Sterne 2000):
#' \deqn{\hat{\upbeta}_k = a + b\, \widehat \mathrm{SE}(\hat{\upbeta}_k) + \epsilon_k \;,\; \epsilon_k \sim \mathcal{N}(0, \phi \; \widehat \mathrm{var}(\hat{\upbeta}_k))}{b = B0 + B1*b.se + e;  e~N(0, P*b.se^2)}
#' 
#' @return a list containing the following entries:
#' \describe{
##'  \item{"pval"}{A two-sided P-value, based on a Student-t distribution, indicating statistical significance of the funnel plot asymettry test}
##'  \item{"model"}{A fitted \code{glm} object, representing the estimated regression model}
##' }
#' @author Thomas Debray
#' 
#' @references Debray TPA, Moons KGM, Riley RD. Detecting small-study effects and funnel plot asymmetry in meta-analysis of survival data: a comparison of new and existing tests. \emph{Res Syn Meth}. 2017 Oct 3.\cr
#' \cr
#' Egger M, Davey Smith G, Schneider M, Minder C. Bias in meta-analysis detected by a simple, graphical test. \emph{BMJ}. 1997;315(7109):629-34. \cr
#' \cr
#' Sterne JA, Gavaghan D, Egger M. Publication and related bias in meta-analysis: power of statistical tests and prevalence in the literature. \emph{J Clin Epidemiol}. 2000;53(11):1119-29. 


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
#' @export
#' 
#'
fat <- function(b, b.se, method="E-UW") 
{
  
  if (method %in% c("E-UW", "E-FIV")) {
    if (length(b) != length(b.se))
      stop("Incompatible vector sizes for 'b' and 'b.se'!")
    #if (length(b) != length(b.se))
    #  stop("Incompatible vector sizes for 'b' and 'b.se'!")
    
    # Identify studies with complete information
    if (method == "E-UW") {
      studies.complete <- c(!is.na(b) & !is.na(b.se))
      w <- rep(NA, length(b))
    } else if (method== "E-FIV") {
      studies.complete <- c(!is.na(b) & !is.na(b.se))
      w <- 1/(b.se**2)
    }
    
    # Identify which studies can be used
    nstudies <- sum(studies.complete)
    
    # Omit sudies with missing information
    b <- b[studies.complete]
    b.se <- b.se[studies.complete]
    w <- w[studies.complete]
    
    if (nstudies < length(studies.complete)) {
      warning("Some studies were removed due to missing data!")
    }
    
  }
  
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
  class(out) <- "fat"
  return(out)
}