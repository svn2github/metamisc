### To add / change:
# Add plot support for P-FPV, D-FIV and D-FAV
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
#' The \code{\link{fat}} function implements several tests for detecting funnel plot asymmetry, 
#' which can be used when the presence of between-study heterogeneity in treatment effect is relatively low.
#' 
#' @param b Vector with the effect size of each study. Examples are log odds ratio, log hazards ratio, 
#' log relative risk. 
#' @param b.se Optional vector with the standard error of the effect size of each study
#' @param n.total Optional vector with the total sample size of each study
#' @param d.total Optional vector with the total number of observed events for each study
#' @param d1 Optional vector with the total number of observed events in the exposed groups
#' @param d2 Optional vector with the total number of observed events in the unexposed groups
#' @param method Method for testing funnel plot asymmetry, defaults to \code{"E-UW"} (Egger's test). 
#' Other options are \code{E-FIV}, \code{M-FIV}, \code{M-FPV}, \code{D-FIV} and \code{D-FAV}. More info in "Details"
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
#' Because the use of \eqn{\widehat \mathrm{SE}(\hat{b}_k)}{b.se} as independet variable is rather problemetic 
#' when the effect sizes \eqn{\hat{b}_k}{b} represent log odds ratios or log hazard ratios, Macaskill et al. 
#' proposed to use the following regression model (\code{method="M-FIV"}, Macaskill 2001):
#' \deqn{\hat{\upbeta}_k = a + b \,n_k + \epsilon_k \;,\; \epsilon_k \sim \mathcal{N}(0, \phi \; \widehat \mathrm{var}(\hat{\upbeta}_k))}{b = B0 + B1*n.total + e;  e~N(0, P*b.se^2)}
#' Macaskill et al. also proposed an alternative test where a 'pooled' estimate of the outcome proportion is used
#' for the variance \eqn{\widehat \mathrm{var}(\hat{b}_k)}{b.se^2} (\code{method="M-FPV"}, Macaskill 2001):
#' \deqn{\hat{\upbeta}_k = a + b \,n_k + \epsilon_k \;,\; \epsilon_k \sim \mathcal{N}\left(0, \phi \; \frac{1}{d_k (1-d_k/n_k)}\right)}{b = B0 + B1*n.total + e;  e~N(0, P/(d.total * (1-d.total/n.total)))}
#' For studies with zero events, a continuity correction is applied by adding 0.5 to cell counts.
#' A modification of Macaskill's test was proposed by Peters et al. to obtain more balanced type-I error rates 
#' in the tail probability areas  (\code{method="P-FPV"}, Peters 2006):
#' \deqn{\hat{\upbeta}_k = a + b \,\frac{1}{n_k} + \epsilon_k \;,\; \epsilon_k \sim \mathcal{N}\left(0, \phi \; \frac{1}{d_k (1-d_k/n_k)}\right)}{b = B0 + B1/n.total + e;  e~N(0, P/(d.total * (1-d.total/n.total)))}
#' Again, 0.5 is added to all cells for studies with zero events.\cr
#' \cr
#' Because the use of aforementioned tests may be less appropriate in the presence of survival data, Debray et al. 
#' proposed using the total number of events as independent variable (\code{D-FIV}, Debray 2017):
#' \deqn{\hat{\upbeta}_k = a + b\, \frac{1}{d_k} + \epsilon_k  \;,\; \epsilon_k \sim \mathcal{N}(0, \phi \; \widehat \mathrm{var}(\hat{\upbeta}_k))}{b = B0 + B1/d.total + e;  e~N(0, P*b.se^2)}
#' For studies with zero events, the total number of observed events is set to 1.
#' Alternatively, when \eqn{\widehat \mathrm{var}(\hat{\upbeta}_k)}{b.se} is unknown or derived from small samples, 
#' Debray at al.proposed to use the following regression model (\code{D-FAV}, Debray 2017):
#' \deqn{\hat{\upbeta}_k = a + b\, \frac{1}{d_k} + \epsilon_k  \;,\; \epsilon_k \sim \mathcal{N}\left(0, \phi \; \left(\frac{1}{d_{k1}}+\frac{1}{d_{k2}}\right)\right)}{b = B0 + B1/d.total + e;  e~N(0, P/(1/d1 + 1/d2))}
#' 
#' @return a list containing the following entries:
#' \describe{
##'  \item{"pval"}{A two-sided P-value indicating statistical significance of the funnel plot asymettry test. 
##'  Values below the significance level (usually defined as 10\%) support the presence of funnel plot asymmetry,
##'  and thus small-study effects.  }
##'  \item{"model"}{A fitted \code{glm} object, representing the estimated regression model used for testing funnel
##'  plot asymmetry.}
##' }
#' @author Thomas Debray
#' 
#' @references Debray TPA, Moons KGM, Riley RD. Detecting small-study effects and funnel plot asymmetry in 
#' meta-analysis of survival data: a comparison of new and existing tests. \emph{Res Syn Meth}. 2017.\cr
#' \cr
#' Egger M, Davey Smith G, Schneider M, Minder C. Bias in meta-analysis detected by a simple, graphical test. 
#' \emph{BMJ}. 1997;315(7109):629--34. \cr
#' \cr
#' Macaskill P, Walter SD, Irwig L. A comparison of methods to detect publication bias in meta-analysis. 
#' \emph{Stat Med}. 2001;20(4):641--54.\cr 
#' \cr
#' Peters JL, Sutton AJ, Jones DR, Abrams KR, Rushton L. Comparison of two methods to detect publication bias 
#' in meta-analysis. \emph{JAMA}. 2006 Feb 8;295(6):676--80.\cr
#' \cr 
#' Sterne JA, Gavaghan D, Egger M. Publication and related bias in meta-analysis: power of statistical tests 
#' and prevalence in the literature. \emph{J Clin Epidemiol}. 2000;53(11):1119--29. 
#' 
#' @seealso \code{\link{plot.fat}}
#'
#' @examples 
#' data(Fibrinogen)
#' b <- log(Fibrinogen$HR)
#' b.se <- ((log(Fibrinogen$HR.975) - log(Fibrinogen$HR.025))/(2*qnorm(0.975)))
#' n.total <- Fibrinogen$N.total
#' d.total <- Fibrinogen$N.events
#' 
#' fat(b=b, b.se=b.se)
#' fat(b=b, n.total=n.total, d.total=d.total, method="P-FPV")
#' fat(b=b, b.se=b.se, d.total=d.total, method="D-FIV")
#'
#' @import stats
#' @importFrom stats pt
#' @importFrom stats qnorm
#' 
#' @export
fat <- function(b, b.se, n.total, d.total, d1, d2, method="E-UW") 
{
  if (missing(b)) {
    stop ("No values given for 'b'")
  }
  
  # Identify studies with complete information
  if (method == "E-UW") {
    if (missing(b.se)) {
      stop ("No values given for 'b.se'")
    }
    if (length(b) != length(b.se)) {
      stop("Incompatible vector sizes for 'b' and 'b.se'!")
    }
    studies.complete <- c(!is.na(b) & !is.na(b.se))
    ds <- as.data.frame(cbind(b, b.se))
    colnames(ds) <- c("y","x")
  } else if (method== "E-FIV") {
    if (missing(b.se)) {
      stop ("No values given for 'b.se'")
    }
    if (length(b) != length(b.se)) {
      stop("Incompatible vector sizes for 'b' and 'b.se'!")
    }
    studies.complete <- c(!is.na(b) & !is.na(b.se))
    ds <- as.data.frame(cbind(b, b.se, (1/(b.se**2))))
    colnames(ds) <- c("y","x","w")
  } else if (method == "M-FIV") {
    if (missing(b.se)) {
      stop ("No values given for 'b.se'")
    }
    if (missing(n.total)) {
      stop ("No values given for 'n.total'")
    }
    if (length(b) != length(b.se)) {
      stop("Incompatible vector sizes for 'b' and 'b.se'!")
    }
    if (length(b) != length(n.total)) {
      stop("Incompatible vector sizes for 'b' and 'n.total'!")
    }
    studies.complete <- c(!is.na(b) & !is.na(b.se) & !is.na(n.total))
    ds <- as.data.frame(cbind(b, n.total, (1/(b.se**2))))
    colnames(ds) <- c("y","x","w")
  } else if (method=="M-FPV") {
    if (missing(n.total)) {
      stop ("No values given for 'n.total'")
    }
    if (missing(d.total)) {
      stop ("No values given for 'd.total'")
    }
    if (length(b) != length(n.total)) {
      stop("Incompatible vector sizes for 'b' and 'n.total'!")
    }
    if (length(b) != length(d.total)) {
      stop("Incompatible vector sizes for 'b' and 'd.total'!")
    }
    studies.complete <- c(!is.na(b) & !is.na(d.total) & !is.na(n.total))
    
    # Consider continuity corrections
    d.total.cc <- d.total
    d.total.cc[d.total==0] <- 1 #0.5 event in exposed group and 0.5 event in non-exposed group
    n.total[d.total==0] <- n.total[d.total==0]+2 #2*0.5 in the events, and 2*0.5 in the non-events
    
    ds <- as.data.frame(cbind(b, n.total, (d.total.cc*(1-d.total.cc/n.total))))
    colnames(ds) <- c("y","x","w")
  } else if (method=="P-FPV") {
    if (missing(n.total)) {
      stop ("No values given for 'n.total'")
    }
    if (missing(d.total)) {
      stop ("No values given for 'd.total'")
    }
    if (length(b) != length(n.total)) {
      stop("Incompatible vector sizes for 'b' and 'n.total'!")
    }
    if (length(b) != length(d.total)) {
      stop("Incompatible vector sizes for 'b' and 'd.total'!")
    }
    studies.complete <- c(!is.na(b) & !is.na(d.total) & !is.na(n.total))
    
    # Consider continuity corrections
    d.total.cc <- d.total
    d.total.cc[d.total==0] <- 1 #0.5 event in exposed group and 0.5 event in non-exposed group
    n.total[d.total==0] <- n.total[d.total==0]+2 #2*0.5 in the events, and 2*0.5 in the non-events
    
    ds <- as.data.frame(cbind(b, 1/n.total, (d.total.cc*(1-d.total.cc/n.total))))
    colnames(ds) <- c("y","x","w")
  } else if (method=="D-FIV") {
    if (missing(b.se)) {
      stop ("No values given for 'b.se'")
    }
    if (missing(d.total)) {
      stop ("No values given for 'd.total'")
    }
    if (length(b) != length(b.se)) {
      stop("Incompatible vector sizes for 'b' and 'b.se'!")
    }
    if (length(b) != length(d.total)) {
      stop("Incompatible vector sizes for 'b' and 'd.total'!")
    }
    studies.complete <- c(!is.na(b) & !is.na(b.se) & !is.na(d.total))
    
    # Consider continuity corrections
    d.total.cc <- d.total
    d.total.cc[d.total==0] <- 1 #0.5 event in exposed group and 0.5 event in non-exposed group

    ds <- as.data.frame(cbind(b, 1/d.total.cc, (1/(b.se**2))))
    colnames(ds) <- c("y","x","w")
  } else if (method=="D-FAV") {
    if (missing(d1)) {
      stop ("No values given for 'd1'")
    }
    if (missing(d2)) {
      stop ("No values given for 'd2'")
    }
    if (length(b) != length(d1)) {
      stop("Incompatible vector sizes for 'b' and 'd1'!")
    }
    if (length(b) != length(d2)) {
      stop("Incompatible vector sizes for 'b' and 'd2'!")
    }
    if (!missing(d.total)) {
      if (sum(d1+d2!=d.total) > 0)
        stop("Incompatible information between 'd.total', 'd1' and 'd2'")
    }
    studies.complete <- c(!is.na(b) & !is.na(d1) & !is.na(d2))
    
    # Consider continuity corrections
    d1.cc <- d1
    d2.cc <- d2
    d1.cc[(d1==0 | d2==0)] <- d1.cc[(d1==0 | d2==0)]+0.5 #0.5 event in exposed group and 0.5 event in non-exposed group
    d2.cc[(d1==0 | d2==0)] <- d2.cc[(d1==0 | d2==0)]+0.5
    
    ds <- as.data.frame(cbind(b, 1/(d1.cc+d2.cc),  1/((1/d1.cc)+(1/d2.cc))))
    colnames(ds) <- c("y","x","w")
  } 
  else {
    stop("Method for testing funnel plot asymmetry not supported")
  }
  
  # Identify which studies can be used
  nstudies <- sum(studies.complete)
  
  # Omit sudies with missing information
  ds <- ds[studies.complete,]
  
  if (nstudies < length(studies.complete)) {
    warning("Some studies were removed due to missing data!")
  }

  
  if (method %in% c("E-FIV", "M-FIV", "M-FPV", "P-FPV", "D-FIV", "D-FAV")) {
    suppressWarnings(m.fat <- try(glm(y~x, weights=ds$w, data=ds), silent=T))
  } else if (method=="E-UW")  {
    suppressWarnings(m.fat <- try(glm(y~x, data=ds), silent=T))
  } else {
    stop("Method for testing funnel plot asymmetry currently not implemented")
  }
  
  if ("try-error" %in% attr(m.fat,"class")) {
    warning("Estimation of the regression model unsuccessful, P-value omitted.")
    z.fat <- NA
    p.fat <- NA
  } else {
    z.fat <- coefficients(m.fat)[2]/sqrt(diag(vcov(m.fat))[2])
    p.fat <- 2*pt(-abs(z.fat),df=(nstudies-2))
  }

  out <- list()
  out$call <- match.call()
  out$method <- method
  out$pval <- p.fat
  out$nstudies <- nstudies
  out$model <- m.fat
  class(out) <- "fat"
  return(out)
}

#' @export
print.fat <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Call: ");                       
  print(x$call); 
  cat("\n")
  cat("Evidence of asymmetry: \n"); 
  cat(c("\tPr(>|t|) =", round(x$pval, digits = digits))); 
  cat("\n")
}

#' Display results from the funnel plot asymmetry test
#' 
#' Generates a funnel plot for a fitted \code{fat} object, with boundaries for the 90\% confidence interval (based on a Student-T distribution)
#' @param x An object of class \code{fat}
#' @param ... Additional arguments for \code{\link{plot}}. The argument \code{funnel.xlab} can be used to define the x-label axis.
#' The argument \code{ref} represents a numeric value indicating the fixed or random effects summary estimate. If no value is provided
#' then it will be retrieved from a fixed effects meta-analysis (if possible). The argument \code{x.rescale} is a character string indicating how effect sizes should be rescaled. Options are \code{exp} 
#' (e.g. for log odds or log hazard ratios).
#' 
#' @examples 
#' data(Fibrinogen)
#' b <- log(Fibrinogen$HR)
#' b.se <- ((log(Fibrinogen$HR.975) - log(Fibrinogen$HR.025))/(2*qnorm(0.975)))
#' n.total <- Fibrinogen$N.total
#' 
#' plot(fat(b=b, b.se=b.se))
#' plot(fat(b=b, b.se=b.se, n.total=n.total, method="M-FIV"), funnel.xlab="Log hazard ratio")
#' plot(fat(b=b, b.se=b.se), funnel.xlab="Hazard ratio", x.rescale="exp")
#' @importFrom metafor rma
#' @importFrom stats qt
#' @importFrom methods hasArg
#' @export
plot.fat <- function(x,  ...) {
  funnel.xlab <- ifelse(hasArg(funnel.xlab), list(...)$funnel.xlab, "Effect size")
  ref <- ifelse(hasArg(ref), list(...)$ref, NA)
  x.rescale <- ifelse(hasArg(x.rescale), list(...)$x.rescale, NA)
  
  plot_fat(object=x, ref=ref, funnel.xlab=funnel.xlab, x.rescale=x.rescale)

}

plot_fat <- function (object, ref, funnel.xlab, x.rescale) {
  xval <-  object$model$data[,"y"]
  
  if (object$method %in% c("E-UW", "E-FIV")) {
    ylab <- "Standard error"
    yval <- (object$model$data[,"x"])
    ylim <- rev(c(0, max(yval, na.rm=T))) #Reverse y axis scale
    yval.min <- -1
    
    # Get the fixed effect summary estimate
    res <- rma(yi=object$model$data[,"y"], sei=object$model$data[,"x"], method="FE")
  } else if (object$method %in% c("M-FIV")) {
    ylab <- "Sample size"
    yval <- (object$model$data[,"x"]) # Sample size
    ylim <- (c(0, max(yval, na.rm=T))) #Reverse y axis scale
    yval.min <- -max(yval) #Generate a minimum value for predictions
    
    # Get the fixed effect summary estimate
    res <- rma(yi=object$model$data[,"y"], sei=(1/sqrt(object$model$data[,"w"])), method="FE")
  } else {
    stop("Plot not supported!")
  }
  
  
  if(is.na(x.rescale)) {
    plot(NULL, xlim=c(min(c(0,xval)), max(xval)), ylim=ylim, 
         ylab=ylab, xlab=funnel.xlab)
  } else if (x.rescale=="exp") {
    plot(NULL, xlim=c(min(c(0,xval)), max((xval))), ylim=ylim, 
         ylab=ylab, xlab=funnel.xlab, xaxt="n")
    axis(1, at=c(log(0.5), log(1), log(2), log(3), log(4)), labels=c(0.5, 1,2,3,4))
  } else {
    stop("Provided argument for 'x.rescale' not supported!")
  }
  
  
  
  
  newdata <- sort(c(yval.min,object$model$data[,"x"], 2*max(object$model$data[,"x"])))
  newdata <- as.data.frame(cbind(newdata,NA))
  colnames(newdata) <- c("x","y")
  predy <- predict(object$model, newdata=newdata, se.fit=T)#
  predy.lowerInt <- as.vector(predy$fit + qt(0.05, df=object$nstudies-2)*predy$se.fit) #90% confidence band
  predy.upperInt <- as.vector(predy$fit + qt(0.95, df=object$nstudies-2)*predy$se.fit) #90% confidence band
  
  polygon(x=c(predy.upperInt,rev(predy.lowerInt)), y=c(newdata[,"x"],rev(newdata[,"x"])), col="skyblue")  
  lines(x=as.vector(predy$fit), y=(newdata[,"x"]), lty=2 )
  points(xval, yval, pch=19)
  box()
  if (is.na(ref)) {
    abline(v=(res$b))
  } else {
    abline(v=ref)
  }
}