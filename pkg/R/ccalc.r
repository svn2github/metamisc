#' Calculate the concordance statistic
#'
#' The function calculates the (transformed) concordance statistic with the corresponding sampling variance
#' from aggregate data. 
#' 
#' @param cstat Optional vector with the estimated c-statistic for each valiation
#' @param cstat.se Optional vector with the standard error of the estimated c-statistics
#' @param cstat.95CI Optional 2-dimensional array with the lower (first column) and upper (second column) boundary 
#' of the 95\% confidence interval of the estimated c-statistics
#' @param sd.LP Optional vector with the standard deviation of the linear predictor (prognostic index) for each validation
#' @param N Optional vector with the total number of participants for each valiation
#' @param O Optional vector with the total number of observed events for each valiation
#' (if specified, during time \code{t.val})
#' @param E Optional vector with the total number of expected events for each valiation 
#' (if specified, during time \code{t.val})
#' @param Po Optional vector with the (cumulative) observed event probability for each valiation
#' @param Pe Optional vector with the (cumulative) expected event probability for each validation
#' @param slab Optional vector with study names
#' @param pars A list with additional arguments: significance level of the confidence intervals (default: \code{level=0.95}), 
#' method for restoring missing estimates for the standard error of the c-statistic 
#' (default: \code{method.restore.c.se="Newcombe.4"}), and likelihood/link for modeling the c-statistic
#' (default: \code{model = "normal/logit"}). see "Details" for more information.
#' @param \ldots Additional arguments.
#' 
#' @details 
#' By default, the logit transformation is applied to the concordance (c-) statistic. If no transformation is needed, 
#' simply set \code{pars$model} equal to \code{"normal/identity"}.
#' 
#' \subsection{Restoring the c-statistic}{
#' For studies where the c-statistic is missing, it is estimated from the standard deviation of the linear predictor 
#' (\code{theta.source="std.dev(LP)"). The corresponding method is described by White et al. (2015). }.
#' }
#' 
#' \subsection{Restoring the standard error of the c-statistic}{
#' When missing, the standard error of the c-statistic can be estimated from the confidence interval. Alternatively, 
#' the standard error can be approximated from a 
#' combination of the reported c-statistic, the total sample size and the total number of events (Debray et al. 2017).
#' The latter strategy can be specified with \code{pars$method.restore.c.se}. In particular, it is possible to adopt
#' the method of Hanley and McNeil (\code{"Newcombe.2"}), also described as method 2 by Newcombe (2006). 
#' Alternatively, it is possible to adopt method 4 (\code{"Newcombe.4"}), a small modification of method 2. 
#' As discussed by Newcombe, both methods may produce zero width intervals in extreme cases.
#' }
#' 
#' @references 
#' \itemize{
#' \item Debray TPA, Damen JAAG, Snell KIE, Ensor J, Hooft L, Reitsma JB, et al. A guide to systematic review 
#' and meta-analysis of prediction model performance. \emph{BMJ}. 2017; 356:i6460.
#' \item Hanley JA, McNeil BJ. The meaning and use of the area under a receiver operating characteristic (ROC) 
#' curve. \emph{Radiology}. 1982; 143(1):29--36.
#' \item Newcombe RG. Confidence intervals for an effect size measure based on the Mann-Whitney statistic. 
#' Part 2: asymptotic methods and evaluation. \emph{Stat Med}. 2006; 25(4):559--73.
#' \item White IR, Rapsomaniki E, the Emerging Risk Factors Collaboration. Covariate-adjusted measures of discrimination 
#' for survival data. \emph{Biom J}. 2015;57(4):592--613. 
#' }
#' 
#' @return An array with the following columns:
#' \describe{
##'  \item{"theta"}{The (transformed) c-statistics. }
##'  \item{"theta.se"}{Standard errors of the (transformed) c-statistics.}
##'  \item{"theta.CIl"}{Lower confidence interval of the (transformed) c-statistics. The level is specified in
##'  \code{pars$level}.}
##'  \item{"theta.CIu"}{Upper confidence interval of the (transformed) c-statistics. The level is specified in
##'  \code{pars$level}.}
##'  \item{"theta.source"}{Method used for calculating the (transformed) c-statistic.}
##'  \item{"theta.se.source"}{Method used for calculating the standard error of the (transformed) c-statistic.}
##' }
#' 
#' @examples 
#' ######### Validation of prediction models with a binary outcome #########
#' data(EuroSCORE)
#' 
#' # Calculate the logit c-statistic and its standard error
#' with(EuroSCORE, ccalc(cstat=c.index, cstat.se=se.c.index, 
#'                       cstat.95CI=cbind(c.index.95CIl,c.index.95CIu), 
#'                       N=n, O=n.events, slab=Study))
#'   
#' # Calculate the c-statistic and its standard error
#' with(EuroSCORE, ccalc(cstat=c.index, cstat.se=se.c.index, 
#'                       cstat.95CI=cbind(c.index.95CIl,c.index.95CIu), 
#'                      N=n, O=n.events, slab=Study, pars=list(model="normal/identity")))
#'                                                             
#' @keywords meta-analysis discrimination  calibration extraction
#' 
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @export
#' 
ccalc <- function(cstat, cstat.se, cstat.95CI, sd.LP, N, O, E, Po, Pe, slab, pars, ...) {
  pars.default <- list(level = 0.95,
                       method.restore.c.se="Newcombe.4",
                       model = "normal/logit") #Alternative: "normal/identity"
  
  if (!missing(pars)) {
    for (i in 1:length(pars)) {
      element <- ls(pars)[i]
      pars.default[[element]] <- pars[[element]]
    }
  }

  
  #######################################################################################
  # Count number of studies
  #######################################################################################
  k <- 0
  
  if (!missing(cstat)) {
    k <- length(cstat)
  } else if (!missing(cstat.se)) {
    k <- length(cstat.se)
  } else if (!missing(cstat.95CI)) {
    k <- dim(cstat.95CI)[2]
  } else if (!missing(sd.LP)) {
    k <- length(sd.LP)
  }
  
  #######################################################################################
  # Prepare data
  #######################################################################################
  if (missing(O)) {
    O <- rep(NA, length=k)
  }
  if (missing(Po)) {
    Po <- rep(NA, length=k)
  }
  if (missing(N)) {
    N <- rep(NA, length=k)
  }
  
  if (missing(cstat.95CI)) {
    cstat.95CI <- array(NA, dim=c(k,2))
  }
  if (is.null(dim(cstat.95CI))) {
    warning("Invalid dimension for 'cstat.95CI', argument ignored.")
    cstat.95CI <- array(NA, dim=c(k,2))
  }
  if (dim(cstat.95CI)[2] != 2 | dim(cstat.95CI)[1] != k) {
    warning("Invalid dimension for 'cstat.95CI', argument ignored.")
    cstat.95CI <- array(NA, dim=c(k,2))
  }
  if (missing(cstat.se)) {
    cstat.se <- array(NA, dim=k)
  }
  if (missing(sd.LP)) {
    sd.LP <- rep(NA, k)
  }
    
  # Calculate O and N from other information if possible
  O <- ifelse(is.na(O), Po*N, O)
  N <- ifelse(is.na(N), O/Po, N)
  
  theta.cil <- theta.ciu <- rep(NA, k)
  
  # Restore c-statistic
  te.method <- c("c-statistic", "std.dev(LP)")
  te.orig  <- restore.c.c(cstat, model=pars.default$model)
  te.white <- restore.c.sdPI(sd.LP, model=pars.default$model)
  
  te.dat <- cbind(te.orig, te.white)
    
  # For each study, find the first colum without missing
  myfun = function(dat) { which.min(is.na(dat)) }
  
  sel.cstat <- apply(te.dat, 1, myfun)
  theta <- te.dat[cbind(seq_along(sel.cstat), sel.cstat)]                            
  theta.source <-  te.method[sel.cstat]
  
  # Define theta
  if (pars.default$model == "normal/identity") {
    if (pars.default$level==0.95) {
      theta.cil <- cstat.95CI[,1]
      theta.ciu <- cstat.95CI[,2]
    }
  } else if (pars.default$model == "normal/logit") {
    if (pars.default$level==0.95) {
      theta.cil <- logit(cstat.95CI[,1])
      theta.ciu <- logit(cstat.95CI[,2])
    }
  } else {
    stop("Undefined link function!")
  }
    
    
    # Calculate all the possible variations of var(theta)
    tv.method <- c("Standard Error", "Confidence Interval", pars.default$method.restore.c.se, pars.default$method.restore.c.se)
    tv.se     <- restore.c.var.se(c.se=cstat.se, cstat=cstat, model=pars.default$model) # Derived from standard error
    tv.ci     <- restore.c.var.ci(ci=cstat.95CI, level=0.95, model=pars.default$model) # Derived from 95% confidence interval
    tv.hanley <- restore.c.var.hanley(cstat=cstat, N.subjects=N, N.events=O, restore.method=pars.default$method.restore.c.se,
                                      model=pars.default$model)
    tv.hanley2 <- restore.c.var.hanley2(sd.LP=sd.LP, N.subjects=N, N.events=O, restore.method=pars.default$method.restore.c.se,
                                        model=pars.default$model)
    
    # Save all estimated variances. The order of the columns indicates the priority             
    dat <-cbind(tv.se, tv.ci, tv.hanley, tv.hanley2)  
    
    
    
    sel.var <- apply(dat, 1, myfun)
    theta.var <- dat[cbind(seq_along(sel.var), sel.var)]                            
    theta.var.source <-  tv.method[sel.var]
    
    # Calculate the desired confidence intervals
    theta.cil[is.na(theta.cil)] <- (theta+qnorm((1-pars.default$level)/2)*sqrt(theta.var))[is.na(theta.cil)]
    theta.ciu[is.na(theta.ciu)] <- (theta+qnorm((1+pars.default$level)/2)*sqrt(theta.var))[is.na(theta.ciu)]
    
    
    # Store results, and method for calculating SE
    ds <- data.frame(theta=theta, theta.se=sqrt(theta.var), theta.CIl=theta.cil, theta.CIu=theta.ciu, 
                     theta.source=theta.source, theta.se.source=theta.var.source)
    
    
    # Assing study labels as rownames
    if(!missing(slab)) {
      slab <- make.unique(as.character(slab))
      rownames(ds) <- slab
    }
    
    return(ds)
}