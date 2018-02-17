#' Calculate the concordance statistic
#'
#' The function calculates the (logit transformed) concordance (c-) statistic with the corresponding sampling variance. 
#' 
#' @param cstat vector to specify the estimated c-statistics.
#' @param cstat.se vector to specify the corresponding standard errors
#' @param cstat.cilb vector to specify the lower limit of the 95\% confidence interval.
#' @param cstat.ciub vector to specify the upper limit of the 95\% confidence interval.
#' @param sd.LP vector to specify the standard deviations of the linear predictor (prognostic index)
#' @param N vector to specify the validation study sizes.
#' @param O vector to specify the total number of observed events.
#' @param Po vector to specify the observed event probabilities.
#' @param data optional data frame containing the variables given to the arguments above.
#' @param slab optional vector with labels for the studies.
#' @param pars optional list to specify additional arguments: 
#' significance level of the confidence interval (default: \code{level=0.95}), 
#' method for calculating the standard error of the c-statistic 
#' (default: \code{method.restore.c.se="Newcombe.4"}). It is possible to apply the logit transformation to 
#' the c-statistics by setting \code{model = "normal/logit"}.
#' @param \ldots Additional arguments.
#' 
#' @details 
#' The c-statistic is a measure of discrimination, and indicates the ability of a prediction model to 
#' distinguish between patients developing and not developing the outcome. The c-statistic typically ranges 
#' from 0.5 (no discriminative ability) to 1 (perfect discriminative ability). 
#' 
#' When performing a meta-analysis of the c-statistic, it is generally recommended to apply the logit transformation 
#' (Debray et al., 2017; Snell et al., 2017). This can be achieved by specifying \code{model} in the additional arguments. 
#' An example is given below.
#' 
#' \subsection{Restoring the c-statistic}{
#' For studies where the c-statistic is missing, it is estimated from the standard deviation of the linear predictor 
#' (\code{theta.source="std.dev(LP)"). The corresponding method is described by White et al. (2015). }.
#' }
#' 
#' \subsection{Restoring the standard error of the c-statistic}{
#' When missing, the standard error of the c-statistic can be estimated from the confidence interval. Alternatively, 
#' the standard error can be approximated from a combination of the reported c-statistic, the total sample size and 
#' the total number of events (Debray et al. 2017). This can be achieved by adopting (a modification of) the method 
#' proposed by Hanley and McNeil, as specified in \code{pars$method.restore.c.se}. Possible options for this argument are
#' \code{"Newcombe.2"} and \code{"Newcombe.4"}, which correspond to method 2 and, respectively, method 4 in Newcombe (2006).
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
#' \item Snell KI, Ensor J, Debray TP, Moons KG, Riley RD. Meta-analysis of prediction model performance across 
#' multiple studies: Which scale helps ensure between-study normality for the C -statistic and calibration measures? 
#' \emph{Statistical Methods in Medical Research}. 2017. 
#' \item White IR, Rapsomaniki E, the Emerging Risk Factors Collaboration. Covariate-adjusted measures of discrimination 
#' for survival data. \emph{Biom J}. 2015;57(4):592--613. 
#' }
#' 
#' 
#' @return An array with the following columns:
#' \describe{
##'  \item{"theta"}{The (transformed) c-statistics. }
##'  \item{"theta.se"}{Standard errors of the (transformed) c-statistics.}
##'  \item{"theta.CIl"}{Lower confidence interval of the (transformed) c-statistics. The level is specified in
##'  \code{pars$level}. Intervals are calculated on the same scale as \code{theta} by assuming a Normal distribution.}
##'  \item{"theta.CIu"}{Upper confidence interval of the (transformed) c-statistics. The level is specified in
##'  \code{pars$level}. Intervals are calculated on the same scale as \code{theta} by assuming a Normal distribution.}
##'  \item{"theta.source"}{Method used for calculating the (transformed) c-statistic.}
##'  \item{"theta.se.source"}{Method used for calculating the standard error of the (transformed) c-statistic.}
##' }
#' 
#' @examples 
#' ######### Validation of prediction models with a binary outcome #########
#' data(EuroSCORE)
#' 
#' # Calculate the c-statistic and its standard error
#' ccalc(cstat=c.index, cstat.se=se.c.index, cstat.cilb=c.index.95CIl, cstat.ciub=c.index.95CIu, 
#'       N=n, O=n.events, data=EuroSCORE, slab=Study)
#'   
#' # Calculate the logit c-statistic and its standard error
#' ccalc(cstat=c.index, cstat.se=se.c.index, cstat.cilb=c.index.95CIl, cstat.ciub=c.index.95CIu, 
#'       N=n, O=n.events, data=EuroSCORE, slab=Study, pars=list(model="normal/logit"))
#'                                                             
#' @keywords meta-analysis discrimination concordance statistic performance
#' 
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @export
#' 
ccalc <- function(cstat, cstat.se, cstat.cilb, cstat.ciub, sd.LP, N, O, Po, data, slab, pars, ...) {
  pars.default <- list(level = 0.95,
                       method.restore.c.se="Newcombe.4",
                       model = "normal/identity") #Alternative: "normal/logit"
  
  #######################################################################################
  # Set default parameters
  #######################################################################################
  if (!missing(pars)) {
    for (i in 1:length(pars)) {
      element <- ls(pars)[i]
      pars.default[[element]] <- pars[[element]]
    }
  }
  
  ### check if data argument has been specified
  if (missing(data))
    data <- NULL
  
  ### need this at the end to check if append=TRUE can actually be done
  no.data <- is.null(data)
  
  ### check if data argument has been specified
  if (is.null(data)) {
    data <- sys.frame(sys.parent())
  } else {
    if (!is.data.frame(data))
      data <- data.frame(data)
  }
  
  #######################################################################################
  # Retrieve all data
  #######################################################################################
  mf <- match.call()
  
  mf.slab       <- mf[[match("slab",   names(mf))]]
  slab          <- eval(mf.slab,   data, enclos=sys.frame(sys.parent()))
  mf.cstat      <- mf[[match("cstat", names(mf))]]
  cstat         <- eval(mf.cstat, data, enclos=sys.frame(sys.parent()))
  mf.cstat.se   <- mf[[match("cstat.se", names(mf))]]
  cstat.se      <- eval(mf.cstat.se, data, enclos=sys.frame(sys.parent()))
  mf.cstat.cilb <- mf[[match("cstat.cilb", names(mf))]]
  cstat.cilb    <- eval(mf.cstat.cilb, data, enclos=sys.frame(sys.parent()))
  mf.cstat.ciub <- mf[[match("cstat.ciub", names(mf))]]
  cstat.ciub    <- eval(mf.cstat.ciub, data, enclos=sys.frame(sys.parent()))
  mf.sd.LP      <- mf[[match("sd.LP", names(mf))]]
  sd.LP         <- eval(mf.sd.LP, data, enclos=sys.frame(sys.parent()))
  mf.N          <- mf[[match("N", names(mf))]]
  N             <- eval(mf.N, data, enclos=sys.frame(sys.parent()))
  mf.O          <- mf[[match("O", names(mf))]]
  O             <- eval(mf.O, data, enclos=sys.frame(sys.parent()))
  mf.Po         <- mf[[match("Po", names(mf))]]
  Po           <- eval(mf.Po, data, enclos=sys.frame(sys.parent()))
  
  #######################################################################################
  # Count number of studies
  #######################################################################################
  k <- 0
  
  if (!no.data) {
    k <- dim(data)[1]
  } else if (!is.null(cstat)) {
    k <- length(cstat)
  } else if (!is.null(cstat.se)) {
    k <- length(cstat.se)
  } else if (!is.null(cstat.cilb)) {
    k <- length(cstat.cilb)
  } else if (!is.null(cstat.ciub)) {
    k <- length(cstat.ciub)
  } else if (!is.null(sd.LP)) {
    k <- length(sd.LP)
  }
  
  #######################################################################################
  # Prepare data
  #######################################################################################
  if (is.null(O)) {
    O <- rep(NA, length=k)
  }
  if (is.null(Po)) {
    Po <- rep(NA, length=k)
  }
  if (is.null(N)) {
    N <- rep(NA, length=k)
  }
  if (is.null(cstat.cilb)) {
    cstat.cilb <- rep(NA, length=k)
  }
  if (is.null(cstat.ciub)) {
    cstat.ciub <- rep(NA, length=k)
  }
  if (is.null(cstat.se)) {
    cstat.se <- array(NA, dim=k)
  }
  if (is.null(sd.LP)) {
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
      theta.cil <- cstat.cilb
      theta.ciu <- cstat.ciub
    }
  } else if (pars.default$model == "normal/logit") {
    if (pars.default$level==0.95) {
      theta.cil <- logit(cstat.cilb)
      theta.ciu <- logit(cstat.ciub)
    }
  } else {
    stop("Supplied model for transforming the c-statistic is not supported!")
  }
  
  cstat.95CI <- cbind(cstat.cilb, cstat.ciub)
  
    
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