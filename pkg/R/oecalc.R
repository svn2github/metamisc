#' Calculate the total O:E ratio
#'
#' This function calculates (transformed versions of) the ratio of total number of observed versus expected events with the 
#' corresponding sampling variance. 
#' 
#' @param OE vector with the estimated ratio of total observed versus total expected events
#' @param OE.se vector with the standard errors of the estimated O:E ratios
#' @param OE.cilb vector to specify the lower limits of the confidence interval for \code{OE}.
#' @param OE.ciub vector to specify the upper limits of the confidence interval for \code{OE}.
#' @param OE.cilv vector to specify the levels of aformentioned confidence interval limits. 
#' (default: 0.95, which corresponds to the 95\% confidence interval).
#' @param citl vector with the estimated calibration-in-the-large statistics
#' @param citl.se vector with the standard error of the calibration-in-the-large statistics
#' @param N vector to specify the validation study sizes.
#' @param O vector to specify the total number of observed events.
#' @param E vector to specify the total number of expected events
#' @param Po vector to specify the (cumulative) observed event probabilities.
#' @param Po.se vector with the standard errors of \code{Po}.
#' @param Pe vector to specify the (cumulative) expected event probabilites
#' (if specified, during time \code{t.val})
#' @param data optional data frame containing the variables given to the arguments above.
#' @param slab optional vector with labels for the studies.
#' @param add a non-negative number indicating the amount to add to zero counts. See `Details'
#' @param g a quoted string that is the function to transform estimates of the total O:E ratio; see the details below.
#' @param level level for confidence interval, default \code{0.95}.
#' @param \ldots Additional arguments.
#' 
#' @examples 
#' ######### Validation of prediction models with a binary outcome #########
#' data(EuroSCORE)
#' 
#' # Calculate the log of the total O:E ratio and its standard error
#' oecalc(O=n.events, E=e.events, N=n, data=EuroSCORE, slab=Study, g="log(OE)")
#' 
#' @keywords meta-analysis calibration performance
#' 
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @export
oecalc <- function(OE, OE.se, OE.cilb, OE.ciub, OE.cilv, citl, citl.se, N, O, E, Po, Po.se, Pe, 
                   data, slab, add=1/2, g=NULL, level=0.95, ...) {
  
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
  mf.OE         <- mf[[match("OE", names(mf))]]
  OE            <- eval(mf.OE, data, enclos=sys.frame(sys.parent()))
  mf.OE.se      <- mf[[match("OE.se", names(mf))]]
  OE.se         <- eval(mf.OE.se, data, enclos=sys.frame(sys.parent()))
  mf.OE.cilb    <- mf[[match("OE.cilb", names(mf))]]
  OE.cilb       <- eval(mf.OE.cilb, data, enclos=sys.frame(sys.parent()))
  mf.OE.ciub    <- mf[[match("OE.ciub", names(mf))]]
  OE.ciub       <- eval(mf.OE.ciub, data, enclos=sys.frame(sys.parent()))
  mf.OE.cilv    <- mf[[match("OE.cilv", names(mf))]]
  OE.cilv       <- eval(mf.OE.cilv, data, enclos=sys.frame(sys.parent()))
  mf.citl       <- mf[[match("citl", names(mf))]]
  citl          <- eval(mf.citl, data, enclos=sys.frame(sys.parent()))
  mf.citl.se    <- mf[[match("citl.se", names(mf))]]
  citl.se       <- eval(mf.citl.se, data, enclos=sys.frame(sys.parent()))
  mf.N          <- mf[[match("N", names(mf))]]
  N             <- eval(mf.N, data, enclos=sys.frame(sys.parent()))
  mf.O          <- mf[[match("O", names(mf))]]
  O             <- eval(mf.O, data, enclos=sys.frame(sys.parent()))
  mf.E          <- mf[[match("E", names(mf))]]
  E             <- eval(mf.E, data, enclos=sys.frame(sys.parent()))
  mf.Po         <- mf[[match("Po", names(mf))]]
  Po            <- eval(mf.Po, data, enclos=sys.frame(sys.parent()))
  mf.Pe         <- mf[[match("Pe", names(mf))]]
  Pe            <- eval(mf.Pe, data, enclos=sys.frame(sys.parent()))
  mf.Po.se      <- mf[[match("Po.se", names(mf))]]
  Po.se         <- eval(mf.Po.se, data, enclos=sys.frame(sys.parent()))
  
  #######################################################################################
  # Count number of studies
  #######################################################################################
  k <- 0
  
  if (!no.data) {
    k <- dim(data)[1]
  } else if (!is.null(OE)) {
    k <- length(OE)
  } else if (!is.null(OE.se)) {
    k <- length(OE.se)
  } else if (!is.null(OE.cilb)) {
    k <- length(OE.cilb)
  } else if (!is.null(OE.ciub)) {
    k <- length(OE.ciub)
  } else if (!is.null(OE.cilv)) {
    k <- length(OE.cilv)
  } else if (!is.null(citl)) {
    k <- length(citl)
  } else if (!is.null(citl.se)) {
    k <- length(citl.se)
  } else if (!is.null(N)) {
    k <- length(N)
  }  else if (!is.null(O)) {
    k <- length(O)
  } else if (!is.null(E)) {
    k <- length(E)
  } else if (!is.null(Po)) {
    k <- length(Po)
  } else if (!is.null(Po.se)) {
    k <- length(Po.se)
  } else if (!is.null(Pe)) {
    k <- length(Pe)
  }

  if (k<1) stop("No data provided!")
  
  if(is.null(OE)) {
    OE <- rep(NA, times=k)
  }
  
  #######################################################################################
  # Assign study labels
  # taken from escalc
  #######################################################################################
  if (!is.null(slab)) {
    
    if (anyNA(slab))
      stop("NAs in study labels.")
    
    if (class(slab)=="factor") {
      slab <- as.character(slab)
    }
    
    ### check if study labels are unique; if not, make them unique
    
    if (anyDuplicated(slab))
      slab <- make.unique(slab)
    
    if (length(slab) != k)
      stop("Study labels not of same length as data.")
  }
  

  
  
  #######################################################################################
  # Restore OE ratio
  #######################################################################################
  t.O.E.N   <- restore.oe.O.E.N(O=O, E=E, N=N, correction = add, g=g) 
  
  
  
  #t.O.Pe.N  <- restore.oe.OPeN(O=O, Pe=Pe, N=N, correction = pars.default$correction, 
  #                             t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, model=pars.default$model.oe)
  #t.E.Po.N  <- restore.oe.EPoN(E=E, Po=Po, N=N, correction = pars.default$correction, 
  #                             t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, model=pars.default$model.oe)
  #T.Po.Pe.N <- restore.oe.PoPeN(Po=Po, Pe=Pe, N=N, correction = pars.default$correction, 
  #                              t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, model=pars.default$model.oe) 
  #t.O.Po.E  <- restore.oe.OPoE(O=O, Po=Po, E=E, correction = pars.default$correction, 
  #                             t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, model=pars.default$model.oe) 
  #t.O.Pe.E  <- restore.oe.OPeE(O=O, Pe=Pe, E=E, correction = pars.default$correction, 
  #                             t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, model=pars.default$model.oe) 
  #t.OE.SE   <- restore.oe.OE(OE=OE, OE.se=OE.se, t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, 
  #                           model=pars.default$model.oe)
  #t.OE.CI   <- restore.oe.OE.95CI(OE=OE, OE.95CI=OE.95CI, t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, 
  #                                model=pars.default$model.oe)
  #t.O.E     <- restore.oe.O.E(O=O, E=E, correction = pars.default$correction, 
  #                            t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, model=pars.default$model.oe) 
  #t.pope   <- restore.oe.PoPe(Po=Po, Pe=Pe, t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, model=pars.default$model.oe) 
  #t.citl   <- restore.oe.citl(citl=citl, citl.se=citl.se, O=O, Po=Po, N=N, t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, 
  #                            model=pars.default$model.oe) 
}

# Calculate OE and its error variance from O, E and N
restore.oe.O.E.N <- function(O, E, N, correction = 0.5, g=NULL) {
  
  k <- length(O)
  out <- array(NA, dim=c(k,2))
  cc <- which(E==0)
  E[cc] <- E[cc]+correction
  O[cc] <- O[cc]+correction
  N[cc] <- N[cc]+correction
  out[,1] <- O/E 
  out[,2] <- sqrt((O*(1-O/N))/(E**2))
  
  if(is.null(g)) {
    return (out)
  }

  logoe <- logoe.var <- rep(NA, k)
  
  for (i in 1:k) {
    oei <- out[i,1]
    logoe[i] <- eval(parse(text=g), list(OE = oei))
    vi  <- out[i,2]**2
    names(oei) <- names(vi) <- "OE"
    logoe.var[i] <- as.numeric((deltaMethod(object=oei, g=g, vcov.=vi))["SE"])**2
  }
  
  out <- cbind(logoe, logoe.var)
  
  return (out)
}
  