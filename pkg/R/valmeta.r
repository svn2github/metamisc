valmeta <- function(x, ...) {
    UseMethod("valmeta")
  #For OE ratio lets use valmeta.data.frame
}
  

# By default try meta-analysis of c-statistic
valmeta.default <- function(cstat, cstat.se, cstat.95CI,
                            N, O, method="REML", knha=TRUE, verbose=FALSE, 
                            method.restore.se="Newcombe.4", scale = "logit", ...) {

  ds <- cbind(cstat, cstat.se)
  
  out <- list()
  out$call <- match.call()
  out$method <- method
  out$method.restore.se <- method.restore.se 
  out$scale <- scale
  class(out) <- "valmeta"
  
  if (missing(cstat.se) & missing(cstat.95CI)) {
    stop("No sampling error was provided for the c-statistic!")
  }
  if (missing(cstat.95CI)) {
    cstat.95CI <- array(NA, dim=c(length(cstat),2))
  }
  if (missing(cstat.se)) {
    cstat.se <- array(NA, dim=length(cstat))
  }
  if (dim(cstat.95CI)[2] != 2 | dim(cstat.95CI)[1] != length(cstat)) {
    stop("Invalid dimension for variable 'cstat.95CI'!")
  }

  inv.logit <- function(x) {  if(is.numeric(x)) 1/(1+exp(-x)) else stop("x is not numeric!") }
  logit <- function(x) { log(x/(1-x)) }
  
  #Update SE(c.index) using Method 4 of Newcombe
  restore.c.var<- function(cstat, N.subjects, N.events, restore.method="Newcombe.4", scale=scale) {
    n <- N.events #Number of events
    m <- N.subjects-N.events #Number of non-events
    
    if (missing(restore.method)) {
      restore.method <- "Newcombe.4"
    }
    
    if (restore.method=="Hanley" | restore.method=="Newcombe.2") {
      mstar <- m-1
      nstar <- n-1
    } else if (restore.method=="Newcombe.4") {
      mstar <- nstar <- N.subjects/2-1
    } else {
      stop ("Method not implemented yet!")
    }
    
    if (scale=="logit") {
      out <- (((1+nstar*(1-cstat)/(2-cstat) + mstar*cstat/(1+cstat)))/(m*n*cstat*(1-cstat)))
    } else {
      out <- ((cstat*(1-cstat)*(1+nstar*(1-cstat)/(2-cstat) + mstar*cstat/(1+cstat)))/(m*n))
    }
    
    return(out)
  }
  
  if (method != "BAYES") { # Use of rma
    num.estimated.var.c <- 0
    
    # Apply necessary data transformations
    if (scale == "identity") {
      theta <- cstat
      theta.var <- (cstat.se)**2
      theta.var.CI <- ((cstat.95CI[,2] - cstat.95CI[,1])/(2*qnorm(0.975)))**2 #Derive from 95% CI
      theta.var <- ifelse(is.na(theta.var), theta.var.CI, theta.var)
    } else if (scale == "logit") {
      theta <- log(cstat/(1-cstat))
      theta.var <- (cstat.se/(cstat*(1-cstat)))**2
      theta.var.CI <- ((logit(cstat.95CI[,2]) - logit(cstat.95CI[,1]))/(2*qnorm(0.975)))**2
      theta.var <- ifelse(is.na(theta.var), theta.var.CI, theta.var)
    } else {
      stop("No appropriate transformation defined")
    }
    
    # Restore missing standard errors
    if (NA %in% theta.var) {
      # Restore missing estimates of the standard error of the c-statistic using information on c, N and O
      if (!missing(O) & !missing(N)) {
        if (verbose) cat("Attempting to restore missing information on the standard error of the c-statistic\n")
        theta.var.hat <- restore.c.var(cstat=cstat, N.subjects=N, N.events=O, restore.method=method.restore.se, scale=scale)
        num.estimated.var.c <- length(which(is.na(theta.var) & !is.na(theta.var.hat)))
        theta.var <- ifelse(is.na(theta.var), theta.var.hat, theta.var)
      }
    }
    ds <- cbind(ds, theta, theta.var)
    
    out$data <- ds
    out$num.estimated.var.c <- num.estimated.var.c
    
    # Apply the meta-analysis
    fit <- rma(yi=theta, vi=theta.var, data=ds, method=method, knha=knha) 
    preds <- predict(fit)
    
    results <- as.data.frame(array(NA, dim=c(1,5)))
    if (scale == "logit") {
      results <- c(inv.logit(coefficients(fit)), inv.logit(c(preds$ci.lb, preds$ci.ub, preds$cr.lb, preds$cr.ub)))
    } else {
      results <- c(coefficients(fit), c(preds$ci.lb, preds$ci.ub, preds$cr.lb, preds$cr.ub))
    }
    names(results) <- c("estimate", "95CIl", "95CIu", "95PIl", "95PIu")
    
    out$rma <- fit
    out$results.c <- results
  }
  
  return(out)
}

print.valmeta <- function(x, ...) {
  if (!is.null(x$results.c)) {
    cat("Model Results for the c-statistic:\n\n")
    print(x$results.c)
    if (x$num.estimated.var.c > 0)
      cat(paste("\nWarning: For ", x$num.estimated.var.c, " validation(s), the standard error was estimated using method '", x$method.restore.se, "'.\n", sep=""))
  }
  
  
}



