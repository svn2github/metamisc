valmeta <- function(cstat, cstat.se, cstat.95CI,
                    N, O, method="REML", knha=TRUE, verbose=FALSE, 
                    method.restore.c.se="Newcombe.4", scale.c = "logit", 
                    n.chains = 4,
                    ...) {

  out <- list()
  out$cstat <- list()
  out$call <- match.call()
  out$method <- method
  out$cstat$method.restore.se <- method.restore.c.se 
  out$cstat$scale <- scale.c
  class(out) <- "valmeta"
  class(out$cstat) <- "vmasum"
  
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
  restore.c.var<- function(cstat, N.subjects, N.events, restore.method="Newcombe.4", scale=scale.c) {
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
  
  
  if(!missing(cstat)) {
    # Apply necessary data transformations
    if (scale.c == "identity") {
      theta <- cstat
      theta.var <- (cstat.se)**2
      theta.var.CI <- ((cstat.95CI[,2] - cstat.95CI[,1])/(2*qnorm(0.975)))**2 #Derive from 95% CI
      theta.var <- ifelse(is.na(theta.var), theta.var.CI, theta.var)
    } else if (scale.c == "logit") {
      theta <- log(cstat/(1-cstat))
      theta.var <- (cstat.se/(cstat*(1-cstat)))**2
      theta.var.CI <- ((logit(cstat.95CI[,2]) - logit(cstat.95CI[,1]))/(2*qnorm(0.975)))**2
      theta.var <- ifelse(is.na(theta.var), theta.var.CI, theta.var)
    } else {
      stop("No appropriate transformation defined")
    }
    
    num.estimated.var.c <- 0
    
    # Restore missing standard errors
    if (NA %in% theta.var) {
      # Restore missing estimates of the standard error of the c-statistic using information on c, N and O
      if (!missing(O) & !missing(N)) {
        if (verbose) cat("Attempting to restore missing information on the standard error of the c-statistic\n")
        theta.var.hat <- restore.c.var(cstat=cstat, N.subjects=N, N.events=O, restore.method=method.restore.c.se, scale=scale.c)
        num.estimated.var.c <- length(which(is.na(theta.var) & !is.na(theta.var.hat)))
        theta.var <- ifelse(is.na(theta.var), theta.var.hat, theta.var)
      }
    }
    ds <- cbind(theta, theta.var)
    out$cstat$data <- ds
    out$cstat$num.estimated.var.c <- num.estimated.var.c
    
    if (method != "BAYES") { # Use of rma
      
      # Apply the meta-analysis
      fit <- rma(yi=theta, vi=theta.var, data=ds, method=method, knha=knha, ...) 
      preds <- predict(fit)
      
      results <- as.data.frame(array(NA, dim=c(1,5)))
      if (scale.c == "logit") {
        results <- c(inv.logit(coefficients(fit)), inv.logit(c(preds$ci.lb, preds$ci.ub, preds$cr.lb, preds$cr.ub)))
      } else {
        results <- c(coefficients(fit), c(preds$ci.lb, preds$ci.ub, preds$cr.lb, preds$cr.ub))
      }
      names(results) <- c("estimate", "95CIl", "95CIu", "95PIl", "95PIu")
      
      out$cstat$rma <- fit
      out$cstat$results <- results
    } else {
      # Perform a Bayesian meta-analysis
      model <- .generateBugsCstat(link=scale.c, ...)
      
      mvmeta_dat <- list(theta = theta,
                         theta.var = theta.var,
                         Nstudies = length(theta))
      jags.model <- run.jags(model=model, 
                             monitor = c("mu.tobs", "mu.obs", "pred.obs", "bsTau", "priorTau"), 
                             data = mvmeta_dat, 
                             n.chains = n.chains,
                             ...)
      fit <- jags.model$summaries
      
      results <- c(fit["mu.obs","Mean"], fit["mu.obs", c("Lower95", "Upper95")], fit["pred.obs", c("Lower95", "Upper95")])
      names(results) <- c("estimate", "95CIl", "95CIu", "95PIl", "95PIu")
      
      out$cstat$runjags <- jags.model
      out$cstat$results <- results
    }
  }
  
  return(out)
}

.generateBugsCstat <- function(link="logit", #Choose between 'log', 'logit' and 'binom'
                               prior="dunif", #Choose between dunif (uniform) or dhalft (half student T)
                               prior.bound=c(0,2), #boundaries for uniform prior
                               prior.sigma=0.5, ...) # standard deviation for student T prior
  {

  prior.prec <- 1/(prior.sigma*prior.sigma)
  out <- "model {\n " 
  out <- paste(out, "for (i in 1:Nstudies)\n  {\n")
  out <- paste(out, "    theta[i] ~ dnorm(alpha[i], wsprec[i])\n")
  out <- paste(out, "    alpha[i] ~ dnorm(mu.tobs, bsprec)\n")
  out <- paste(out, "    wsprec[i] <- 1/(theta.var[i])\n")
  out <- paste(out, " }\n")
  out <- paste(out, " bsprec <- 1/(bsTau*bsTau)\n")
  
  if (prior=="dunif") {
    out <- paste(out, "  priorTau ~ dunif(", prior.bound[1], ",", prior.bound[2], ")\n", sep="") 
    out <- paste(out, "  bsTau ~ dunif(", prior.bound[1], ",", prior.bound[2], ")\n", sep="") 
  } else if (prior=="dhalft") {
    out <- paste(out, "  priorTau ~ dt(0,", prior.prec, ",3)T(0,10)\n", sep="") 
    out <- paste(out, "  bsTau ~ dt(0,", prior.prec, ",3)T(0,10)\n", sep="") 
    
  } else {
    stop("Specified prior not implemented")
  }
  
  if (link == "logit") {
    out <- paste(out, "  mu.tobs ~ dnorm(0.0,1.0E-6)\n", sep="")
    out <- paste(out, "  mu.obs <- 1/(1+exp(-mu.tobs))\n", sep="")
    out <- paste(out, "  pred.obs <- 1/(1+exp(-pred.tobs))\n", sep="")
    out <- paste(out, "  pred.tobs ~ dnorm(mu.tobs, bsprec)\n", sep="")
  } else {
    stop("Specified link function not implemented")
  }
  out <- paste(out, "}", sep="")
  return(out)
}

print.valmeta <- function(x, ...) {
  if (!is.null(x$cstat$results)) {
      print(x$cstat)}
}

print.vmasum <- function(x, ...) {
  cat("Model Results for the c-statistic:\n\n")
  print(x$results)
  if (x$num.estimated.var.c > 0)
    cat(paste("\nNote: For ", x$num.estimated.var.c, " validation(s), the standard error was estimated using method '", x$method.restore.se, "'.\n", sep=""))
  if (!is.null(x$runjags)) {
    # Check if model converged
    psrf.ul <-  x$runjags$psrf$psrf[,"97.5% quantile"]
    psrf.target <- x$runjags$psrf$psrf.target
    
    if(sum(psrf.ul > psrf.target)>1) {
      warning(paste("Model did not properly converge! The upper bound of the convergence diagnostic (psrf) exceeds", 
                    psrf.target, "for the parameters", 
              paste(rownames(x$runjags$psrf$psrf)[which(psrf.ul > psrf.target)], " (psrf=", 
                    round(x$runjags$psrf$psrf[which(psrf.ul > psrf.target),2],2), ")", collapse=", ", sep=""),
              ". Consider re-running the analysis by increasing the optional arguments 'adapt', 'burnin' and/or 'sample'."  ))
    }
  }
}

plot.valmeta <- function(x, ...) {
  inv.logit <- function(x) {1/(1+exp(-x)) }
  
  if (!is.null(x$cstat)) {
    if (!is.null(x$cstat$rma)) {
      # Forest plot for the c-statistic
      if (x$cstat$scale=="logit") {
        forest(x$cstat$rma, transf=inv.logit, xlab="c-statistic", addcred=T, ...)
      } else {
        forest(x$cstat$rma, transf=NULL, xlab="c-statistic", addcred=T, ...)
      }
    } else {
      #mu.scaled <- x$cstat$runjags$summary$statistics["mu.tobs","Mean"]
      #slab <- seq(1, nrow(x$cstat$data))

      #rma <- list()
      #class(rma) <- c("rma.uni", "rma")
      #rma$int.only <- T
      #rma$slab.null <- T
      #rma$yi <- rma$yi.f <- x$cstat$data[,"theta"]
      #attr(rma$yi,"slab") <- attr(rma$yi.f,"slab") <- slab
      #attr(rma$yi,"measure") <- attr(rma$yi.f,"measure") <- "GEN"
      #rma$slab <- slab
      #rma$coef.na <- F
      #names(rma$coef.na) <- "X"
      #rma$k.f <- rma$k <- nrow(x$cstat$data)
      #rma$vi <- rma$vi.f <- x$cstat$data[,"theta.var"]
      #rma$X.f <- rma$X <- array(1, dim=c(nrow(x$cstat$data),1))
      #names(rma$X.f) <- names(rma$X) <- c("intrcpt")
      #rma$b <- array(mu.scaled, dim=c(1,1))
      #rownames(rma$b) <- "intrcpt"
      #rma$tau2 <- x$cstat$runjags$summary$statistics["bsTau","Mean"]**2
      #rma$weighted <- T
      #forest(rma, transf=inv.logit, xlab="c-statistic", addcred=T, ...)
      warning("Forest plot not implemented yet for the Bayesian meta-analysis!")
      #create own version of forest plot
    }
  }
}

