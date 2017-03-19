valmeta <- function(cstat, cstat.se, cstat.95CI, OE, OE.se, OE.95CI, citl, citl.se,
                    N, O, E, Po, Po.se, Pe, t.val, t.ma, method="REML", knha=TRUE, verbose=FALSE, 
                    slab, n.chains = 4, pars, ...) {
  pars.default <- list(hp.mu.mean = 0, 
                       hp.mu.var = 1E6,
                       hp.tau.min = 0,
                       hp.tau.max = 2,
                       hp.tau.sigma = 0.5,
                       hp.tau.dist = "dunif",
                       method.restore.c.se="Newcombe.4",
                       model.cstat = "normal/logit", #Alternative: "normal/identity"
                       model.oe = "normal/log") #Alternative: "poisson/log" or "normal/identity"
  
  if (!missing(pars)) {
    for (i in 1:length(pars)) {
      element <- ls(pars)[i]
      pars.default[[element]] <- pars[[i]]
    }
  }
  
  
  out <- list()
  out$call <- match.call()
  out$method <- method
  class(out) <- "valmeta"
  
  N.studies.OE <- 0
  if (!missing(OE)) {
    N.studies.OE <- length(OE)
  } else if (!missing(E)) {
    N.studies.OE <- length(E)
  } else if (!missing(Pe)) {
    N.studies.OE <- length(Pe)
  }else if (!missing(citl)) {
    N.studies.OE <- length(citl)
  }
  
  t.ma <- ifelse(missing(t.ma), NA, t.ma)
  t.val <- ifelse(missing(t.val), NA, t.val)
  
  

  inv.logit <- function(x) {  if(is.numeric(x)) 1/(1+exp(-x)) else stop("x is not numeric!") }
  logit <- function(x) { log(x/(1-x)) }
  
  # Check if we need to load runjags
  if (method=="BAYES") {
    if (!requireNamespace("runjags", quietly = TRUE)) {
      stop("The package 'runjags' is currently not installed!")
    } 
    if (!requireNamespace("rjags", quietly = TRUE)) {
      stop("The package 'rjags' is currently not installed!")
    } 
  }
  
  if(!missing(cstat)) {
    if (missing(cstat.se) & missing(cstat.95CI)) {
      stop("No sampling error was provided for the c-statistic!")
    }
    if (missing(cstat.95CI)) {
      cstat.95CI <- array(NA, dim=c(length(cstat),2))
    }
    if (is.null(dim(cstat.95CI))) {
      warning("Invalid dimension for 'cstat.95CI', argument ignored.")
      cstat.95CI <- array(NA, dim=c(length(cstat),2))
    }
    if (dim(cstat.95CI)[2] != 2 | dim(cstat.95CI)[1] != length(cstat)) {
      warning("Invalid dimension for 'cstat.95CI', argument ignored.")
      cstat.95CI <- array(NA, dim=c(length(cstat),2))
    }
    if (missing(cstat.se)) {
      cstat.se <- array(NA, dim=length(cstat))
    }
    
    out$cstat <- list()
    out$cstat$method.restore.se <- pars.default$method.restore.c.se 
    out$cstat$model <- pars.default$model.cstat
    class(out$cstat) <- "vmasum"
    
    if(missing(slab)) {
      out$cstat$slab <- paste("Study",seq(1, length(cstat)))
    } else {
      out$cstat$slab <- slab
    }
    
    # Apply necessary data transformations
    if (pars.default$model.cstat == "normal/identity") {
      theta <- cstat
      theta.var <- (cstat.se)**2
      theta.cil <- cstat.95CI[,1]
      theta.ciu <- cstat.95CI[,2]
      theta.var.CI <- ((theta.ciu - theta.cil)/(2*qnorm(0.975)))**2 #Derive from 95% CI
      theta.var <- ifelse(is.na(theta.var), theta.var.CI, theta.var) #Prioritize reported SE
    } else if (pars.default$model.cstat == "normal/logit") {
      theta <- log(cstat/(1-cstat))
      theta.var <- (cstat.se/(cstat*(1-cstat)))**2
      theta.cil <- logit(cstat.95CI[,1])
      theta.ciu <- logit(cstat.95CI[,2])
      theta.var.CI <- ((theta.ciu - theta.cil)/(2*qnorm(0.975)))**2
      theta.var <- ifelse(is.na(theta.var.CI), theta.var, theta.var.CI) #Prioritize variance from 95% CI
    } else {
      stop(paste("No appropriate meta-analysis model defined: '", pars.default$model.cstat, "'", sep=""))
    }
    
    num.estimated.var.c <- 0
    
    # Restore missing standard errors
    if (NA %in% theta.var) {
      # Restore missing estimates of the standard error of the c-statistic using information on c, N and O
      if (!missing(O) & !missing(N)) {
        if (verbose) cat("Attempting to restore missing information on the standard error of the c-statistic\n")
        theta.var.hat <- restore.c.var(cstat=cstat, N.subjects=N, N.events=O, 
                                       restore.method=pars.default$method.restore.c.se, model=pars.default$model.cstat)
        num.estimated.var.c <- length(which(is.na(theta.var) & !is.na(theta.var.hat)))
        theta.var <- ifelse(is.na(theta.var), theta.var.hat, theta.var)
      }
      # Replace remaining missing values in theta.var by very large values
      theta.var <- ifelse(is.na(theta.var), 10e6, theta.var)
    }
    
    #Only calculate 95% CI for which no original values were available
    theta.cil[is.na(theta.cil)] <- (theta+qnorm(0.025)*sqrt(theta.var))[is.na(theta.cil)]
    theta.ciu[is.na(theta.ciu)] <- (theta+qnorm(0.975)*sqrt(theta.var))[is.na(theta.ciu)]
    

    ds <- cbind(theta, sqrt(theta.var), theta.cil, theta.ciu)
    colnames(ds) <- c("theta", "theta.se", "theta.95CIl", "theta.95CIu")
    out$cstat$data <- ds
    out$cstat$num.estimated.var.c <- num.estimated.var.c
    
    if (method != "BAYES") { # Use of rma
      
      # Apply the meta-analysis
      fit <- rma(yi=theta, vi=theta.var, data=ds, method=method, knha=knha, slab=out$cstat$slab, ...) 
      preds <- predict(fit)
      
      results <- as.data.frame(array(NA, dim=c(1,5)))
      if (pars.default$model.cstat == "normal/logit") {
        results <- c(inv.logit(coefficients(fit)), inv.logit(c(preds$ci.lb, preds$ci.ub, preds$cr.lb, preds$cr.ub)))
      } else if (pars.default$model.cstat == "normal/identity") {
        results <- c(coefficients(fit), c(preds$ci.lb, preds$ci.ub, preds$cr.lb, preds$cr.ub))
      } else {
        stop ("Meta-analysis model not implemented yet")
      }
      names(results) <- c("estimate", "95CIl", "95CIu", "95PIl", "95PIu")
      
      out$cstat$rma <- fit
      out$cstat$results <- results
    } else {
      # Perform a Bayesian meta-analysis
      model <- .generateBugsCstat(model=pars.default$model.cstat, pars=pars.default, ...)
      
      # Generate initial values from the relevant distributions
      model.pars <- list()
      model.pars[[1]] <- list(param="mu.tobs", param.f=rnorm, param.args=list(n=1, mean=pars.default$hp.mu.mean, sd=sqrt(pars.default$hp.mu.var)))
      model.pars[[2]] <- list(param="bsTau", param.f=runif, param.args=list(n=1, min=pars.default$hp.tau.min, max=pars.default$hp.tau.max))
      inits <- generateMCMCinits(n.chains=n.chains, model.pars=model.pars)
      
      mvmeta_dat <- list(theta = theta,
                         theta.var = theta.var,
                         Nstudies = length(theta))
      jags.model <- runjags::run.jags(model=model, 
                             monitor = c("mu.tobs", "mu.obs", "pred.obs", "bsTau", "PED"), 
                             data = mvmeta_dat, 
                             n.chains = n.chains,
                             silent.jags = !verbose,
                             inits=inits,
                             ...)
      fit <- jags.model$summaries
      
      
      #Extract PED
      fit.dev <- runjags::extract(jags.model,"PED")
      
      results <- c(fit["mu.obs","Mean"], fit["mu.obs", c("Lower95", "Upper95")], fit["pred.obs", c("Lower95", "Upper95")])
      names(results) <- c("estimate", "95CIl", "95CIu", "95PIl", "95PIu")
      
      out$cstat$runjags <- jags.model
      out$cstat$PED <- sum(fit.dev$deviance)+sum(fit.dev$penalty)
      out$cstat$results <- results
    }
  }
  
  ##### Prepare data for OE ratio
  if(N.studies.OE > 0) {
    if (missing(N)) {
      N <- rep(NA, length=N.studies.OE)
    }
    if (missing(O)) {
      O <- rep(NA, length=N.studies.OE)
    }
    if (missing(E)) {
      E <- rep(NA, length=N.studies.OE)
    }
    if (missing(Po)) {
      Po <- rep(NA, length=N.studies.OE)
    }
    if (missing(Po.se)) {
      Po.se <- rep(NA, length=N.studies.OE)
    }
    if (missing(Pe)) {
      Pe <- rep(NA, length=N.studies.OE)
    }
    if (missing(OE)) {
      OE <- rep(NA, length=N.studies.OE)
    }
    if (missing(OE.se)) {
      OE.se <- rep(NA, length=N.studies.OE)
    }
    if (missing(citl)) {
      citl <- rep(NA, length=N.studies.OE)
    }
    if (missing(citl.se)) {
      citl.se <- rep(NA, length=N.studies.OE)
    }
    if (missing(OE.95CI)) {
      OE.95CI <- array(NA, dim=c(N.studies.OE,2))
    }
    if (is.null(dim(OE.95CI))) {
      warning("Invalid dimension for 'OE.95CI', argument ignored.")
      OE.95CI <- array(NA, dim=c(N.studies.OE,2))
    }
    if (dim(OE.95CI)[2] != 2 | dim(OE.95CI)[1] != N.studies.OE) {
      warning("Invalid dimension for 'OE.95CI', argument ignored.")
      OE.95CI <- array(NA, dim=c(N.studies.OE,2))
    }
    
    # Check if the length of all relevant arguments is consistent
    if (length(unique(c(length(N), length(O), length(E), length(Po), length(Po.se), 
                        length(Pe), length(OE), length(OE.se), length(citl),
                        length(citl.se), dim(OE.95CI)[1]))) > 1) {
      stop("Dimension mismatch")
    }
    
    out$oe$model <- pars.default$model.oe
    class(out$oe) <- "vmasum"
    
    if(missing(slab)) {
      out$oe$slab <- paste("Study",seq(1, N.studies.OE))
    } else {
      out$oe$slab <- slab
    }
    
    ds <- generateOEdata(O=O, E=E, Po=Po, Po.se=Po.se, Pe=Pe, OE=OE, OE.se=OE.se, OE.95CI=OE.95CI, 
                         citl=citl, citl.se=citl.se, N=N, t.ma=t.ma, t.val=t.val,
                         model=pars.default$model.oe)
    out$oe$data <- ds
      
    if (method != "BAYES") { # Use of rma
      
      if (pars.default$model.oe=="normal/identity") {
        fit <- rma(yi=ds$theta, sei=ds$theta.se, data=ds, method=method, knha=knha, slab=out$oe$slab, ...) 
        preds <- predict(fit)
        results <- c(coefficients(fit), c(preds$ci.lb, preds$ci.ub, preds$cr.lb, preds$cr.ub))
      } else if (pars.default$model.oe=="normal/log") {
        fit <- rma(yi=ds$theta, sei=ds$theta.se, data=ds, method=method, knha=knha, slab=out$oe$slab, ...) 
        preds <- predict(fit)
        results <- c(exp(coefficients(fit)), exp(c(preds$ci.lb, preds$ci.ub, preds$cr.lb, preds$cr.ub)))
      } else {
        stop("Model not implemented yet!")
      }
      names(results) <- c("estimate", "95CIl", "95CIu", "95PIl", "95PIu")
      
      out$oe$rma <- fit
      out$oe$results <- results
    } else {
      if (pars.default$model.oe=="normal/identity") {
        stop("Model not implemented yet!")
      } else if (pars.default$model.oe=="normal/log") {
        mvmeta_dat <- list(theta=ds$theta,
                           theta.var=(ds$theta.se)**2,
                           Nstudies = N.studies.OE)
      } else if (pars.default$model.oe =="poisson/log") {
        
        # Truncate hyper parameter variance
        pars.default$hp.mu.var = min(pars.default$hp.mu.var, 100)
        
        mvmeta_dat <- list(obs=O,
                           exc=E,
                           Nstudies = N.studies.OE)
      } else {
        stop("Model not implemented yet!")
      }
      
      # Perform a Bayesian meta-analysis
      model <- generateBugsOE(model=pars.default$model.oe, pars=pars.default, ...)
      
      # Generate initial values from the relevant distributions
      model.pars <- list()
      model.pars[[1]] <- list(param="mu.tobs", param.f=rnorm, param.args=list(n=1, mean=pars.default$hp.mu.mean, sd=sqrt(pars.default$hp.mu.var)))
      model.pars[[2]] <- list(param="bsTau", param.f=runif, param.args=list(n=1, min=pars.default$hp.tau.min, max=pars.default$hp.tau.max))
      inits <- generateMCMCinits(n.chains=n.chains, model.pars=model.pars)
      
      
      jags.model <- runjags::run.jags(model=model, 
                                      monitor = c("mu.tobs", "mu.obs", "pred.obs", "bsTau", "PED"), 
                                      data = mvmeta_dat, 
                                      n.chains = n.chains,
                                      silent.jags = !verbose,
                                      inits=inits,
                                      ...)
      fit <- jags.model$summaries
      

      #Extract PED
      fit.dev <- runjags::extract(jags.model,"PED")
      
      results <- c(fit["mu.obs","Mean"], fit["mu.obs", c("Lower95", "Upper95")], fit["pred.obs", c("Lower95", "Upper95")])
      names(results) <- c("estimate", "95CIl", "95CIu", "95PIl", "95PIu")
      
      out$oe$runjags <- jags.model
      out$oe$PED <- sum(fit.dev$deviance)+sum(fit.dev$penalty)
      out$oe$results <- results
    }
  }
  
  
  
  return(out)
}

.generateBugsCstat <- function(model="normal/logit", #Choose between 'log', 'logit' and 'binom'
                               pars, 
                                ...) # standard deviation for student T prior
  {

  hp.tau.prec <- 1/(pars$hp.tau.sima**2)
  hp.mu.prec <- 1/pars$hp.mu.var
  
  out <- "model {\n " 
  out <- paste(out, "for (i in 1:Nstudies)\n  {\n")
  out <- paste(out, "    theta[i] ~ dnorm(alpha[i], wsprec[i])\n")
  out <- paste(out, "    alpha[i] ~ dnorm(mu.tobs, bsprec)\n")
  out <- paste(out, "    wsprec[i] <- 1/(theta.var[i])\n")
  out <- paste(out, " }\n")
  out <- paste(out, " bsprec <- 1/(bsTau*bsTau)\n")
  
  if (pars$hp.tau.dist=="dunif") {
    out <- paste(out, "  bsTau ~ dunif(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep="") 
  } else if (pars$hp.tau.dist=="dhalft") {
    out <- paste(out, "  bsTau ~ dt(0,", hp.tau.prec, ",3)T(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep="") 
    
  } else {
    stop("Specified prior not implemented")
  }
  
  if (model == "normal/logit") {
    out <- paste(out, "  mu.tobs ~ dnorm(", pars$hp.mu.mean, ",", hp.mu.prec, ")\n", sep="")
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
    cat("Model results for the c-statistic:\n\n")
    print(x$cstat)
    if (x$cstat$num.estimated.var.c > 0)
      cat(paste("\nNote: For ", x$cstat$num.estimated.var.c, " validation(s), the standard error was estimated using method '", x$cstat$method.restore.se, "'.\n", sep=""))
    
    if (!is.null(x$oe$results)) {
      cat("\n\n")
    }
  }
  if (!is.null(x$oe$results)) {
    cat("Model results for the total O:E ratio:\n\n")
    print(x$oe)
  }
}

print.vmasum <- function(x, ...) {
  print(x$results)
  if (!is.null(x$runjags)) {
    #Print penalized expected deviance
    cat(paste("\nPenalized expected deviance: ", round(x$PED,2), "\n"))
    
    # Check if model converged
    psrf.ul <-  x$runjags$psrf$psrf[,2]
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
  if (!is.null(x$cstat)) {
    plotForest(x$cstat, xlab="c-statistic", refline=NULL, ...)
    if (!is.null(x$oe)) {
      readline(prompt="Press [enter] to continue")
    }
  }
  if (!is.null(x$oe)) {
    plotForest(x$oe, xlab="O:E ratio", refline=1, ...)
  }
}

