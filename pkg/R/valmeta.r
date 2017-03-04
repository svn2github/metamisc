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
      # Replace remaining missing values in theta.var by very large values
      theta.var <- ifelse(is.na(theta.var), 10e6, theta.var)
    }
    theta.cil <- theta+qnorm(0.025)*sqrt(theta.var)
    theta.ciu <- theta+qnorm(0.975)*sqrt(theta.var)
    ds <- cbind(theta, sqrt(theta.var), theta.cil, theta.ciu)
    colnames(ds) <- c("theta", "theta.se", "theta.95CIl", "theta.95CIu")
    out$cstat$data <- ds
    out$cstat$slab <- paste("Study",seq(1, length(theta)))
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
                             monitor = c("mu.tobs", "mu.obs", "pred.obs", "bsTau", "priorTau", "PED"), 
                             data = mvmeta_dat, 
                             n.chains = n.chains,
                             silent.jags = !verbose,
                             ...)
      fit <- jags.model$summaries
      
      
      #Extract PED
      fit.dev <- extract(jags.model,"PED")
      
      results <- c(fit["mu.obs","Mean"], fit["mu.obs", c("Lower95", "Upper95")], fit["pred.obs", c("Lower95", "Upper95")])
      names(results) <- c("estimate", "95CIl", "95CIu", "95PIl", "95PIu")
      
      out$cstat$runjags <- jags.model
      out$cstat$PED <- sum(fit.dev$deviance)+sum(fit.dev$penalty)
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
  cat("Model results for the c-statistic:\n\n")
  print(x$results)
  if (!is.null(x$runjags)) {
    #Print penalized expected deviance
    cat(paste("\nPenalized expected deviance: ", round(x$PED,2), "\n"))
    
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
  if (x$num.estimated.var.c > 0)
    cat(paste("\nNote: For ", x$num.estimated.var.c, " validation(s), the standard error was estimated using method '", x$method.restore.se, "'.\n", sep=""))
  
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
      col <- c("black", "gray50")
      border <- "black"
      lty <- c("solid", "dotted", "solid")
      cex <- 0.8
      efac <- 1
      efac <- rep(efac, 2)
      xlim <- c(-0.5, 1.5)
      
      par.usr <- par("usr")
      height <- par.usr[4] - par.usr[3]
      
      k <- dim(x$cstat$data)[1]
      slab <- c(x$cstat$slab, "RE Model")
      yi <- x$cstat$data[,"theta"]
      ci.lb <- x$cstat$data[,"theta.95CIl"]
      ci.ub <- x$cstat$data[,"theta.95CIu"]
      
      if (x$cstat$scale=="logit") {
        yi <- sapply(yi, inv.logit)
        ci.lb <- sapply(ci.lb, inv.logit)
        ci.ub <- sapply(ci.ub, inv.logit)
      }
      
      #Add the meta-analysis summary to the results
      #Note that no transormations are needed here, as summaries are always presented on original scale
      yi <- c(yi, x$cstat$results["estimate"])
      ci.lb <- c(ci.lb, x$cstat$results["95CIl"])
      ci.ub <- c(ci.ub, x$cstat$results["95CIu"])
      
      rows <- c(seq(k,1),-1)
      
      annotext <- round(cbind(yi, ci.lb, ci.ub), 2)
      annotext <- matrix(apply(annotext, 2, format, nsmall = 2), ncol = 3)
      annotext <- paste(annotext[,1], "[", annotext[,2], ",", annotext[,3], "]")
      
      
      par.mar <- par("mar")
      par.mar.adj <- par.mar - c(0, 3, 1, 1)
      par.mar.adj[par.mar.adj < 0] <- 0
      par(mar = par.mar.adj)
      on.exit(par(mar = par.mar))
      
      par.usr <- par("usr")
      height <- par.usr[4] - par.usr[3]
      lheight <- strheight("O")
      cex.adj <- ifelse(k * lheight > height * 0.8, height/(1.25 * k * lheight), 1)
      cex <- par("cex") * cex.adj
      
      plot(NA, NA, xlim=xlim, ylim=c(-2,k), ylab="", xlab="c-statistic",yaxt = "n", xaxt = "n", xaxs = "i", bty = "n", ...)
      for (i in 1:k) {
        points(yi[i], rows[i], pch = 15, ...)
        
        segments(ci.lb[i], rows[i], ci.ub[i], rows[i], ...)
        
        segments(ci.lb[i], rows[i] - (height/150) * cex * 
                   efac[1], ci.lb[i], rows[i] + (height/150) * cex * 
                   efac[1], ...)
        
        segments(ci.ub[i], rows[i] - (height/150) * cex * 
                   efac[1], ci.ub[i], rows[i] + (height/150) * cex * 
                   efac[1], ...)
      }
      text(xlim[1], rows, slab, pos = 4, cex = cex, ...)
      text(x = xlim[2], rows, labels = annotext, pos = 2, cex = cex, ...)
      
      # Add meta-analysis summary
      abline(h = 0, lty = 1, ...)

      axis(side = 1, at = c(0,0.2,0.4,0.6,0.8,1), labels = c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis = 1, ...)
    }
  }
}

