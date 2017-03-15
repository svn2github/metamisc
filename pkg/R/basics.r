logit <- function(x) { if(is.numeric(x))  log(x/(1-x)) else stop("x is not numeric!") }

inv.logit <- function(x) {  if(is.numeric(x)) 1/(1+exp(-x)) else stop("x is not numeric!") }

generateMCMCinits <- function(n.chains, model.pars)
{
  inits <- list()
  for (i in 1:n.chains) {
    inits.i <- list()
    for (j in 1:length(model.pars)) {
      parname <- model.pars[[j]]$param
      fprior <- model.pars[[j]]$param.f
      fargs <- model.pars[[j]]$param.args
      inits.i[[parname]] = do.call(fprior, fargs)
    }
    inits[[i]] <- inits.i
  }
  return(inits)
}


#Update SE(c.index) using Method 4 of Newcombe
restore.c.var<- function(cstat, N.subjects, N.events, restore.method="Newcombe.4", scale="logit") {
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
  } else if (scale=="identity") {
    out <- ((cstat*(1-cstat)*(1+nstar*(1-cstat)/(2-cstat) + mstar*cstat/(1+cstat)))/(m*n))
  } else {
    stop ("Scale not implemented!")
  }
  
  return(out)
}

restore.oe.var <- function(citl, citl.se, Po) {
  nom <- ((Po-1)**2)*((Po**2)+1)*((exp(Po+citl))**2)*(citl.se**2)
  denom <- (Po*(-exp(citl))+Po+exp(citl))**2
  out <- nom/denom
  return(out)
}

extrapolateOE <- function(Po, Pe, var.Po, t.val, t.ma, N, scale="log") {
  Po.new <- 1-exp(t.ma*log(1-Po)/t.val)
  Pe.new <- 1-exp(t.ma*log(1-Pe)/t.val)
  theta <- theta.var <- NA
  
  if (scale=="identity") {
    theta <- Po.new/Pe.new
    
    if (missing(var.Po)) {
      # Approximate SE of Kaplan-Meier using binomial distribution
      theta.var <- ((t.ma**2)*exp(2*t.ma*log(1-Po)/t.val)*Po)/((t.val**2)*N*(1-Po)*(1-exp(t.ma*log(1-Pe)/t.val))**2)
    } else {
      # Use error variance of Po, which is equal to the error variance of the KM estimate
      theta.var <- ((t.ma**2)*exp(2*t.ma*log(1-Po)/t.val)*var.Po)/((t.val**2)*((1-Po)**2)*(1-exp(t.ma*log(1-Pe)/t.val))**2)
    }
  } else if (scale=="log") {
    theta <- log(Po.new/Pe.new)
    
    if (missing(var.Po)) {
      # Approximate SE of Kaplan-Meier using binomial distribution
      theta.var <- ((t.ma**2)*Po*exp(2*t.ma*log(1-Po)/t.val))/((t.val**2)*N*(1-Po)*((1-exp(t.ma*log(1-Po)/t.val))**2))
    } else {
      # Use error variance of Po, which is equal to the error variance of the KM estimate
      theta.var <- ((t.ma**2)*exp(2*t.ma*log(1-Po)/t.val)*var.Po)/((t.val**2)*((1-Po)**2)*((1-exp(t.ma*log(1-Po)/t.val))**2))
    }
  } else {
    stop ("Scale not implemented")
  }
  out <- c(theta, theta.var)
  return(out)
}

