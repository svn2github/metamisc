# Meta-analysis is a statistical technique by which information from
# independent studies is assimilated. This function allows to perform
# fixed-effects and random-effects meta-analysis.
# r    : Vector of the effect sizes
# vars : Vector of the effect variances
###############################################################################
# Example
# example.r = c(0.10,0.30,0.35,0.65,0.45,0.15)
# example.var = c(0.03,0.03,0.05,0.01,0.05,0.02)
# uvmeta(example.r,example.var)
###############################################################################

#TODO: allow data transformations
uvmeta <- function(r, vars, model="random", method="MOM", na.action,
                   pars=list(quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975), 
                             n.chains=4, n.adapt=5000, n.init=1000, 
                             n.iter=10000), verbose=FALSE, ...) 
  UseMethod("uvmeta")

uvmeta.default <- function(r,vars, model="random", method="MOM", na.action, 
                           pars=list(quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
                                     n.chains=4, n.adapt=5000, n.init=1000, 
                                     n.iter=10000), verbose=FALSE, ...)
{
  if (length(r)!=length(vars)) {
    stop("The vectors 'r' and 'vars' have different lengths!")
  }

  ds <- as.data.frame(cbind(as.vector(r),as.vector(vars)))
  colnames(ds) <- c("theta","v")
  
  
  if (missing(na.action)) 
    na.action <- "na.fail"
  if (length(na.action)) 
    ds <- do.call(na.action, list(ds))

  est <- NA 
  
  
  #############################################################################
  # Start analyses
  #############################################################################
  numstudies = dim(ds)[1]
  dfr = numstudies-1
  
  if(numstudies < 3) {
    warning("There are very few primary studies!")
  }
  
  if (method == "MOM") { 
    results = as.data.frame(array(NA,dim=c(4, length(pars$quantiles)+2)))
    colnames(results) = c("Estimate","Var",paste(pars$quantiles*100,"%",sep=""))
    rownames(results) = c("mu","tausq","Q","Isq")
    
    # FIXED EFFECTS MODEL
    w = 1/ds$v
    
    #Combined effect
    weighted_Tbar = sum(ds$theta*w)/sum(w)
    
    # Variance of the combined effect
    var_T = 1/sum(w)
    
    # Standard error combined effect
    se_T = sqrt(var_T)
    
    # The Z-value
    z_T = weighted_Tbar/se_T
    
    # RANDOM EFFECTS MODEL
    Q = sum(w*(ds$theta-weighted_Tbar)**2)
    results["Q",] = c(Q,NA,qchisq(pars$quantiles,df=dfr))
    
    
    # Between-study variance
    if (model=="random" & Q > dfr) {
      re_C =  sum(w) - sum(w**2)/sum(w)
      between_study_var = (Q - dfr)/re_C
    } else {
      between_study_var = 0
    }

    # Within-study plus between-study variance
    re_v = vars + between_study_var
    
    # Updated weights
    re_w = 1/re_v
    
    # Combined effect
    re_weighted_Tbar =  sum(ds$theta*re_w)/sum(re_w)
    
    # Variance of the combined effect
    re_var_T  = 1/sum(re_w)
    
    # Standard error of combined effect
    re_se_T = sqrt(re_var_T)
    
    # The Z-value
    re_z_T = re_weighted_Tbar/re_se_T
    
    
    
    if (model=="random") {
      results["mu",] = c(re_weighted_Tbar,re_var_T,re_weighted_Tbar+qnorm(pars$quantiles)*sqrt(re_var_T))
      results["tausq",] = c(between_study_var,NA,rep(NA,length(pars$quantiles)))
      
      # Calculate I2 and its confidence limits
      Isq <- (results["Q",]-dfr)/results["Q",]
      Isq[which(Isq>1)] <- 1
      Isq[which(Isq<0)] <- 0
      results["Isq",] = Isq
    } else if (model=="fixed") {
      results["mu",] = c(weighted_Tbar,var_T,weighted_Tbar+qnorm(pars$quantiles)*sqrt(var_T))
      results["tausq",] = c(0,0,rep(NA,length(pars$quantiles)))
    }
    pred.int <- results["mu","Estimate"] + qt(pars$quantiles,df=(numstudies-2))*sqrt(results["tausq","Estimate"]+results["mu","Var"])
    names(pred.int) <- paste(pars$quantiles*100,"%",sep="")
    
    est <- list(results=results,model=model,df=dfr,numstudies=numstudies, pred.int=pred.int)
    
  } else if (method == "bayes") { 
    quiet = !verbose
    
    modelfile <-  if (model=="random") system.file(package="metamisc", "model", "uvmeta_ranef.bug") else system.file(package="metamisc", "model", "uvmeta_fixef.bug")
    jags <- jags.model(modelfile,
                       data = list('r' = ds$theta,
                                   'vars' = ds$v,
                                   'k' = numstudies), #prior precision matrix
                       n.chains = pars$n.chains,
                       n.adapt = pars$n.adapt,
                       quiet = quiet)
    update(jags, pars$n.init) #initialize
    samples <- coda.samples(jags, c('mu','tausq','Q','Isq'),n.iter=pars$n.iter)
    
    results <- summary(samples,quantiles=pars$quantiles) 
    
    #TODO: calculate prediction interval
    
    results.overview = as.data.frame(array(NA,dim=c(dim(results[[1]])[1], length(pars$quantiles)+2)))
    colnames(results.overview) = c("Estimate","Var",paste(pars$quantiles*100,"%",sep=""))
    rownames(results.overview) = rownames(results[[2]])
    results.overview[,1] = (results[[1]])[,"Mean"]
    results.overview[,2] = (results[[1]])[,"SD"]**2
    for (i in 1:length(pars$quantiles)) {
      results.overview[,(i+2)] = (results[[2]])[,i]
    }
    
    est <- list(results=results.overview,model=model,df=dfr,numstudies=numstudies)
  } else {
    stop("Invalid meta-analysis method!")
  }
  # } else if (method=="pl") {
  #   results = as.data.frame(array(NA,dim=c(4, length(pars$quantiles)+2)))
  #  colnames(results) = c("Estimate","Var",paste(pars$quantiles*100,"%",sep=""))
  #  rownames(results) = c("mu","tausq","Q","Isq")
  #  
  #  mle.loglik <- function(theta, tausq, ds) {
  #    loglik <- sum(dnorm(x=ds$theta,mean=theta, sd=sqrt(tausq+ds$v),log=T))
  #    return (-loglik)
  #  }
  #      
  #  #### 5.99 ===> qchisq (0.95,df=2)
  #  if (model == "random") 
  #  {
  #    mle <- mle2(minuslogl=mle.loglik, start=list(theta=0,tausq=0), data=list(ds=ds),method="L-BFGS-B",lower=list(theta=-Inf,tausq=0))
  #    mle.cov <- vcov(mle)
  #    p0 <- profile(mle)
  #    
  #    levels = pars$quantiles
  #    levels[which(pars$quantiles<0.5)]  = 1-(pars$quantiles[which(pars$quantiles<0.5)]*2)
  #    levels[which(pars$quantiles>=0.5)] = 1-(1-pars$quantiles[which(pars$quantiles>=0.5)])*2
  #    pci = array(NA,dim=c(2,length(levels)))
  #    colnames(pci) = paste(pars$quantiles*100,"%",sep=" ")
  #    
  #    for (i in 1:length(levels)) {
  #      pcint <- confint(p0,level=levels[i])
  #      cols.select <- which(colnames(pcint) %in% colnames(pci))
  #      pci[,colnames(pcint)[cols.select]] <- pcint[,cols.select]
  #    }
  #    
  #    wt <- 1/(coef(mle)["tausq"]+ds$v)
  #    Q <- sum(wt*(ds$theta-coef(mle)["theta"])**2)
  #    results["Q",] = c(Q,NA,qchisq(pars$quantiles,df=dfr))
  #    
  #    
  #    # Use profile log-likelihood to calculate confidence intervals
  #    results["mu",] = c(coef(mle)["theta"],mle.cov[1,1],pci[1,])
  #    results["tausq",] = c(mle.tausq,mle.cov[2,2],pci[2,])
  #    
  #    # Calculate I2 and its confidence limits
  #    Isq <- (results["Q",]-dfr)/results["Q",]
  #    Isq[which(Isq>1)] <- 1
  #    Isq[which(Isq<0)] <- 0
  #    results["Isq",] = Isq
  # }
  #  pred.int <- results["mu","Estimate"] + qt(pars$quantiles,df=(numstudies-2))*sqrt(results["tausq","Estimate"]+results["mu","Var"])
  #  names(pred.int) <- paste(pars$quantiles*100,"%",sep="")

  #  est <- list(results=results,model=model,df=dfr,numstudies=numstudies, pred.int=pred.int)
 
  

  est$na.action <- na.action
  est$method <- method
  est$call <- match.call()
  class(est) <- "uvmeta"
  return(est)
}



print.uvmeta <- function(x, ...)
{
  out <- (x$results)
  text.model <- if (x$model=="fixed") "Fixed" else "Random"
  text.method <- if(x$method=="bayes") "credibility" else "confidence"
  cat(paste(text.model,"effects estimates with corresponding", text.method, "intervals:\n\n"))
	print(out)
  if (x$model=="random") {
    cat(paste("\n\nPrediction interval for mu:\n\n"))
    print(x$pred.int)
  }
  
	out
}


summary.uvmeta <- function(object, ...)
{
    cat("Call:\n")
    print(object$call)
    if (object$model=="fixed")  cat(paste("\nFixed effects summary:\t",round(object$results["mu","Estimate"],5))," (SE: ",round(sqrt(object$results["mu","Var"]),5), ")",sep="")
    if (object$model=="random") {
        cat(paste("\nRandom effects summary:\t",round(object$results["mu","Estimate"],5))," (SE: ",round(sqrt(object$results["mu","Var"]),5), ")",sep="")
        cat(paste("\n\nTau squared: \t\t",round(object$results["tausq","Estimate"],5),sep=""))
    }
    Q_p = 1-pchisq(object$results["Q","Estimate"],df=object$df)
    cat(paste("\nCochran's Q statistic: \t",round(object$results["Q","Estimate"],5)," (p-value: ",round(Q_p,5),")",sep=""))
    cat(paste("\nI-square index: \t", round(object$results["Isq","Estimate"]*100,3)," %\n",sep=""))
}



