# Meta-analysis is a statistical technique by which information from
# independent studies is assimilated. This function allows to perform
# fixed-effects and random-effects meta-analysis.
# r    : Vector of the effect sizes
# vars : Vector of the effect variances
###############################################################################
# Author  : Thomas Debray
# Version : 10 May 2011
###############################################################################
# Example
# example.r = c(0.10,0.30,0.35,0.65,0.45,0.15)
# example.var = c(0.03,0.03,0.05,0.01,0.05,0.02)
# uvmeta(example.r,example.var)
###############################################################################
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

  ds <- as.data.frame(cbind(as.vector(r),as.vector(vars)))
  colnames(ds) <- c("theta","v")
  est <- NA    
  
  if (length(x)!=length(y)) {
    stop("The vectors 'r' and 'vars' have a different length!")
  }
  
  if (missing(na.action)) 
    na.action <- "na.fail"
  if (length(na.action)) 
    ds <- do.call(na.action, list(ds))
  
  
  #############################################################################
  # Start analyses
  #############################################################################
  numstudies = dim(ds)[1]
  dfr = numstudies-1
  
  if(numstudies < 3) {
    warning("There are very few primary studies!")
  }
  
  if (method == "MOM") { 
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
    I_sq = 0
    H_sq = 1    #2
    se_lnH = 0  #SE(ln(H2))
    
    # Between-study variance
    if (Q > dfr) {
      re_C =  sum(w) - sum(w**2)/sum(w)
      between_study_var = (Q - dfr)/re_C
      I_sq = (Q-dfr)/Q
      H_sq = Q/dfr
      se_lnH = (log(Q)-log(dfr))/(2*(sqrt(2*Q)-sqrt((2*length(ds$theta))-3)))
    } else {
      between_study_var = 0
      se_lnH = sqrt((1/(2*(length(ds$theta)-2)))*(1-(1/(3*((length(ds$theta)-2)**2)))))
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
    
    #Q.critical = qchisq(0.95,df=(length(r)-1))
    Q_p = 1-pchisq(Q,df=dfr)
    
    results = as.data.frame(array(NA,dim=c(4, length(pars$quantiles)+2)))
    colnames(results) = c("Mean","Var",paste(pars$quantiles*100,"%",sep=""))
    rownames(results) = c("mu","tausq","Q","Isq")
    results[3,] = c(Q,NA,qchisq(pars$quantiles,df=dfr))
    
    if (model=="random") {
      results[1,] = c(re_weighted_Tbar,re_var_T,re_weighted_Tbar+qnorm(pars$quantiles)*sqrt(re_var_T))
      results[2,] = c(between_study_var,NA,rep(NA,length(pars$quantiles)))
      results[4,] = c(I_sq,NA,rep(NA,length(pars$quantiles)))
    } else if (model=="fixed") {
      results[1,] = c(weighted_Tbar,var_T,weighted_Tbar+qnorm(pars$quantiles)*sqrt(var_T))
      results[2,] = c(0,0,rep(NA,length(pars$quantiles)))
    }
    est <- list(results=results,model=model,df=dfr,numstudies=numstudies)
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
    
    results <- summary(samples,quantiles=pars$quantiles) #summary.mcmc(samples,quantiles=pars$quantiles)
    
    results.overview = as.data.frame(array(NA,dim=c(dim(results[[1]])[1], length(pars$quantiles)+2)))
    colnames(results.overview) = c("Mean","Var",paste(pars$quantiles*100,"%",sep=""))
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
  
  est$na.action <- na.action
    est$call <- match.call()
    class(est) <- "uvmeta"
    return(est)
}



print.uvmeta <- function(x, ...)
{
	out <- (x$results)
	print(out)
	out
}


predict.uvmeta <- function(object, level = 0.95, ...)
{
  alpha = (1-level)/2

  #The correct number of degrees of freedom for this t distribution is complex, and we use a value of k–2 largely for pragmatic reasons. (Riley 2011)
  df = 2 
  
  pred.mean  <- object$results["mu","Mean"]
  pred.lower <- object$results["mu","Mean"] + qt(alpha,df=(object$numstudies-df))*sqrt(object$results["tausq","Mean"]+object$results["mu","Var"])
  pred.upper <- object$results["mu","Mean"] + qt((1-alpha),df=(object$numstudies-df))*sqrt(object$results["tausq","Mean"]+object$results["mu","Var"])
  predint <- c(pred.mean,pred.lower,pred.upper)
  names(predint) <- c("Estimate", paste((alpha*100),"%"),paste(((1-alpha)*100),"%"))
  predint
}

summary.uvmeta <- function(object, ...)
{
    cat("Call:\n")
    print(object$call)
    if (object$model=="fixed")  cat(paste("\nFixed effects summary:\t",round(object$results["mu","Mean"],5))," (SE: ",round(sqrt(object$results["mu","Var"]),5), ")",sep="")
    if (object$model=="random") {
        cat(paste("\nRandom effects summary:\t",round(object$results["mu","Mean"],5))," (SE: ",round(sqrt(object$results["mu","Var"]),5), ")",sep="")
        cat(paste("\n\nTau squared: \t\t",round(object$results["tausq","Mean"],5),sep=""))
    }
    Q_p = 1-pchisq(object$results["Q","Mean"],df=object$df)
    cat(paste("\nCochran's Q statistic: \t",round(object$results["Q","Mean"],5)," (p-value: ",round(Q_p,5),")",sep=""))
    cat(paste("\nI-square index: \t", round(object$results["Isq","Mean"]*100,3)," %\n",sep=""))
}



