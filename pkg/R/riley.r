riley <- function(X, ...) UseMethod("riley")

# Estimated parameters
# mu1 : logit(sensitivity)
# mu2 : logit(specificity)
# psi1 : standard deviation of mu1
# psi2 : standard deviation of mu2
# rhoT : (logit(rho+1)/2) where rho is the correlation between mu1 and mu2
riley.default <-
  function(X = NULL, TP, FN, FP, TN, correction = 0.5, 
           correction.control = "all", ...)
  {
      if(!is.null(X)){
        X <- as.data.frame(X)
        origdata <- X
        TP <- X$TP
        FN <- X$FN
        FP <- X$FP
        TN <- X$TN
      }
      
      ## apply continuity correction to _all_ studies if one contains zero
      if(correction.control == "all"){if(any(c(TP,FN,FP,TN) == 0)){TP <- TP + correction;
                                                                   FN <- FN + correction;
                                                                   FP <- FP + correction;
                                                                   TN <- TN + correction}}
      if(correction.control == "single"){
        correction = ((((TP == 0)|(FN == 0))|(FP == 0))| (TN == 0))*correction
        TP <- correction + TP
        FN <- correction + FN
        FP <- correction + FP
        TN <- correction + TN
      }
      
      # Numstudies
      numstudies = length(TP)
      df <- numstudies - 5  
      if(df < 0){warning("There are very few primary studies!")}
      
     
      #Calculate sensitivities and specificities (original scale)
      number.of.pos <- TP + FN
      number.of.neg <- FP + TN
      sens<-TP/number.of.pos
      fpr <- FP/number.of.neg
      var.sens = sens*(1-sens)/number.of.pos
      var.fpr = fpr*(1-fpr)/number.of.neg
      
      logit.sens <- logit(sens)
      logit.fpr <- logit(fpr)
      var.logit.sens <- 1/(sens*(1-sens)*number.of.pos)
      var.logit.fpr <- 1/(fpr*(1-fpr)*number.of.neg)
      
      vars = cbind(var.logit.sens,var.logit.fpr)
      
      Y = array(NA,dim=c((length(logit.sens)*2),1))
      for (i in 1:length(logit.sens))
      {
        Y[((i-1)*2+1)] = logit.sens[i]
        Y[((i-1)*2+2)] = logit.fpr[i]
      }

      #Calculate starting values for optim
      sumlsens <- uvmeta(r=logit.sens, v=var.logit.sens, method="MOM")
      sumlfpr  <- uvmeta(r=logit.fpr, v=var.logit.fpr, method="MOM")

      negfullloglik <- function(pars,Y,vars)
      {
        beta1 = pars[1]
        beta2 = pars[2]
        psisq1 = pars[3]**2 #ensure variance is positive
        psisq2 = pars[4]**2 #ensure variance is positive
        rho = inv.logit(pars[5])*2-1 #ensure correlation is in [-1,1], and values in that interval move symmetric from -1 to 0 and from 1 to 0
        k = 2 #2 endpoints
        n = dim(Y)[1]/2
        
        #Beta vector
        Beta = rbind(beta1,beta2)
        
        #Design matrix
        X = array(NA,dim=c(n*2,2))
        X[,1] = rep(c(1,0),n)
        X[,2] = rep(c(0,1),n)
        
        #Create Phi matrix
        Phi = array(0,dim=c((n*2),(n*2)))
        for (i in 1:n) {
          Phi[((i-1)*2+1),((i-1)*2+1)] = vars[i,1]+psisq1
          Phi[((i-1)*2+2),((i-1)*2+2)] = vars[i,2]+psisq2
          Phi[((i-1)*2+1),((i-1)*2+2)] = rho*sqrt((vars[i,1]+psisq1)*(vars[i,2]+psisq2))
          Phi[((i-1)*2+2),((i-1)*2+1)] = rho*sqrt((vars[i,1]+psisq1)*(vars[i,2]+psisq2))
        }
        
        #Minimize the negative of the restricted log-lkh
        0.5*((n-k)*log(2*pi)-log(det(t(X)%*%X))+log(det(Phi))+log(det(t(X)%*%solve(Phi)%*%X))+(t(Y-X%*%Beta)%*%solve(Phi)%*%(Y-X%*%Beta)))
      }
      
      pars.start = c(sumlsens$ranef$mean,sumlfpr$ranef$mean,sqrt(sumlsens$ranef$var),sqrt(sumlfpr$ranef$var),0)
      fit = optim(pars.start,negfullloglik,Y=Y,vars=vars,hessian=T)
      
      mu1 = fit$par[1]
      mu2 = fit$par[2]
      psi1 = abs(fit$par[3])
      psi2 = abs(fit$par[4])
      #rho = inv.logit(fit$par[5])*2-1
      rhoT = fit$par[5]
      
      coefficients = c(mu1,mu2,psi1,psi2,rhoT)
      names(coefficients) = c("mu1","mu2","psi1","psi2","rhoT")
      
      Sigma = solve(fit$hessian)
      colnames(Sigma) = c("mu1","mu2","psi1","psi2","rhoT")
      rownames(Sigma) = c("mu1","mu2","psi1","psi2","rhoT")
      
      iterations <- fit$iterations
      logLik <- -fit$value
      
      output <- list(coefficients = coefficients, vcov = Sigma, df = df, nobs = 2*numstudies, logLik = -logLik,
                   iterations = (iterations+1), call = match.call(), data = origdata)
      class(output) <- "riley"
      return(output)
  }

print.riley <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
}

# Calculate prediction interval (not identical interpretation to random effects!)
predict.riley <- function(object, level = 0.95, ...)
{
  alpha = (1-level)/2
  
  #Prediction interval
  predint=array(NA,dim=c(2,3))
  rownames(predint) = c("Sens","FPR")
  colnames(predint) = c("Estimate", paste((alpha*100),"%"),paste(((1-alpha)*100),"%"))
  df = object$df
  if (df > 0)
  {
    predint[1,] = inv.logit(c(coefficients(object)["mu1"],(qt(alpha,df=df)*sqrt((coefficients(object)["psi1"]**2)+diag(vcov(object))[1])+coefficients(object)["mu1"]),(qt((1-alpha),df=df)*sqrt((coefficients(object)["psi1"]**2)+diag(vcov(object))[1])+coefficients(object)["mu1"])))
    predint[2,] = inv.logit(c(coefficients(object)["mu2"],(qt(alpha,df=df)*sqrt((coefficients(object)["psi2"]**2)+diag(vcov(object))[2])+coefficients(object)["mu2"]),(qt((1-alpha),df=df)*sqrt((coefficients(object)["psi2"]**2)+diag(vcov(object))[2])+coefficients(object)["mu2"])))
  } else {
    predint[1,] = c(inv.logit(mu1),0,1)
    predint[2,] = c(inv.logit(mu2),0,1)
  }
  predint
}


summary.riley <- function(object, level = 0.95, ...)
{
  confints <- cbind(object$coefficients, confint(object,level=level))
  colnames(confints)[1] <- "Estimate"
  confints <- rbind(inv.logit(confints[1:2,]),confints)
  rownames(confints)[1:2] <-  c("Sens", "FPR") 
  #Transform last parameter back to rho
  confints[7,] =  inv.logit(confints[7,])*2-1
  rownames(confints)[7] = "rho"
  
  res <- list(call=object$call, confints = confints)
  class(res) <- "summary.riley"
  res
}

print.summary.riley <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  print(x$confints)
}

vcov.riley <- function(object, ...){
  object$vcov
  }

logLik.riley <- function(object, ...){object$logLik}


plot.riley <- function(x, plotsumm = TRUE, plotnumerics = TRUE, level = 0.95, main="",
                       ylim = c(0,1), xlim = c(0,1), pch = 1, lty = 1, lwd = 1, cex.numerics=0.45,
                       add=F, ...)
{

  if (!add) plot(-500,-500, type = "l", xlim = xlim, ylim = ylim, xlab="1-Specificity",ylab="Sensitivity",main=main)
  
  alpha = (1-level)/2
  FP <- x$data$FP
  negatives <- FP + x$data$TN
  FPR <- FP/negatives
  mu = x$coefficients[c("mu1","mu2")]
  Sigma = x$vcov[c(1,2),c(1,2)] 
  mu.ellipse <- ellipse(Sigma, centre = mu, level = level)           
  polygon(inv.logit(mu.ellipse[,2]),inv.logit(mu.ellipse[,1]),lty=lty, lwd=lwd)
  
  if(plotsumm) points(inv.logit(mu[2]),inv.logit(mu[1]),pch=pch) # add the point estimate of the mean
  
  if(plotnumerics) {
    ci = summary(x,level=level)[2]$confints
    text(0.8,0.15,labels="Estimate",pos=2,cex=cex.numerics)
    text(0.9,0.15,labels=paste((alpha*100),"% CI",sep=""),pos=2,cex=cex.numerics)
    text(1.0,0.15,labels=paste(((1-alpha)*100),"% CI",sep=""),pos=2,cex=cex.numerics)
    text(0.5,0.10,labels= "Sensitivity",pos=4, cex=cex.numerics)
    text(0.8,0.10,labels=paste("",formatC(round( ci["Sens",1],2),2,format="f",flag="0")),pos=2, cex=cex.numerics)
    text(0.9,0.10,labels=paste("",formatC(round( ci["Sens",2],2),2,format="f",flag="0")),pos=2, cex=cex.numerics)
    text(1.0,0.10,labels=paste("",formatC(round( ci["Sens",3],2),2,format="f",flag="0")),pos=2, cex=cex.numerics)
    text(0.5,0.05,labels= "1-Specificity",pos=4, cex=cex.numerics)
    text(0.8,0.05,labels=paste("",formatC(round( ci["FPR",1],2),2,format="f",flag="0")),pos=2, cex=cex.numerics)
    text(0.9,0.05,labels=paste("",formatC(round( ci["FPR",2],2),2,format="f",flag="0")),pos=2, cex=cex.numerics)
    text(1.0,0.05,labels=paste("",formatC(round( ci["FPR",3],2),2,format="f",flag="0")),pos=2, cex=cex.numerics)
  }
}
