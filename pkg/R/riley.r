#TODO: incoroporate level in meta-analysis, and omit it from plots and summaries
#TODO: use generic forest method for generating plots
#TODO: check prediction intervals
#TODO: Add cholesky decomposition in loglikelihood
#TODO: Alter design matrix in LogLik to set entries with missing data to zero (in which case #df needs to be altered)


#' Fit the alternative model for bivariate random-effects meta-analysis
#' 
#' This function fits the alternative model for bivariate random-effects meta-analysis when the within-study correlations 
#' are unknown. This bivariate model was proposed by Riley et al. (2008) and is similar to the general bivariate 
#' random-effects model (van Houwelingen et al. 2002), but includes an overall correlation parameter rather than 
#' separating the (usually unknown) within- and between-study correlation. As a consequence, the alternative model 
#' is not fully hierarchical, and estimates of additional variation beyond sampling error (\code{psi}) are not 
#' directly equivalent to the between-study variation (\code{tau}) from the general model. This model is particularly 
#' useful when there is large within-study variability, few primary studies are available or the general model 
#' estimates the between-study correlation as 1 or -1. 
#' 
#' @param  X data frame containing integer variables \code{Y1}, \code{vars1}, \code{Y2} and \code{vars2},
#' where the columns \code{Y1} and \code{Y2} represent the effect sizes of outcome 1 and, respectively, outcome 2. The columns 
#' \code{vars1} and \code{vars2} represent the error variances of \code{Y1} and, respectively, \code{Y2}. Alternatively, 
#' when considering a meta-analysis of diagnostic test accuracy data, the columns \code{TP}, \code{FN}, \code{FP} and 
#' \code{TN} may be specified. Corresponding values then represent the number of true positives, the number of false negatives,
#' the number of false positives and, respectively, the number of true negatives.
#' @param optimization The optimization method that should be used for minimizing the negative (restricted) 
#' log-likelihood function. The default method is an implementation of that of Nelder and Mead (1965), 
#' that uses only function values and is robust but relatively slow. Other methods are described in \link[stats]{optim}.
#' @param control A list of control parameters to pass to \link[stats]{optim}.
#' @param pars List with additional arguments. The width of confidence, credibility and prediction intervals is 
#' defined by \code{level} (defaults to 0.95). 
#' @param \dots Arguments to be passed on to other functions. See "Details" for more information.
#' 
#' @details Parameters are estimated by iteratively maximizing the restriced log-likelihood using the Newton-Raphson procedure. 
#' The results from a univariate random-effects meta-analysis with a method-of-moments estimator are used as starting 
#' values for \code{beta1}, \code{beta2}, \code{psi1} and \code{psi2} in the \code{optim} command. 
#' The starting value for the transformed correlation \code{rhoT} is 0. Standard errors for all parameters are obtained 
#' from the inverse Hessian matrix. The original correlation is given as \code{inv.logit(rhoT)*2-1}. 
#' 
#' \subsection{Meta-analysis of effect sizes}{
#' The following parameters are estimated by iteratively maximizing the restriced log-likelihood using the Newton-Raphson 
#' procedure: pooled effect size for outcome 1 (\code{beta1}), pooled effect size for outcome 2 (\code{beta2}), 
#' additional variation of \code{beta1} beyond sampling error (\code{psi1}), additional variation of \code{beta2} 
#' beyond sampling error (\code{psi2}) and the (transformed) correlation \code{rhoT} between \code{psi1} and \code{psi2}. 
#' 
#' }
#' 
#' \subsection{Meta-analysis of diagnostic test accuracy}{
#' Although the model can also be used for diagnostic test accuracy data when substantial within-study correlations 
#' are expected, assuming zero within-study correlations (i.e. applying Reitsma's approach) is usually justified 
#' (Reitsma et al. 2005, Daniels and Hughes 1997, Korn et al. 2005, Thompson et al. 2005, Van Houwelingen et al. 2002).
#' 
#' A transformation is applied to the sensitivities ans false positive rates of each study, in order to meet the normality 
#' assumptions. When zero cell counts occur, continuity corrections may be required. The correction value can be specified using
#' \code{correction} (defaults to 0.5). Further, when the argument \code{correction.control} is set to \code{"all"} 
#' (the default) the continuity correction is added to the whole data if only one cell in one study is zero. 
#' If \code{correction.control="single"} the correction is only applied to rows of the data which have a zero.
#' 
#' The following parameters are estimated: logit of sensitivity (\code{beta1}), logit of false positive rate (\code{beta2}), 
#' additional variation of \code{beta1} beyond sampling error (\code{psi1}), additional variation of \code{beta2} beyond 
#' sampling error (\code{psi2}) and the (transformed) correlation \code{rhoT} between \code{psi1} and \code{psi2}. 
#' }
#' 
#' @note Algorithms for dealing with missing data are currently not implemented, but Bayesian approaches will become 
#' available in later versions.
#' 
#' @references
#' \itemize{
#' \item Korn EL, Albert PS, McShane LM. Assessing surrogates as trial endpoints using mixed models. 
#' \emph{Statistics in Medicine} 2005; \bold{24}: 163--182.
#' \item Nelder JA, Mead R. A simplex algorithm for function minimization. \emph{Computer Journal} (1965); \bold{7}: 308--313.
#' \item Reitsma J, Glas A, Rutjes A, Scholten R, Bossuyt P, Zwinderman A. Bivariate analysis of sensitivity and 
#' specificity produces informative summary measures in diagnostic reviews. \emph{Journal of Clinical Epidemiology} 2005; 
#' \bold{58}: 982--990.
#' \item Riley RD, Thompson JR, Abrams KR. An alternative model for bivariate random-effects meta-analysis when 
#' the within-study correlations are unknown. \emph{Biostatistics} 2008; \bold{9}: 172--186.
#' \item Thompson JR, Minelli C, Abrams KR, Tobin MD, Riley RD. Meta-analysis of genetic studies using mendelian 
#' randomization--a multivariate approach. \emph{Statistics in Medicine} 2005; \bold{24}: 2241--2254.
#' \item van Houwelingen HC, Arends LR, Stijnen T. Advanced methods in meta-analysis: multivariate approach and 
#' meta-regression. \emph{Statistics in Medicine} 2002; \bold{21}: 589--624.
#' }
#' 
#' @examples 
#' data(Scheidler)
#' data(Daniels)
#' data(Kertai)
#' 
#' # Meta-analysis of potential surrogate markers data
#' # The results obtained by Riley (2008) were as follows:
#' # beta1 = -0.042 (SE = 0.063),
#' # beta2 = 14.072 (SE = 4.871)
#' # rho   = -0.759
#' \dontrun{
#' fit1 <- riley(Daniels) #maxit reached, try again with more iterations
#' }
#' fit1 <- riley(Daniels, control=list(maxit=10000))
#' summary(fit1)
#' 
#' # Meta-analysis of prognostic test studies
#' fit2 <- riley(Kertai)
#' fit2
#' 
#' # Meta-analysis of computed tomography data 
#' ds <- Scheidler[which(Scheidler$modality==1),]
#' fit3 <- riley(ds)
#' fit3
#' 
#' @return An object of the class \code{riley} for which many standard methods are available. A warning message is 
#' casted when the Hessian matrix contains negative eigenvalues, which implies that the identified solution is a 
#' saddle point and thus not optimal.
#' 
#' 
#' @keywords regression multivariate bivariate riley meta-analysis
#' 
#' @author Thomas Debray <thomas.debray@gmail.com>
#'
#' @export
riley <- function(X, optimization = "Nelder-Mead", control = list(), pars, ...) UseMethod("riley")

#' @export
riley.default <- function(X, optimization = "Nelder-Mead", control = list(), pars, ...)
{
  
  pars.default <- list(level = 0.95)
  
  # Load default parameters
  if (!missing(pars)) {
    for (i in 1:length(pars)) {
      element <- ls(pars)[i]
      pars.default[[element]] <- pars[[element]]
    }
  }
  
  # Check parameter values
  if (pars.default$level < 0 | pars.default$level > 1) {
    stop ("Invalid value for significance level!")
  } 
  
  out <- NA

  
  if(missing(X)) stop("X is missing!")
  
  # Assess which type of meta-analysis is needed
  if (sum(c("Y1","vars1","Y2","vars2","TP","FN","FP","TN") %in% colnames(X))==8) {
    stop(paste("Too many variables specified in X! Please choose whether to perform a meta-analysis of effect sizes",
               "or a meta-analysis of diagnostic test accuracy data!"))
  }
  if (sum(c("Y1","vars1","Y2","vars2") %in% colnames(X))==4) {
    out <- rileyES(X, optimization = optimization, control=control, pars=pars.default, ...)
  } else if (sum(c("Y1","vars1","Y2","vars2") %in% colnames(X))>0) {
    stop ("Missing column(s) in X!")
  } else if (sum(c("TP","FN","FP","TN") %in% colnames(X))==4) {
    out <- rileyDA(X, optimization = optimization, control=control, pars=pars.default, ...)
  } else if (sum(c("TP","FN","FP","TN") %in% colnames(X))>0) {
    stop ("Missing column(s) in X!")
  } else {
    stop ("Provided data not supported, please verify column names in X!")
  }
  
	class(out) <- "riley"
	out
}

# X is the design matrix
negfullloglikRiley <- function(parsll, X, Y,vars)
{

  Beta = rbind(parsll[1], parsll[2]) #Beta vector
  psisq1 = parsll[3]**2 #ensure variance is positive
  psisq2 = parsll[4]**2 #ensure variance is positive
  rho = inv.logit(parsll[5])*2-1 #ensure correlation is in [-1,1], and values in that interval move symmetric from -1 to 0 and from 1 to 0
  k = 2 #2 endpoints
  n = dim(Y)[1]/2
  
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

# effect sizes data meta-analysis
rileyES <- function(X = NULL, Y1, Y2, vars1, vars2, optimization = "Nelder-Mead", control = list(), pars=list(level=0.95), ...)
{
	if(!is.null(X)){
		X <- as.data.frame(X)
		origdata <- X
		Y1 <- X$Y1
		Y2 <- X$Y2
		vars1 <- X$vars1
		vars2 <- X$vars2
	} else {
		origdata <- cbind(Y1,vars1,Y2,vars2)
		colnames(origdata) <- c("Y1","vars1","Y2","vars2")
	}
  
	numstudies <- length(Y1)
	nobs <- length(which(!is.na(Y1)))+length(which(!is.na(Y2)))
	
	if(nobs != numstudies*2){warning("There are missing observations in the data!")}
	
	df <- 5 #There are 5 parameters to estimate
	if(numstudies-df < 0){warning("There are very few primary studies!")}
	
	# Set up the design matrix
	Xdesign <- array(0,dim=c(numstudies*2,2))
	Xdesign[seq(from=1, to=(numstudies*2), by=2)[which(!is.na(Y1))],1] <- 1 #Indicate which variables are observed for Y1
	Xdesign[seq(from=2, to=(numstudies*2), by=2)[which(!is.na(Y2))],2] <- 1
	
	
	vars <- cbind(vars1, vars2)
	Y <- array(NA,dim=c((length(Y1)*2),1))
	Y[seq(from=1, to=(numstudies*2), by=2)] <- Y1
	Y[seq(from=2, to=(numstudies*2), by=2)] <- Y2
	
	#Calculate starting values for optim
	pars.start = c(0,0,0,0,0)
	if (numstudies >= 2) {
		sumlY1 <- metafor::rma(yi=Y1, vi=vars1, method="DL")
		sumlY2 <- metafor::rma(yi=Y2, vi=vars2, method="DL")
		pars.start = c(sumlY1$beta[1], sumlY2$beta[1], sqrt(vcov(sumlY1)[1,1]), sqrt(vcov(sumlY2)[1,1]), 0)
	}
	
	
	
	fit = optim(pars.start, negfullloglikRiley, X=Xdesign, Y=Y, vars=vars, method=optimization, hessian=T, control=control)
	
	if(fit$convergence != 0) { 
		if(fit$convergence == 1) warning ("Iteration limit had been reached.")
		else if (fit$convergence == 10) warning("Degeneracy of the Nelder-Mead simplex.")
		else if (fit$convergence == 51 | fit$convergence == 52) warning(fit$message)
    else warning("Unspecified convergence error in optim.")
	}
	
	beta1 = fit$par[1]
	beta2 = fit$par[2]
	psi1 = abs(fit$par[3])
	psi2 = abs(fit$par[4])
	rhoT = fit$par[5]
	coefficients = c(beta1,beta2,psi1,psi2,rhoT)
	names(coefficients) = c("beta1","beta2","psi1","psi2","rhoT")
	
	hessian = fit$hessian
	colnames(hessian) = c("beta1","beta2","psi1","psi2","rhoT")
	rownames(hessian) = c("beta1","beta2","psi1","psi2","rhoT")
	
	if (length(which(eigen(fit$hessian,symmetric=TRUE)$values<0))>0) warning("The Hessian contains negative eigenvalues!")
	
	iterations <- fit$iterations
	logLik <- -fit$value
	
	output <- list(coefficients = coefficients, hessian = hessian, df = df, numstudies = numstudies, nobs = nobs, logLik = logLik,
			   iterations = (iterations+1), call = match.call(), data = origdata, type="effect.size", level=pars$level)  
	return(output)
}

# Diagnostic test accuracy data meta-analysis
rileyDA <-
  function(X = NULL, TP, FN, FP, TN, correction = 0.5, 
           correction.control = "all", optimization = "Nelder-Mead", control = list(), pars=list(level=0.95), ...)
  {
      if(!is.null(X)){
        X <- as.data.frame(X)
        origdata <- newdata <- X
      } else {
        origdata <- newdata <- as.data.frame(cbind(TP,FN,FP,TN))
        colnames(origdata) <- c("TP","FN","FP","TN")
      }
      
	  ## The following corrections are copied from the "mada" package to facilitate comparison of results
      ## apply continuity correction to _all_ studies if one contains zero
      if(correction.control == "all"){if(any(origdata == 0)){newdata$TP <- origdata$TP + correction;
                                                             newdata$FN <- origdata$FN + correction;
                                                             newdata$FP <- origdata$FP + correction;
                                                             newdata$TN <- origdata$TN + correction}}
      if(correction.control == "single"){
        correction = ((((origdata$TP == 0)|(origdata$FN == 0))|(origdata$FP == 0))| (origdata$TN == 0))*correction
        newdata$TP <- correction + origdata$TP
        newdata$FN <- correction + origdata$FN
        newdata$FP <- correction + origdata$FP
        newdata$TN <- correction + origdata$TN
      }
      
      
      #Calculate sensitivities and specificities (original scale)
      number.of.pos <- newdata$TP + newdata$FN
      number.of.neg <- newdata$FP + newdata$TN
      sens <-newdata$TP/number.of.pos
      fpr <- newdata$FP/number.of.neg
      var.sens = sens*(1-sens)/number.of.pos
      var.fpr = fpr*(1-fpr)/number.of.neg
      
      logit.sens <- logit(sens)
      logit.fpr <- logit(fpr)
      var.logit.sens <- 1/(sens*(1-sens)*number.of.pos)
      var.logit.fpr <- 1/(fpr*(1-fpr)*number.of.neg)
      
	    #Apply ordinary bivariate meta-analysis on transformed data
      output = rileyES(X=NULL, Y1=logit.sens,Y2=logit.fpr,vars1=var.logit.sens,vars2=var.logit.fpr,optimization = optimization, control = control, ...)
      output$type = "test.accuracy"
      output$data = newdata
      output$correction = correction 
      output$correction.control = correction.control
      output$level <- pars$level
      
      return(output)
  }

#' @author Thomas Debray <thomas.debray@gmail.com>
#' @method print riley
#' @export
print.riley <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  if (length(which(eigen(x$hessian,symmetric=TRUE)$values<0))>0) cat("\nWarning: the Hessian matrix contains negative eigenvalues, parameter estimates are thus not optimally fitted!\n")
}

#' Prediction Interval
#' Calculates a prediction interval for the summary parameters of Riley's alternative model for bivariate random-effects 
#' meta-analysis. This interval predicts in what range future observations will fall given what has already been observed.
#' 
#' @param object  A \code{riley} object.
#' @param \dots Additional arguments (currently ignored)
#' 
#' @details Prediction intervals are based on Student's t-distribution with (numstudies - 5) degrees of freedom. The width of 
#' the interval is specified by the significance level chosen during meta-analysis.
#' 
#' @return Data frame containing prediction intervals with the summary estimates \code{beta1} and \code{beta2} 
#' (for effect size data), or with the mean sensitivity and false positive rate (for diagnostic test accuracy data).
#' 
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @method predict riley
#' @export
predict.riley <- function(object,  ...)
{
  alpha <- (1-object$level)/2
  
  predint		<- array(NA,dim=c(2,3))
  colnames(predint) <- c("Estimate", paste((alpha*100),"%"),paste(((1-alpha)*100),"%"))
  df 			<- object$df
  numstudies 		<- object$numstudies
  
  if (object$type=="test.accuracy")
  {
	rownames(predint) = c("Sens","FPR")
	if ((numstudies - df) > 0)
	{
		predint[1,] = inv.logit(c(coefficients(object)["beta1"],(qt(alpha,df=(numstudies-df))*sqrt((coefficients(object)["psi1"]**2)+diag(vcov(object))[1])+coefficients(object)["beta1"]),(qt((1-alpha),df=(numstudies-df))*sqrt((coefficients(object)["psi1"]**2)+diag(vcov(object))[1])+coefficients(object)["beta1"])))
		predint[2,] = inv.logit(c(coefficients(object)["beta2"],(qt(alpha,df=(numstudies-df))*sqrt((coefficients(object)["psi2"]**2)+diag(vcov(object))[2])+coefficients(object)["beta2"]),(qt((1-alpha),df=(numstudies-df))*sqrt((coefficients(object)["psi2"]**2)+diag(vcov(object))[2])+coefficients(object)["beta2"])))
	} else {
		predint[1,] = c(inv.logit(coefficients(object)["beta1"]),0,1)
		predint[2,] = c(inv.logit(coefficients(object)["beta2"]),0,1)
	}
  } else {
	rownames(predint) = c("beta1","beta2")
	predint[1,] = c(coefficients(object)["beta1"],(qt(alpha,df=(numstudies-df))*sqrt((coefficients(object)["psi1"]**2)+diag(vcov(object))[1])+coefficients(object)["beta1"]),(qt((1-alpha),df=(numstudies-df))*sqrt((coefficients(object)["psi1"]**2)+diag(vcov(object))[1])+coefficients(object)["beta1"]))
	predint[2,] = c(coefficients(object)["beta2"],(qt(alpha,df=(numstudies-df))*sqrt((coefficients(object)["psi2"]**2)+diag(vcov(object))[2])+coefficients(object)["beta2"]),(qt((1-alpha),df=(numstudies-df))*sqrt((coefficients(object)["psi2"]**2)+diag(vcov(object))[2])+coefficients(object)["beta2"]))
  }
  predint
}

#' @author Thomas Debray <thomas.debray@gmail.com>
#' @method summary riley
#' @export
summary.riley <- function(object, level = 0.95, ...)
{
	confints <- cbind(object$coefficients, confint(object,level=level))
	colnames(confints)[1] <- "Estimate"
	
	if (object$type=="test.accuracy") {
		confints <- rbind(inv.logit(confints[1:2,]),confints)
		rownames(confints)[1:2] <-  c("Sens", "FPR") 
	} 
	
	#Transform last parameter back to rho
	confints["rhoT",] =  inv.logit(confints["rhoT",])*2-1
	rownames(confints)[which(rownames(confints)=="rhoT")] = "rho"
	
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

#' @author Thomas Debray <thomas.debray@gmail.com>
#' @method vcov riley
#' @export
vcov.riley <- function(object, ...) {
  if (length(which(eigen(object$hessian,symmetric=TRUE)$values<0))>0) warning("The Hessian contains negative eigenvalues!")
 
  # It is known that 'optim' has problems.  Perhaps the simplest thing to do is to call 'optim' with each of 
  # the 'methods' in sequence, using the 'optim' found by each 'method' as the starting value for the next.  
  # When I do this, I often skip 'SANN', because it typically takes so much more time than the other methods.  
  # However, if there might be multiple local minima, then SANN may be the best way to find a global minimum, 
  # though you may want to call 'optim' again with another method, starting from optimal solution returned by 'SANN'. 

  Sigma = solve(object$hessian)
  Sigma
}

#' Print the log-likelihood
#' 
#' This function provides the (restricted) log-likelihood of a fitted model.
#' 
#' @param  object A \code{riley} object, representing a fitted alternative model for bivariate random-effects 
#' meta-analysis when the within-study correlations are unknown.
#' @param \dots Additional arguments to be passed on to other functions, currently ignored.
#' @return Returns an object of class \code{logLik}. This is the (restricted) log-likelihood of the model represented 
#' by \code{object} evaluated at the estimated coefficients. It contains at least one attribute, 
#' "\code{df}" (degrees of freedom), giving the number of (estimated) parameters in the model.
#' 
#' @references Riley RD, Thompson JR, Abrams KR. An alternative model for bivariate random-effects meta-analysis when 
#' the within-study correlations are unknown. \emph{Biostatistics} 2008; \bold{9}: 172--186.
#' 
#' @examples 
#' data(Daniels)
#' fit <- riley(Daniels,control=list(maxit=10000))
#' logLik(fit)
#' 
#' @keywords likelihood riley bivariate meta-analysis
#' 
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @method logLik riley
#' @export
logLik.riley <- function(object, ...) {
	val 				      <- object$logLik
	attr(val, "nobs") <- object$nobs
	attr(val, "df") 	<- object$df
	class(val) 			  <- "logLik"
	return(val)
}


plot.riley <- function(x, plotsumm = TRUE, plotnumerics = TRUE, level = 0.95, main="",
                       ylim = c(0,1), xlim = c(0,1), pch = 1, lty = 1, lwd = 1, cex.numerics=0.45,
                       add=FALSE, ...)
{
	alpha = (1-level)/2
	
	if (x$type=="test.accuracy") {
		xlab = "1-Specificity"
		ylab = "Sensitivity"
		
		FP <- x$data$FP
		negatives <- FP + x$data$TN
		FPR <- FP/negatives
		mu = x$coefficients[c("beta2","beta1")]
		Sigma = vcov(x)[c("beta2","beta1"),c("beta2","beta1")] 
		mu.ellipse <- ellipse(Sigma, centre = mu, level = level) 
		summary1 = inv.logit(mu[1])
		summary2 = inv.logit(mu[2])
		ellipse1 = inv.logit(mu.ellipse[,1])
		ellipse2 = inv.logit(mu.ellipse[,2])
	} else {
		plotnumerics = FALSE
		xlab = "Y1"
		ylab = "Y2"
		mu = x$coefficients[c("beta1","beta2")]
		Sigma = vcov(x)[c(1,2),c(1,2)] 
		mu.ellipse <- ellipse(Sigma, centre = mu, level = level) 
		summary1 = mu[1]
		summary2 = mu[2]
		ellipse1 = mu.ellipse[,1]
		ellipse2 = mu.ellipse[,2]
	}
	
	if (!add) plot(-500,-500, type = "l", xlim = xlim, ylim = ylim, xlab=xlab,ylab=ylab,main=main, ...)
	#if (!add) NextMethod("plot")	
	polygon(ellipse1,ellipse2,lty=lty, lwd=lwd)
	if(plotsumm) points(summary1,summary2,pch=pch) # add the point estimate of the mean
	
	if(plotnumerics) {
		ci = summary(x,level=level)[2]$confints
		text(0.8,0.15,labels="Estimate",pos=2,cex=cex.numerics)
		text(0.9,0.15,labels=paste((alpha*100),"% CI",sep=""),pos=2,cex=cex.numerics)
		text(1.0,0.15,labels=paste(((1-alpha)*100),"% CI",sep=""),pos=2,cex=cex.numerics)
		text(0.5,0.10,labels= "Sensitivity",pos=4, cex=cex.numerics)
		text(0.8,0.10,labels=paste("",formatC(round( ci["Sens",1],2),2,format="f",flag="0")),pos=2, cex=cex.numerics)
		text(0.9,0.10,labels=paste("",formatC(round( ci["Sens",2],2),2,format="f",flag="0")),pos=2, cex=cex.numerics)
		text(1.0,0.10,labels=paste("",formatC(round( ci["Sens",3],2),2,format="f",flag="0")),pos=2, cex=cex.numerics)
		text(0.5,0.05,labels= "Specificity",pos=4, cex=cex.numerics)
		text(0.8,0.05,labels=paste("",formatC(round( 1-ci["FPR",1],2),2,format="f",flag="0")),pos=2, cex=cex.numerics)
		text(0.9,0.05,labels=paste("",formatC(round( 1-ci["FPR",3],2),2,format="f",flag="0")),pos=2, cex=cex.numerics)
		text(1.0,0.05,labels=paste("",formatC(round( 1-ci["FPR",2],2),2,format="f",flag="0")),pos=2, cex=cex.numerics)
	}
}
