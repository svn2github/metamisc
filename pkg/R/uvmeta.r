# Meta-analysis is a statistical technique by which information from
# independent studies is assimilated. This function allows to perform
# fixed-effects and random-effects meta-analysis.
# r    : Vector of the effect sizes
# vars : Vector of the effect variances
################################################################################
# Author  : Thomas Debray
# Version : 10 May 2011
################################################################################
# Example
# example.r = c(0.10,0.30,0.35,0.65,0.45,0.15)
# example.var = c(0.03,0.03,0.05,0.01,0.05,0.02)
# macc(example.r,example.var)
################################################################################
# References
# * Tania Huedo-Medina; Julio Sanchez-Meca; Fulgencio Marin-Martinez;
#   Juan Botella,"Assessing heterogeneity in meta-analysis: Q statistci
#   or I2 index?" (2006). CHIP Documents . Paper 19
#   http://digitalcommons.uconn.edu/chip_docs/19
# * Michael Borenstein; Larry Hedges; Hahhah Rothstein, "Meta-Analysis .
#   Fixed effect vs. random effects" (2007) . www.Meta-Analysis.com
################################################################################

uvmeta <- function(x, ...) UseMethod("uvmeta")

uvmetaMOM <- function(r,vars) {
    ############################################################################
    # FIXED EFFECTS MODEL
    ############################################################################
    w = 1/vars

    #Combined effect
    weighted_Tbar = sum(r*w)/sum(w)

    # Variance of the combined effect
    var_T = 1/sum(w)

    # Standard error combined effect
    se_T = sqrt(var_T)

    # The Z-value
    z_T = weighted_Tbar/se_T

    ############################################################################
    # RANDOM EFFECTS MODEL
    ############################################################################
    # Q-statistic
    Q = sum(w*(r-weighted_Tbar)**2)
    I_sq = 0
    
    # Higgins and Thompson (2002) have also developed a confidence interval for
    # I2. The interval is formulated by calculating another of their proposed
    # measures of heterogeneity, the H2 index obtained by
    # (Higgins & Thompson, 2002, p. 1545, eq. 6) also known as Birge?s
    # ratio (Birge, 1932)
    H_sq = 0    #H2
    se_lnH = 0 #SE(ln(H2))

    # Degrees of freedom
    dfr = length(r)-1

    # Between-study variance
    if (Q > dfr) {
        re_C =  sum(w) - sum(w**2)/sum(w)
        between_study_var = (Q - dfr)/re_C
        I_sq = (Q-dfr)/Q
        H_sq = Q/dfr
        se_lnH = (log(Q)-log(dfr))/(2*(sqrt(2*Q)-sqrt((2*length(r))-3)))
    } else {
        between_study_var = 0
        se_lnH = sqrt((1/(2*(length(r)-2)))*(1-(1/(3*((length(r)-2)**2)))))
    }

    # Within-study plus between-study variance
    re_v = vars + between_study_var

    # Updated weights
    re_w = 1/re_v

    # Combined effect
    re_weighted_Tbar =  sum(r*re_w)/sum(re_w)

    # Variance of the combined effect
    re_var_T  = 1/sum(re_w)

    # Standard error of combined effect
    re_se_T = sqrt(re_var_T)

    # The Z-value
    re_z_T = re_weighted_Tbar/re_se_T
    
    ############################################################################
    # Fixed effect statistics  (2-tailed)
    ############################################################################
    fe_p1t = pnorm(weighted_Tbar/sqrt(var_T),lower=F)
    fe_p2t = 2*pnorm(weighted_Tbar/sqrt(var_T),lower=F)
    fe.lowerconf =  weighted_Tbar + qnorm(0.025)*sqrt(var_T)
    fe.upperconf =  weighted_Tbar + qnorm(0.975)*sqrt(var_T)
    fixef.results = list(mean=weighted_Tbar,var=var_T,p2t=fe_p2t,
        lowerconf=fe.lowerconf,upperconf=fe.upperconf)

    ############################################################################
    # Random effect statistics (2-tailed)
    ############################################################################
    re_p1t = pnorm(re_weighted_Tbar/sqrt(re_var_T),lower=F)
    re_p2t = 2*pnorm(re_weighted_Tbar/sqrt(re_var_T),lower=F)
    re.lowerconf = re_weighted_Tbar + qnorm(0.025)*sqrt(re_var_T)
    re.upperconf = re_weighted_Tbar + qnorm(0.975)*sqrt(re_var_T)
    ranef.results = list(mean=re_weighted_Tbar,var=re_var_T,p2t=re_p2t,
        lowerconf=re.lowerconf,upperconf=re.upperconf)
    
    ############################################################################
    # H statistics
    ############################################################################
    lnH = log(sqrt(H_sq))
    lnH_p2t = 2*pnorm(lnH/se_lnH,lower=F) # P-value of ln(H)
    lnH.lowerconf = lnH + qnorm(0.025)*se_lnH
    lnH.upperconf = lnH + qnorm(0.975)*se_lnH
    H2.lowerconf = (exp(lnH.lowerconf))**2 #95% CI
    H2.upperconf = (exp(lnH.upperconf))**2
    
    H2.results =  list(H2=H_sq,lnH=lnH,lnH_p2t=lnH_p2t,lnH.lowerconf=lnH.lowerconf,
        lnH.upperconf=lnH.upperconf,H2.lowerconf=H2.lowerconf,H2.upperconf=H2.upperconf)

    ############################################################################
    # I square statistics
    ############################################################################
    I2.lowerconf = (H2.lowerconf-1)/H2.lowerconf
    I2.upperconf = (H2.upperconf-1)/H2.upperconf
    I2.results = list(I2=I_sq,I2.lowerconf=I2.lowerconf,I2.upperconf=I2.upperconf)
    
    ############################################################################
    # Q statistics
    ############################################################################
    # The Q statistic has a chi-square distribution with k ? 1 degrees of
    # freedom, k being the number of studies. Thus, Q values higher than the
    # critical point for a given significance level alfa enable us to reject the
    # null hypothesis and conclude that there is statistically significant
    # between-study variation
    Q.critical = qchisq(0.95,df=(length(r)-1))
    Q_p = 1-pchisq(Q,df=(length(r)-1))
    Q.results = list(Q=Q,critical=Q.critical,p.value=Q_p)
    
    ############################################################################
    # Output
    ############################################################################
    out <- list(fixef=fixef.results,ranef=ranef.results,Q=Q.results,H2=H2.results,
        I2=I2.results,tau_sq=between_study_var)
    return (out)
}

uvmeta.default <- function(r,v, method="MOM", ...)
{
    x = as.vector(r)
    y = as.vector(v)
    est <- NA
    if (method == "MOM") {
	est <- uvmetaMOM(x, y)
	est$call <- match.call()
        class(est) <- "uvmeta"
    }
    est
}

print.uvmeta <- function(x, ...)
{
    cat("Call:\n")
    print(x$call)
    cat("\nCoefficients:\n\n")
    fe = array(NA,dim=c(2,3))
    fe[1,] = cbind(round(x$fixef$mean,5),round(sqrt(x$fixef$var),5),x$fixef$p2t)
    fe[2,] = cbind(round(x$ranef$mean,5),round(sqrt(x$ranef$var),5),x$ranef$p2t)
    colnames(fe) = c("Estimate","StdErr","p.value (2-tailed)")
    rownames(fe) = c("Fixed Effects","Random Effects")
    print(fe,quote=F)
    cat(paste("\nTau squared: \t\t",round(x$tau_sq,5),sep=""))
    cat(paste("\nCochran's Q statistic: \t",round(x$Q$Q,5)," (p.value: ",round(x$Q$p.value,5),")",sep=""))
    cat(paste("\nH-square index: \t", round(x$H2$H2,5)," (95% CI: ",round(x$H2$H2.lowerconf,5),"; ",round(x$H2$H2.upperconf,5),")",sep=""))
    cat(paste("\nI-square index: \t", round(x$I2$I2*100,5),"% (95% CI: ",round(x$I2$I2.lowerconf*100,5),"% ; ",round(x$I2$I2.upperconf*100,5),"%)\n",sep=""))
}

summary.uvmeta <- function(object, ...)
{
    se <- c(sqrt(object$fixef$var),sqrt(object$ranef$var),NA,NA,NA)
    lower.conf = c(object$fixef$lowerconf,object$ranef$lowerconf,NA,object$H2$H2.lowerconf,object$I2$I2.lowerconf)
    upper.conf = c(object$fixef$upperconf,object$ranef$upperconf,NA,object$H2$H2.upperconf,object$I2$I2.upperconf)
    
    TAB = cbind(Estimate = c(object$fixef$mean,object$ranef$mean,object$Q$Q,object$H2$H2,object$I2$I2),StdErr=se,lower.conf=lower.conf,upper.conf=upper.conf)
    rownames(TAB) = c("Fixed Effects","Random Effects","Cochran Q","H-square index","I-square index")
    res = list(call=object$call,estimates=TAB)
    class(res) = "summary.macc"
    res
}



