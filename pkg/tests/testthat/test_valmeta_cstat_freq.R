context("Testing frequentist valmeta functionalities (c-statistic)")
library(metafor)

data(EuroSCORE)
skip_on_cran()


test_that("Fixed effect meta-analysis of c-statistic works", {
  
  ## Prepare data
  dat.cstat <- ccalc(cstat=c.index, cstat.se=se.c.index, 
                     cstat.cilb=c.index.95CIl, cstat.ciub=c.index.95CIu, 
                     N=n, O=n.events, data=EuroSCORE, 
                     g="log(cstat/(1-cstat))")
  
  fit1.metafor <- rma(yi=theta, sei=theta.se, method="FE", data=dat.cstat, test="z")
  
  fit1.valmeta <- with(EuroSCORE, valmeta(cstat=c.index, cstat.se=se.c.index, 
                                          cstat.95CI=cbind(c.index.95CIl,c.index.95CIu), 
                                          N=n, O=n.events, slab=Study, method="FE", test="z"))
  
  expect_equal(fit1.valmeta$est, inv.logit(as.numeric(fit1.metafor$beta)))
})
