context("Testing frequentist valmeta functionalities (O:E ratio)")
library(lme4)

data(EuroSCORE)
skip_on_cran()

test_that("Random effects meta-analysis of total O:E ratio works (Poisson distribution)", {
  
  ds <- EuroSCORE
  ds$Study <- c(1:dim(ds)[1])

  fit1.lme4 <- glmer(n.events ~ 1 | Study, offset = log(e.events), 
                     family = poisson(link = "log"), data = ds)
  
  fit1.valmeta <- with(EuroSCORE, valmeta(measure="OE", 
                                          O=n.events, 
                                          E=e.events, 
                                          method="ML", pars=list(model.oe="poisson/log")))
  
  expect_identical(fit1.valmeta$est, exp(fixef(fit1.lme4)))
})


