context("metapred and its S3 methods.")

### TODO
# predFUN.
# Tests for options of predict.metapred

### Some stuff necessary for testing
# The data
set.seed(8092017)
n <- 100
n.cov <- 3
td <- data.frame(matrix(rbinom(n * (n.cov + 1), 1, .5), ncol = n.cov + 1, nrow = n))
td.ig <- td + 1 # For inverse gaussian and Gamma.

# Arguments
f <- X1 ~ X2 + X3
gl <- glm(f, family = binomial, data = td)
st.i <- td[["X4"]]
st.u <- sort(unique(st.i))
folds <- l1o(st.u)

# Ideally these would be separated into 3 tests, but they use each other's objects, which are cleaned up by test_that.
test_that("Stratified models can be estimated and MA.", {
  # Stratified estimation
  expect_is(stratum.fit <- mp.stratum.fit(gl), "mp.stratum.fit")
  expect_is(stratified.fit <- mp.stratified.fit(formula = f, data = td, 
                                                           st.i = st.i, st.u = st.u, estFUN = glm, 
                                      family = binomial)  , "mp.stratified.fit")
  
  # MA 
  # Note: this essentially yields the global model
  expect_is(meta.fit <- mp.meta.fit(stratified.fit, urma), "mp.meta.fit")
  expect_length(coef(meta.fit), 3)
  expect_is(cv.meta.fit <- mp.cv.meta.fit(stratified.fit = stratified.fit, folds = folds, metaFUN = urma), "mp.cv.meta.fit")
  expect_equal(dim(coef(cv.meta.fit)), c(length(unique(td[["X4"]])), ncol(td) - 1) )
  
  # Recal of MA
  expect_is(recal.meta.fit <- mp.recal.meta.fit(meta.fit = meta.fit, formula = f, newdata = td, estFUN = glm), "mp.meta.fit")
  expect_gt(coef(recal.meta.fit)[1], coef(meta.fit)[1])
  expect_identical(coef(recal.meta.fit)[-1], coef(meta.fit)[-1])
})

test_that("Stratified models can be cross-validated", {
  # CV: development
  expect_is(cv.dev <- mp.cv.dev(formula = f, data = td, st.i = st.i, st.u = st.u, folds = folds, 
                                estFUN = glm, metaFUN = urma, meta.method = "DL", family = binomial), "mp.cv.dev")
  expect_equal(family(cv.dev), binomial())
  expect_equal(getPredictMethod(fit = cv.dev, two.stage = TRUE), metamisc:::predictGLM)
  
  # CV: recalibration
  # Note: double recalibration sadly removes the previous original coefs and should not be done.
  expect_is(cv.recal  <- mp.cv.recal(cv.dev = cv.dev, newdata = td, folds = folds, estFUN = glm), "mp.cv.dev")
  expect_is(cv.recal2 <- mp.cv.recal(cv.recal, td, folds, glm), "mp.cv.dev")
  
  # CV: validation
  expect_is(cv.val <- mp.cv.val(cv.dev = cv.dev, data = td, st.i = st.i, folds = folds), "mp.cv.dev")
  expect_is(cv.val, "mp.cv.val")
  
  # CV: validation with recalibration
  # Note: with recalibration, gen and perf should always appear to be smaller (=better)
  expect_is(cv.val.recal <- mp.cv.val(cv.dev = cv.dev, data = td, st.i = st.i, folds = folds, 
                                      recal.int = TRUE, estFUN = glm), "mp.cv.dev")
  expect_gt(cv.val$gen, cv.val.recal$gen)
  expect_true(all(cv.val$perf$perf > cv.val.recal$perf$perf))
  
  # Whole sequence
  expect_is(cv <- mp.cv(formula = f, data = td, st.i = st.i, st.u = st.u, folds = folds, family = binomial), "mp.cv")
  expect_true(all(class(cv) == c("mp.cv", "mp.cv.val", "mp.cv.dev")))
  
  # Global model
  expect_is(global <- mp.global(cv.val, urma), "mp.global")
  expect_is(global <- mp.global(cv, urma), "mp.global")
  expect_is(global <- mp.global(cv.val.recal, urma), "mp.global")
})

test_that("A stepwise stratified model can be fitted", {
  # No stepwise
  expect_is(step0 <- mp.step(formula = f, data = td, remaining.changes = c(""), st.i = st.i, 
                             st.u = st.u, folds = folds, family = binomial), "mp.step")
  expect_length(step0$cv, 1)
  expect_is(mp.step.get.best(step0), "mp.cv")
  
  # Main effects
  change.main <- c("X2", "X3")
  expect_is(step <- mp.step(formula = f, data = td, remaining.changes = change.main, st.i = st.i, 
                            st.u = st.u, folds = folds, family = binomial), "mp.step")
  expect_is(mp.step.get.best(step), "mp.cv")
  
  # Interaction effect
  change.interaction <- c("X2:X3")
  expect_is(step2 <- mp.step(formula = f, data = td, remaining.changes = change.interaction, st.i = st.i, 
                             st.u = st.u, folds = folds, family = binomial), "mp.step")
  expect_is(mp.step.get.best(step2), "mp.cv")
  
  # Entire fit
  expect_is(fit <- mp.fit(formula = f, data = td, remaining.changes = change.main, st.i = st.i, 
                st.u = st.u, folds = folds, family = binomial, max.steps = 3), "mp.fit")
  expect_equal(fit$best.step, "s1") # for this data and seed
})

test_that("metapred produces a model.", {
  # Stepwise, default
  expect_is(mp <- metapred(data = td, strata = "X4", scope = f, formula = X1 ~ 1, family = binomial), "metapred")
  expect_identical(mp$step.count, 2) # and: stop.reason == no improvement was possible, but that may change.
  
  # Stepwise, stop because of step.count
  expect_is(mp <- metapred(data = td, strata = "X4", scope = f, formula = X1 ~ 1, family = binomial, max.steps = 1), "metapred")
  expect_identical(mp$step.count, 1)
  
  # No stepwise
  expect_is(mp <- metapred(data = td, strata = "X4", scope = f, formula = f, family = binomial), "metapred")
  expect_identical(mp$step.count, 0)
  
  # Possible frequent user mistake: wrong strata variable or it is missing from data.
  # The error below is the intended message. Somehow it is not shown in this test.
  expect_error(mp <- metapred(data = td, strata = "X5", scope = f, formula = f, family = binomial) ) #,
               # "Error in `[.data.frame`(data, , strata) : undefined columns selected")
})

test_that("metapred can handle different distributions.", {
  expect_true(is.list(mp <- metapred(data = td, strata = "X4", family = binomial))) # binomial
  expect_true(inherits(mp, "metapred"))
  expect_true(is.list(mp$stepwise))
  expect_true(is.list(mp$FUN))
  expect_true(is.call(mp$call))

  expect_true(is.list(mp <- metapred(data = td, strata = "X4", family = binomial(link = "log")))) # binomial, loglink
  expect_true(inherits(mp, "metapred"))
  expect_true(is.list(mp$stepwise))
  expect_true(is.list(mp$FUN))
  expect_true(is.call(mp$call))

  expect_true(is.list(mp <- metapred(data = td, strata = "X4"))) # gaussian
  expect_true(inherits(mp, "metapred"))
  expect_true(is.list(mp$stepwise))
  expect_true(is.list(mp$FUN))
  expect_true(is.call(mp$call))

  expect_true(is.list(mp <- metapred(data = td.ig, strata = "X4", family = Gamma))) # Gamma
  expect_true(inherits(mp, "metapred"))
  expect_true(is.list(mp$stepwise))
  expect_true(is.list(mp$FUN))
  expect_true(is.call(mp$call))

  expect_true(is.list(mp <- metapred(data = td.ig, strata = "X4", family = inverse.gaussian))) # inverse.gaussian
  expect_true(inherits(mp, "metapred"))
  expect_true(is.list(mp$stepwise))
  expect_true(is.list(mp$FUN))
  expect_true(is.call(mp$call))

  expect_true(is.list(mp <- metapred(data = td, strata = "X4", family = poisson))) # poisson
  expect_true(inherits(mp, "metapred"))
  expect_true(is.list(mp$stepwise))
  expect_true(is.list(mp$FUN))
  expect_true(is.call(mp$call))

  expect_true(is.list(mp <- metapred(data = td, strata = "X4", family = quasi))) # quasi
  expect_true(inherits(mp, "metapred"))
  expect_true(is.list(mp$stepwise))
  expect_true(is.list(mp$FUN))
  expect_true(is.call(mp$call))

  expect_true(is.list(mp <- metapred(data = td, strata = "X4", family = quasibinomial))) # quasibinomial
  expect_true(inherits(mp, "metapred"))
  expect_true(is.list(mp$stepwise))
  expect_true(is.list(mp$FUN))
  expect_true(is.call(mp$call))

  expect_true(is.list(mp <- metapred(data = td, strata = "X4", family = quasipoisson))) # quasipoisson
  expect_true(inherits(mp, "metapred"))
  expect_true(is.list(mp$stepwise))
  expect_true(is.list(mp$FUN))
  expect_true(is.call(mp$call))
})

test_that("metapred's stepwise is WAD.", {
  expect_true(is.list(mp <- metamisc:::metapred(data = td, strata = "X4" ))) 
  expect_true(inherits(mp, "metapred"))
  expect_true(is.list(mp$stepwise))
  expect_true(is.list(mp$FUN))
  expect_true(is.call(mp$call))
  expect_length(coef(mp), 2) # One is selected due to random fluctuation.
  
  set.seed(324234)
  td.none <- data.frame(matrix(rbinom(n * (n.cov + 1), 1, .5), ncol = n.cov + 1, nrow = n))
  expect_true(is.list(mp <- metamisc:::metapred(data = td.none, strata = "X4" ))) 
  expect_true(inherits(mp, "metapred"))
  expect_true(is.list(mp$stepwise))
  expect_true(is.list(mp$FUN))
  expect_true(is.call(mp$call))
  expect_length(coef(mp), 1) # None are selected, because the data is pure noise.
  
  td.all <- data.frame(matrix(rbinom(n * (n.cov + 1), 1, .5), ncol = n.cov + 1, nrow = n))
  td.all[ , 1] <- rowSums(td.all)
  expect_true(is.list(mp <- metamisc:::metapred(data = td.all, strata = "X4" ))) 
  expect_true(inherits(mp, "metapred"))
  expect_true(is.list(mp$stepwise))
  expect_true(is.list(mp$FUN))
  expect_true(is.call(mp$call))
  expect_length(coef(mp), 3) # All are selected, as predictors are good predictors.
  
  expect_true(is.list(mp <- metamisc:::metapred(data = td.none, strata = "X4", 
                                                formula = f, scope = f )))
  expect_true(inherits(mp, "metapred"))
  expect_true(is.list(mp$stepwise))
  expect_true(is.list(mp$FUN))
  expect_true(is.call(mp$call))
  expect_length(coef(mp), 3)  # All noise predictors are selected, because stepwise = F.
})

test_that("coef.metapred gets the coefficients", {
  mp <- metapred(data = td, strata = "X4", family = binomial)
  expect_true(is.numeric(coef(mp)))

  gl.b <- glm(formula = X1 ~ X2 + X3, data = td, family  = binomial)
  expect_equal(length(coef(mp)), n.cov - 1)
})

# This one can be a little annoying
# test_that("print.metapred prints a metapred object", {
#   mp <- metapred(data = td, strata = "X4", family = binomial)
#   cat("\n")
#   print(mp)
# })

test_that("metapred.predict predicts.", {
  mp <- metapred(data = td, strata = "X4", family = binomial)
  p <- predict(mp, newdata = td)
  expect_true(is.numeric(p))
  expect_true(all(p <= 1))
  expect_true(all(p >= 0))
})

test_that("metapred.family and gets the family.", {
  expect_true(is.list(mp <- metapred(data = td, strata = "X4", family = binomial(link = "log"))))
  gl <- glm(formula = X1 ~ X2 + X3, data = td, family = binomial(link = "log"))
  expect_equal(family(gl), family(mp))
})

test_that("metapred.formula gets the formula (test-dependent).", {# formula of glm is specific for this data set!
  expect_true(is.list(mp <- metapred(data = td, strata = "X4", family = binomial(link = "log"))))
  gl <- glm(formula = X1 ~ X3, data = td, family = binomial(link = "log"))
  expect_equal(formula(gl), formula(mp)) 
})

test_that("metapred.subset gets the right object.", {
  mp <- metapred(data = td, strata = "X4", family = binomial)
  expect_is(subset(mp), "mp.cv") # needs better test: test wether it is best.
  expect_is(subset(mp, select = "global"), "mp.global")
})

# This is to prevent a previous bug from reappearing
# If the order of formla does not match that of data, it needs to be reordered internally
# This can give issues when centering.
test_that("metapred uses formula correctly.", {
  expect_is(mp <- metapred(data = td, strata = "X2", formula = X3 ~ X1 + X4, family = binomial), "metapred")
})