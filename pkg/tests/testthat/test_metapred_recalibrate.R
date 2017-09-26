### Some stuff necessary for testing
set.seed(8092017)
n <- 100
n.cov <- 3
d  <- data.frame(matrix(rbinom(n * (n.cov + 1), 1, .7), ncol = n.cov + 1, nrow = n))
td <- data.frame(matrix(rbinom(n * (n.cov + 1), 1, .5), ncol = n.cov + 1, nrow = n))
d3 <- data.frame(matrix(rbinom(n * (n.cov + 1), 1, .5), ncol = n.cov + 1, nrow = n))



# These tests are rather weak. They mostly just test whether the code runs.

test_that("computeInt recalibrates", {
  g <- glm(d, family = binomial)
  expect_true(is.numeric(r.int <- metamisc:::computeInt(g, newdata = td, estFUN = glm, family = binomial)))
  expect_true(r.int[[1]] <= coef(g)[1])
})

test_that("recalibrate recalibrates", {
  ### basically recalibrate.glm
  g <- glm(d, family = binomial)
  expect_true(is.list(r <- recalibrate(g, td, family = binomial)))
  expect_true(inherits(r, "glm"))
  r.int <- metamisc:::computeInt(r, newdata = td, estFUN = glm, family = binomial)
  expect_equal(r.int[[1]], coef(r)[[1]]) # ComputeInt after recalibrate should have no added effect.
  
  ### basically recalibrate.metapred
  mp <- metapred(td, strata = "X4", recal.int = FALSE)
  mp.r <- recalibrate(mp, newdata = td) # ignores clustering. Hence, intercept changes.
  expect_true(is.list(mp.r))
  expect_true(inherits(mp.r, "metapred"))
})

test_that("metapred with recal.int = T works.", {
  mp.t <- metapred(d, strata = "X4", recal.int = TRUE)
  expect_true(is.list(mp.t$stepwise$coefficients.recal))
  
  mp.f <- metapred(d, strata = "X4", recal.int = FALSE)
  expect_true(is.null(mp.f$stepwise$coefficients.recal))
})


test_that("predict.metapred with recal.int = T works.", {
  mp.f <- metapred(d, strata = "X4", recal.int = FALSE)
  expect_true(is.numeric(p <- predict(object = mp.f, newdata = td, recal.int = TRUE)))
  expect_true(all(p <= 1))
  expect_true(all(p >= 0))
  
  mp.f <- metapred(d3, strata = "X4", recal.int = FALSE)
  expect_true(is.numeric(p <- predict(object = mp.f, newdata = td, recal.int = TRUE)))
  expect_true(all(p <= 1))
  expect_true(all(p >= 0))
  
  expect_true(is.numeric(p2 <- predict(object = mp.f, newdata = td, recal.int = FALSE)))
  expect_true(all(p2 <= 1))
  expect_true(all(p2 >= 0))
  
  expect_false(all(p == p2))
  expect_true(all(p != p2)) # With this seed.
})
