context("formula functions for metapred")

f0    <- formula(y ~ 1)
fa    <- formula(y ~ a)
fab   <- formula(y ~ a + b)
fabc  <- formula(y ~ a + b + c)
fab0  <- formula(y ~ a + b - 1)
faxb  <- formula(y ~ a * b)
fab2  <- formula(y ~ a + I(b^2))
fabb2 <- formula(y ~ a + I(b^2) + b) # Reordered, for testing purposes
fapb  <- formula(y ~ I(a+b))

# utility function for obtaining terms (predictors, or interactions etc) of formula
test_that("terms can be obtained", {
  expect_true(all(f2tl(fab0) == c("a", "b")))
  expect_false(any(f2tl(fab) == c("b", "a"))) # sadly gives false, 
  expect_true(all(c("b", "a") %in% f2tl(fab))) # that's why %in% is used internally.
})


# Add predictors
test_that("terms can be added", {
  expect_equal(updateFormula(f0, "a"), fa)
  expect_equal(updateFormula(fa, "b"), fab)
  expect_equal(updateFormula(fa, c("b", "c")), fabc)
  expect_equal(updateFormula(fab0, 1), fab)
  expect_equal(updateFormula(fab2, "b"), fabb2)
  expect_equal(updateFormula(f0, "I(a+b)"), fapb)
})

# Remove predictors
test_that("terms can be removed", {
  expect_equal(updateFormula(fab, c("a", "b")), f0)
  expect_equal(updateFormula(fab, "b"), fa)
  expect_equal(updateFormula(fa, "a"), f0)
  expect_equal(updateFormula(fabb2, "b"), fab2)
  expect_equal(updateFormula(fabb2, "I(b^2)"), fab)
  #expect_equal(updateFormula(fab0, c("a", "b")))  # a no intercept-model without predictors gets a bugged formula :(
})

# Error flagging
test_that("error is raised when terms are added and removed", {
  expect_error(updateFormula(fab, c("z", "b"))) # error, as intended
})

















