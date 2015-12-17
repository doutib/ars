library(ars)
context("abscissae test")

test_that ("invalid domain", {
  expect_warning(abscissae(function(x) log(dnorm(x, mean = 0, sd = 1)), c(Inf, -Inf)))
})
