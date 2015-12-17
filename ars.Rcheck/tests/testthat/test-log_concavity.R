library(ars)
context("log concavity check")

h <- function (x) {
  dcauchy(x, location = 0, scale = 1, log = T)
}

u <- envelop(h,c(-10, -5, 0, 5, 10), c(-Inf, Inf))
l <- squeezing(h, c(-10, -5, 0, 5, 10))

test_that ("detect non-log-concave function", {
  expect_false(log_concavity(u, l))
})

h <- function (x) {
  df(x, df1=1, df2=1, log = T)
}

u <- envelop(h,c(0, 1, 2, 3, 4, 5, 6), c(0, Inf))
l <- squeezing(h, c(0, 1, 2, 3, 4, 5, 6))

test_that ("detect non-log-concave function", {
  expect_false(log_concavity(u, l))
})
