library(ars)
## Consistency

test_that("Consistency", {
  h <- function(x)
  {
    log(dnorm(x,mean=10,sd=2))
  }
  n<-10000
  domain <- c(-Inf,Inf)
  samples <- ars(h,n,domain)
  
  expect_that( samples, is_a("numeric") )
  expect_that( length(samples), equals(n) )
})
