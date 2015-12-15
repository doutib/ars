
## Normal Distribution

test_that("N(9737.67,0.1) Unbounded domain", {
  domain <- c(-Inf,Inf)
  n<-10000
  g1 <- dnorm
  h1 <- function(x)
  {
    log(dnorm(x,mean=9737.67,sd=0.1))
  }
  # Run
  samples <- ars(h1,n,domain)
  ks_test=ks.test(samples,"pnorm",mean=9737.67,sd=0.1)
  
  # Test
  expect_that( ks_test$p.value>=0.05, is_true() )
  
  # Plot
  hist(samples,freq=F,breaks=100)
  lines(sort(samples),dnorm(sort(samples),mean=9737.67,sd=0.1))
  ks_test
})

