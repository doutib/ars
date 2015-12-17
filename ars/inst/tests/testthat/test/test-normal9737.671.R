
source("test-ks.R")

## Normal Distribution
distibution.test("N(9737.67,1)",c(-Inf,Inf),dnorm,pnorm,mean=9737.67,sd=1)
