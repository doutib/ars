
source("test-ks.R")

## Normal Distribution
distibution.test("N(2,0.1)",c(-Inf,Inf),dnorm,pnorm,mean=2,sd=0.1)
