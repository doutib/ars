library(ars)
context("Beta Distribution")
source("ks.R")

## Beta Distribution
distibution.test("Beta(2,2)",c(0,1),dbeta,pbeta,shape1 = 2,shape2 = 2)
## Normal Distribution
distibution.test("N(9737.67,0.1)",c(-Inf,Inf),dnorm,pnorm,mean=9737.67,sd=0.1)
## Gamma Distribution
distibution.test("Gamma(shape=4)",c(0,Inf),dgamma,pgamma,shape=4)
## Chi Squared Distribution
distibution.test("chisq(7)",c(0,Inf),dchisq,pchisq,df=7)
