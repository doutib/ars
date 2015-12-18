
source("test-ks.R")

## Gamma Distribution
distibution.test("Gamma(shape=8)",c(0,Inf),dgamma,pgamma,shape=8)