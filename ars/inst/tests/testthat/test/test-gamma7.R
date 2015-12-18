
source("test-ks.R")

## Gamma Distribution
distibution.test("Gamma(shape=7)",c(0,Inf),dgamma,pgamma,shape=7)