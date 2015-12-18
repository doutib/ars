
source("test-ks.R")

## Gamma Distribution
distibution.test("Gamma(shape=5)",c(0,Inf),dgamma,pgamma,shape=5)