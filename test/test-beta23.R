
source("test-ks.R")

## Beta Distribution
distibution.test("Beta(2,3)",c(0,1),dbeta,pbeta,shape1 = 2,shape2 = 3)
