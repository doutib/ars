
source("test-ks.R")

## Beta Distribution
distibution.test("Beta(4,2)",c(0,1),dbeta,pbeta,shape1 = 4,shape2 = 2)
