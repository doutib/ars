
source("test-ks.R")

## Beta Distribution
distibution.test("Beta(3,2)",c(0,1),dbeta,pbeta,shape1 = 3,shape2 = 2)
