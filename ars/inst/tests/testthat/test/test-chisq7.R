
source("test-ks.R")

## Chi Squared Distribution
distibution.test("chisq(7)",c(0,Inf),dchisq,pchisq,df=7)
