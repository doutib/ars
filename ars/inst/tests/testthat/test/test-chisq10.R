
source("test-ks.R")

## Chi Squared Distribution
distibution.test("chisq(10)",c(0,Inf),dchisq,pchisq,df=10)
