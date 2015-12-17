
source("test-ks.R")

## Chi Squared Distribution
distibution.test("chisq(8)",c(0,Inf),dchisq,pchisq,df=8)
