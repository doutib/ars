
## Kolmogorov-smirnov test

# Test function --------------------------------------------
  # name: description of the test
  # domain: supprt of the cdf
  # distribution: dnorm, dgamma...
  # true_cdf: pnorm, pgamma...
  # n: number of samples
  # pvalue: threshold of the pvalue of the test 
distibution.test = function(name,domain,distribution,true_cdf,n=10000,pvalue=0.05,...)
{
  test_that(name, {
    h <- function(x){
      log(distribution(x,...))
    }
    # Run
    cat('\nComputing test: ',name)
    samples <- ars(h,n,domain)
    ks_test=ks.test(samples,true_cdf,...)
    
    # Test
    expect_that( ks_test$p.value >= pvalue, is_true() )
    
    # Plot
    hist(samples,freq = F)
    lines(sort(samples),distribution(sort(samples),...))
    
    # Print Result
    cat("\nResult:\nD = ",ks_test$statistic,"\np-value = ",ks_test$p.value,"\n------------")
  })
  
}






