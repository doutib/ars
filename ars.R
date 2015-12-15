
## ars

# Load functions
require(pracma)
source("Initialization.R")
source("Sampling.R")
source("Update.R")
source("Envelop.R")

# Ars

ars <- function(h,n,domain,x1 = NULL,xk = NULL,x0=0,nmesh=5,
                x_start=0,max_x=10000,min_step=0.01,relax_factor=10)
{
  n_sampled <- 0
  x_sample <- vector("numeric")
  first_run <-TRUE
  while(n_sampled < n)
  {
    if(first_run)
    { 
      if(is.infinite(domain[1])& is.infinite(domain[2]))
      {
        tmp <- search(h,c(-Inf,Inf),x_start,max_x,min_step,relax_factor)
        x <- abscissae(h,c(-Inf,Inf),x1,xk,x0=tmp,nmesh,min_step,relax_factor)
      } else{x <- abscissae(h,domain,x1,xk,x0,nmesh,min_step,relax_factor)}
      u <- envelop(h,x,domain)
      l <- squeezing(h,x)
      if(log_concavity(u,l))
      {
        warning('local log non-concavity detected,check input density function!')
      }
      s <- envelop_density(h,x,domain)
    }
    x_sampled <-sample_one_point(s)
    
    pass_squeezing <- squeezing_test(x_sampled,l,u)
    if(!pass_squeezing)
    {
      pass_rejection<-rejection_test(x_sampled,u,h)
    }
    if (pass_squeezing || pass_rejection)
    {
      n_sampled=n_sampled+1
      x_sample[n_sampled]<-x_sampled[1]
      if (!pass_squeezing)
      {
        x <-update_grid(x,x_sampled[1])
        u <- envelop(h,x,domain)
        l <- squeezing(h,x)
        if(log_concavity(u,l))
        {
          warning('local log non-concavity detected,check input density function!')
        }
        s <- envelop_density(h,x,domain)  
      }
    }
    first_run <- FALSE
  }
  return(x_sample)
}