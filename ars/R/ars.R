#' Adaptive rejection sampler
#' The main ars() function
#' @param h The log of the density function
#' @param n The number of points to be sampled
#' @param domain A vector of length 2 giving the domain (left bound, right bound) of the density function
#' @param x1 The leftmost node of the abscissae, if not provided, the abscissae() will try to find one which has a suitable slope, leading to a unbiased and numerical-issues friendly sampling process
#' @param xk The the rightmost node of the abscissae, if not provided, the abscissae() will try to find one which has a suitable slope, leading to a unbiased and numerical-issues friendly sampling process
#' @param max_x description to be added
#' @param min_step description to be added
#' @param relax_factor description to be added
#' @param x0 description to be added
#' @param nmesh description to be added
#' @param x_start description to be added
#' @return samples description to be completed
#' @import pracma
#' @examples
#' ars(function (x){log (dnorm (x, 0, 1))}, 10000, c(-Inf, Inf))
#' @export

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
