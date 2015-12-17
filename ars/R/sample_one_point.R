#' Sampling one point
#'
#' Can sample one point from the given piece-wise exponential density using the upper hull. The Inverse-CDF method is used here. A random number is generated, then it's used with the help of the calculated inverse CDF function of the piece-wise exponential density to generate an x.
#' @usage sample_one_point(s)
#' @param s A matrix with each row (h'(x), h(x_i),x_i,xmin,xmax,cumulative probability, normalized cumulative probability) at each node. This can be calculated using envelop_density(h,x,domain)
#' @export
sample_one_point <-
function(s)
{
  tmp<- runif(1)
  i <- sum(tmp > s[,7]) + 1 # Find i such that the ith interval contains x*.

  tmp2 <- ifelse(i==1,0,s[i-1,7]) # Cumulative probability of the first i-1 intervals.

  # If the h_prime is very close to zero, then avoid using the CDF and inv_CDF,
  # otherwise use inv_CDF to sample one x analytically.
  if (abs(s[i,1])<1e-12){
    x_sampled <- s[dim(s)[1],6]*(tmp-tmp2)/exp(s[i,2])+s[i,4]
  } else {
    CDF <- function(s,x) {
      exp(s[i,2]+(x-s[i,3])*s[i,1])/s[i,1]
    }
    inv_CDF <- function(s,y) {
      (log(y*s[i,1])-s[i,2])/s[i,1]+s[i,3]
    }
    x_sampled <- inv_CDF(s,(tmp-tmp2)*s[dim(s)[1],6]+CDF(s,s[i,4]))
  }
  return(c(x_sampled[[1]],as.integer(i)))
}
