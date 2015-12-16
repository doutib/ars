
#' Sampling one point
#' 
#' Can sample one point from the given piece-wise exponential density using the upper hull. The Inverse-CDF method is used here. A random number is generated, then it's used with the help of the calculated inverse CDF function of the piece-wise exponential density to generate an x.
#' @usage sample_one_point(s)
#' @param s A matrix with each row (h'(x), h(x_i),x_i,xmin,xmax,cumulative probability, normalized cumulative probability) at each node. This can be calculated using envelop_density(h,x,domain)

sample_one_point <- function(s)
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


#' Squeezing test
#' 
#' Performs the squeezing test for a newly sampled x
#' @usage squeezing_test(x_sampled,l,u)
#' @param x_sampled a newly sampled x, from sample_one_point()
#' @param u the upper hull calculated using envelop(h,x,domain)
#' @param l the lower hull calculated using squeezing(h,x)

squeezing_test <- function(x_sampled,l,u)
{
  w <- runif(1)
  i <- x_sampled[2]
  j <- ifelse(x_sampled[1]>u[i,3],i+1,i ) ## choose the squeezing line for x_sampled
  u_x_sampled <-u[i,2]+(x_sampled[1]-u[i,3])*u[i,1]
  l_x_sampled <- l[j,2]+(x_sampled[1]-l[j,3])*l[j,1]
  x_accept <- ifelse (w <=exp(l_x_sampled-u_x_sampled),TRUE,FALSE)
  return(x_accept)
}

#' Rejection test
#' 
#' Performs the rejection test for a newly sampled x if it fails the squeezing test.
#' @usage rejection_test(x_sampled,u,h)
#' @param x_sampled a newly sampled x, from sample_one_point()
#' @param u the upper hull calculated using envelop(h,x,domain)
#' @param h the log density function

rejection_test <- function(x_sampled,u,h)
{
  w <- runif(1)
  i <- x_sampled[2]
  u_x_sampled <-u[i,2]+(x_sampled[1]-u[i,3])*u[i,1]
  h_x_sampled <- h(x_sampled[1])
  x_accept <- (w <=exp(h_x_sampled-u_x_sampled))
  return(x_accept)
}


