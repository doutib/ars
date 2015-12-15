
## Sampling step

# -------------------------------------------------------------
# Use the inverse of CDF of the envelop density function to sample x*.
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

# Following the hints given by Chris for code efficiency, sampling, squeezing test and 
# rejection test are isolated so that we can sample many points first and then perform 
# tests
# note that the input x_sampled has 2 elements, x itself and interval that x locates
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

rejection_test <- function(x_sampled,u,h)
{
  w <- runif(1)
  i <- x_sampled[2]
  u_x_sampled <-u[i,2]+(x_sampled[1]-u[i,3])*u[i,1]
  h_x_sampled <- h(x_sampled[1])
  x_accept <- (w <=exp(h_x_sampled-u_x_sampled))
  return(x_accept)
}


