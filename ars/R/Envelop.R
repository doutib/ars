#' Envelope function
#' 
#' Calculates the upper hull of h. The returned value is a matrix with each row storing (h'(x), h(x_i),x_i,xmin,xmax) for each node. In this way, each upper hull is uniquely defined, as well as its corresponding domain.
#' @usage envelop(h,x,domain)
#' @param h The log of the density function
#' @param x The nodes, found by abscissae()
#' @param domain A vector of length 2 giving the domain (left bound, right bound) of the density function


envelop <- function(h, x, domain)
{
  # h'(x) and h(x_k):
  temp <- cbind ( sapply(x, function(a) grad(h,a)), h(x)) #col 1 thru 2 
  # intersections:
  l=length(x)
  intersect <- (temp[-1,2] - temp[-l,2] + x[-l]*temp[-l,1] - x[-1]*temp[-1,1]) / (temp[-l, 1] - temp[-1, 1])
  # xmin/max from intersections:
  xmin <- c(domain[1], intersect)
  xmax <- c(intersect, domain[2]) 
  #format, output:
  upper_hull <- cbind (temp, x, xmin, xmax)
  return(upper_hull)
}

#' Squeezing function
#' 
#' Calculates the lower hull of h. The returned value is a matrix with each row storing (slope, h(x_i),x_i) for each node. In this way, each lower hull is uniquely defined. Note that the number of squeezing lines is that of the envelop lines plus one.
#' @usage squeezing(h,x)
#' @param h The log of the density function
#' @param x The nodes, found by abscissae()

squeezing <- function(h,x)
{
  n <- length(x)
  lower_hull <- matrix(0,n+1,3)
  lower_hull[1,3] <- x[1]
  lower_hull[2:n,1] <- (h(x[-1])-h(x[-n]))/(x[-1]-x[-n])
  lower_hull[2:n,2] <- h(x[-n])
  lower_hull[-1,3] <- x
  return(lower_hull)
}

#' Log concavity test
#' 
#' Checks that at each node, the upper hull and lower hull would bound the h function. The slope of the envelop lines and squeezing lines are examined.
#' @usage log_concavity(u,l)
#' @param u the upper hull calculated using envelop(h,x,domain)
#' @param l the lower hull calculated using squeezing(h,x)


log_concavity <-function (u,l)
{
  n<- dim(u)[1]
  tmp=0+ifelse(u[1,1]>l[2,1],FALSE,TRUE)
  tmp=tmp+sum(sapply(2:(n-1),function(i) ifelse(u[i,1]<l[i,1] & u[i,1]>l[i+1,1],FALSE,TRUE)))
  tmp=tmp+ifelse(u[n,1]<l[n,1],FALSE,TRUE)
  return(as.logical(tmp))
}

#' Envelope density
#' 
#' Returns a normalized density function using the upper hull. A matrix with each row (h'(x_i), h(x_i),x_i,xmin,xmax,cumulative probability, normalized cumulative probability) at each node is returned. With this infomation we can sample with sk.
#' @usage envelop_density(h,x,domain)
#' @param h The log of the density function
#' @param x The nodes, found by abscissae()
#' @param domain A vector of length 2 giving the domain (left bound, right bound) of the density function

envelop_density <- function(h,x,domain)
{
  n <- length(x)
  u <- envelop(h,x,domain)
  integrate <- function (u){
    exp(u[2]+(u[5]-u[3])*u[1])/u[1]-exp(u[2]+(u[4]-u[3])*u[1])/u[1]
  }
  integral <- sapply (1:nrow(u), function(i) integrate (u[i,]))
  
  #take care of zeros 
  zeros <- which (abs(u[,1]) < 1e-12)
  integral[zeros]<-exp(u[zeros,2])*(u[zeros,5]-u[zeros,4])
  
  cum_probability <- cumsum(integral)
  norm_probability<- cum_probability/sum(integral)
  return(cbind(u,cum_probability,norm_probability))
}

