
## Envelop

# -------------------------------------------------------------
# Using the abscissae to calculate the tangent line, return a matrix with each 
# row (h'(x), h(x_i),x_i,xmin,xmax).
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

# -------------------------------------------------------------
# Using the abscissae to calculate the squeezzing lines, return a matrix with
# each row (slope, h(x_i),x_i).
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

# -------------------------------------------------------------
#  Test log concavity using u and l at abscissae
log_concavity <-function (h,x)
{
  second <- sapply(x, function(a){hessian(h,a)})
  return(!all(second < 0))
}

log_concavity <-function (u,l)
{
  n<- dim(u)[1]
  tmp=0+ifelse(u[1,1]>l[2,1],FALSE,TRUE)
  tmp=tmp+sum(sapply(2:(n-1),function(i) ifelse(u[i,1]<l[i,1] & u[i,1]>l[i+1,1],FALSE,TRUE)))
  tmp=tmp+ifelse(u[n,1]<l[n,1],FALSE,TRUE)
  return(as.logical(tmp))
}

# -------------------------------------------------------------
# Return a matrix with each row (h'(x), h(x_k),x_k,xmin,xmax,cumulative 
# probability, normalized cumulative probability).
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

