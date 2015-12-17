#' Envelope density
#' 
#' Returns a normalized density function using the upper hull. A matrix with each row (h'(x_i), h(x_i),x_i,xmin,xmax,cumulative probability, normalized cumulative probability) at each node is returned. With this infomation we can sample with sk.
#' @usage envelop_density(h,x,domain)
#' @param h The log of the density function
#' @param x The nodes, found by abscissae()
#' @param domain A vector of length 2 giving the domain (left bound, right bound) of the density function
envelop_density <-
function(h,x,domain)
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
