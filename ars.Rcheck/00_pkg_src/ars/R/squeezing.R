#' Squeezing function
#'
#' Calculates the lower hull of h. The returned value is a matrix with each row storing (slope, h(x_i),x_i) for each node. In this way, each lower hull is uniquely defined. Note that the number of squeezing lines is that of the envelop lines plus one.
#' @usage squeezing(h,x)
#' @param h The log of the density function
#' @param x The nodes, found by abscissae()
#' @export
squeezing <-
function(h,x)
{
  n <- length(x)
  lower_hull <- matrix(0,n+1,3)
  lower_hull[1,3] <- x[1]
  lower_hull[2:n,1] <- (h(x[-1])-h(x[-n]))/(x[-1]-x[-n])
  lower_hull[2:n,2] <- h(x[-n])
  lower_hull[-1,3] <- x
  return(lower_hull)
}
