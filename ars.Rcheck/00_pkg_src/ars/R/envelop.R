#' Envelope function
#'
#' Calculates the upper hull of h. The returned value is a matrix with each row storing (h'(x), h(x_i),x_i,xmin,xmax) for each node. In this way, each upper hull is uniquely defined, as well as its corresponding domain.
#' @usage envelop(h,x,domain)
#' @param h The log of the density function
#' @param x The nodes, found by abscissae()
#' @param domain A vector of length 2 giving the domain (left bound, right bound) of the density function
#' @export
envelop <-
function(h, x, domain)
{
  # h'(x) and h(x_k):
  temp <- cbind ( sapply(x, function(a) pracma::grad(h,a)), h(x)) #col 1 thru 2
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
