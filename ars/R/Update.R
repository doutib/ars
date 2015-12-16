#' Updating step
#' 
#' UPdates the grid if both h and the gradient of h is evaluated at a sampled point.#' @usage envelop_density(h,x,domain)
#' @param x The abscissae in the previous step
#' @param x_add An x value to be added to the abscissae

update_grid <- function(x,x_add)
{
  i <- sum(x_add > x)
  x_new <- c(x[1:i],x_add,x[(i+1):length(x)])
  return(x_new)
}