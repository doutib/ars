#' Search function
#'
#' Searches for a x0 such that h(x0) is finite and the gradient of h at x0 is attained. This is particular useful when the domain is unbounded, the mode is far away from 0 (yet unknown) and the variance is very small, i.e. the density function is very narrow. In such cases, it's even difficult to make a guess for an initial value to find x1 and xk.
#' @usage search(h,domain,x_start,max_x,min_step,relax_factor)
#' @param h The log of the density function
#' @param domain A vector of length 2 giving the domain (left bound, right bound) of the density function
#' @param x_start The starting x value for the search
#' @param max_x The maximum search limit for x, if the mode is shifted beyond this limit, the user should define this parameter
#' @param relax_factor The relaxation factor in the adaptive-step search algorithm
#' @param min_step The minimum step allowed in the adaptive-step search algorithm, if the variance is extremely small, then the used can decrease the min_step, yet it would be slower
#' @export
search <-
function(h,domain,x_start,max_x,min_step,relax_factor){
  if(is.finite(domain[1]))
  {
    x_start <- domain[1]
  }
  x_left <- x_start
  x_right <- x_start
  hvalue1 <- h(x_left)
  hvalue2 <- h(x_right)
  num_iter <-0
  step =10
  finite_value <- FALSE
  while(!finite_value)
  {
    while(is.infinite(hvalue1) & is.infinite(hvalue2)& num_iter< max_x/step)
    {
      x_left <- x_left - step
      hvalue1 <- h(x_left)
      x_right <- x_right + step
      hvalue2<- h(x_right)
      num_iter = num_iter+1
    }
    if (is.finite(hvalue1) || is.finite(hvalue2))
    {
      finite_value <- TRUE
    } else{
      x_left <-0
      x_right <-0
      num_iter <-0
      step <-max(step/relax_factor,min_step)
    }
  }
  x0 <- ifelse(is.finite(hvalue1),x_left,x_right)
  return(x0)
}
