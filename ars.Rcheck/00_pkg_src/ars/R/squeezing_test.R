#' Squeezing test
#'
#' Performs the squeezing test for a newly sampled x
#' @usage squeezing_test(x_sampled,l,u)
#' @param x_sampled a newly sampled x, from sample_one_point()
#' @param u the upper hull calculated using envelop(h,x,domain)
#' @param l the lower hull calculated using squeezing(h,x)
#' @export
squeezing_test <-
function(x_sampled,l,u)
{
  w <- runif(1)
  i <- x_sampled[2]
  j <- ifelse(x_sampled[1]>u[i,3],i+1,i ) ## choose the squeezing line for x_sampled
  u_x_sampled <-u[i,2]+(x_sampled[1]-u[i,3])*u[i,1]
  l_x_sampled <- l[j,2]+(x_sampled[1]-l[j,3])*l[j,1]
  x_accept <- ifelse (w <=exp(l_x_sampled-u_x_sampled),TRUE,FALSE)
  return(x_accept)
}
