#' Rejection test
#'
#' Performs the rejection test for a newly sampled x if it fails the squeezing test.
#' @usage rejection_test(x_sampled,u,h)
#' @param x_sampled a newly sampled x, from sample_one_point()
#' @param u the upper hull calculated using envelop(h,x,domain)
#' @param h the log density function
#' @export
rejection_test <-
function(x_sampled,u,h)
{
  w <- runif(1)
  i <- x_sampled[2]
  u_x_sampled <-u[i,2]+(x_sampled[1]-u[i,3])*u[i,1]
  h_x_sampled <- h(x_sampled[1])
  x_accept <- (w <=exp(h_x_sampled-u_x_sampled))
  return(x_accept)
}
