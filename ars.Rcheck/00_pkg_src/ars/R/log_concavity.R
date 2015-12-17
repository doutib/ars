#' Log concavity test
#'
#' Checks that at each node, the upper hull and lower hull would bound the h function. The slope of the envelop lines and squeezing lines are examined.
#' @usage log_concavity(u,l)
#' @param u the upper hull calculated using envelop(h,x,domain)
#' @param l the lower hull calculated using squeezing(h,x)
#' @export
log_concavity <-
function (u,l)
{
  n<- dim(u)[1]
  tmp=0+ifelse(u[1,1]>l[2,1],FALSE,TRUE)
  tmp=tmp+sum(sapply(2:(n-1),function(i) ifelse(u[i,1]<l[i,1] & u[i,1]>l[i+1,1],FALSE,TRUE)))
  tmp=tmp+ifelse(u[n,1]<l[n,1],FALSE,TRUE)
  return(as.logical(tmp))
}
