
## Part 4: updating step

update_grid <- function(x,x_add)
{
  i <- sum(x_add > x)
  x_new <- c(x[1:i],x_add,x[(i+1):length(x)])
  return(x_new)
}