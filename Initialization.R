
## Initialization step.


# -------------------------------------------------------------
# The default max_x is the 10000, the default min_step is 0.1.
# They define the search area. search() search for a 
# x0 with finite value to start, the default x_start is 0 for
# unbounded functions.
search <- function(h,domain,x_start,max_x,min_step,relax_factor){
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


# -------------------------------------------------------------
# abscissae() initialize the grid and check for input validity
abscissae <- function(h,domain,x1 = NULL,xk = NULL,x0=0,nmesh=5,min_step=0.001,relax_factor=5)
{
  if(!is.null(x1) && !is.null(xk)){
    return( x1+(c(1:nmesh)-1)*(xk-x1)/(nmesh-1))
  } else{
    # Check the domain is valid.
    if (domain[1]>=domain[2])
    {
      warning('invalid domain');
    }
    xk_found <- FALSE
    x1_found <- FALSE
    step2<-1
    tol1 <- 0.001 
    tol2<- 1000 # tol1 and tol2 give a bound for the slopes, for numerical considerations.  
    if (is.infinite(domain[1])&is.infinite(domain[2]))
    {  
      hprime <- grad(h,x0)
      hvalue <- h(x0)
      x_min <- x0
      x_max <- x0
      hprime_min <- hprime
      hprime_max <- hprime
      hvalue_max <- hvalue
      hvalue_min <- hvalue
      region1 <- (hprime>= tol2)
      region2 <- (hprime> -tol2 && hprime< tol2)
      region3 <- (hprime<= -tol2)
      num_iter <-0
      # Find a xk such that the slope is bounded by -tol1 and -tol2, if not found
      # after 1000 iterations, then relax the step by dividing relax_factor, then
      # restart looking for xk.
      while(!xk_found)
      {
        while((hprime_max <=-tol2|| hprime_max >= -tol1 ||is.na(hprime_max)||is.infinite(hvalue_max))& num_iter <1000)
        {
          x_max <- x_max - step2*region3 + step2*region1 + step2*region2
          hprime_max <- grad(h,x_max)
          hvalue_max <- h(x_max)
        }
        if (!is.na(hprime_max) & ! is.infinite(hvalue_max)& hprime_max > -tol2 & hprime_max < -tol1)
        {
          xk_found <- TRUE
        } else{
          num_iter <-0
          step2 <-max(step2/relax_factor,min_step)
          x_max <-x0
        }
      }
      # Find a x1 such that the slope is bounded by tol1 and tol2, if not found
      # after 1000 iterations, then relax the step by dividing relax_factor, then
      # restart looking for x1.
      num_iter <-0
      while(!x1_found)
      {
        while((hprime_min >= tol2 || hprime_min <= tol1 ||is.na(hprime_min)||is.infinite(hvalue_min))& num_iter <1000)
        {
          x_min <- x_min - step2*region3 + step2*region1-step2*region2
          hprime_min <- grad(h,x_min)
          hvalue_min <- h(x_min)
        }
        if (!is.na(hprime_min) & ! is.infinite(hvalue_min)& hprime_min < tol2 & hprime_min > tol1)
        {
          x1_found <- TRUE
        } else{
          num_iter <-0
          step2 <-max(step2/relax_factor,min_step)
          x_min <-x0
        }
      }
    }
    # -------------------------------------------------------------
    # case2: bounded from left.
    if (is.finite(domain[1])&is.infinite(domain[2]))
    {
      x_min <- domain[1]
      x_max <- domain[1]
      hprime_min <- grad(h,x_min)
      hprime_max <- grad(h,x_max)
      hvalue_max <- h(x_min)
      hvalue_min <- h(x_max)
      num_iter <-0
      # Find a xk such that the slope is bounded by -tol1 and -tol2, if not found
      # after 1000 iterations, then relax the step by dividing relax_factor, then
      # restart looking for xk. In this case, xk should not be too small to avoid
      # sampling biase, nor should it be too large to avoid numerical issues, 
      # therefore the average of the minimum and the maximum possible xk is used.
      # However, the choice of xk really depends on how flat the density is, it's 
      # not easy to tune a algorithm for xk that works fine for all distributions.
      
      while(!xk_found)
      {
        while((hprime_max <=-tol2|| hprime_max >= -tol1 ||is.na(hprime_max)||is.infinite(hvalue_max))& num_iter <1000)
        {
          x_max <- x_max +step2
          hprime_max <- grad(h,x_max)
          hvalue_max <- h(x_max)
        }
        if (!is.na(hprime_max) & ! is.infinite(hvalue_max)& hprime_max > -tol2 & hprime_max < -tol1)
        {
          xk_found <- TRUE
        } else{
          num_iter <-0
          step2 <-max(step2/relax_factor,min_step)
          x_max <-domain[1]
        }
      }
      # Save the minimum possible xk first, then search for the maximum possible one
      x_max_min <-x_max  
      num_iter=0
      while((!is.na(hprime_max)&is.finite(hvalue_max)&hprime_max >-tol2 & hprime_max < -tol1)& num_iter <1000)
      {
        x_max <- x_max + 1
        hprime_max <- grad(h,x_max)
        hvalue_max <- h(x_max)
        num_iter<-num_iter + 1
      }
      # Find a xk lie between the min and max possible ones
      x_max <-(x_max+x_max_min)/2 
      num_iter=0
      while(!x1_found)
      {
        while((abs(hprime_min) <= tol1 || abs(hprime_min) >= tol2 ||is.na(hprime_min)||is.infinite(hvalue_min))& num_iter <1000)
        {
          x_min <- x_min +step2
          hprime_min <- grad(h,x_min)
          hvalue_min <- h(x_min)
        }
        if (!is.na(hprime_min) & ! is.infinite(hvalue_min)& abs(hprime_min) > tol1 & abs(hprime_min) < tol2)
        {
          x1_found <- TRUE
        } else{
          num_iter <-0
          step2 <-max(step2/relax_factor,min_step)
          x_min <-domain[1]
        }
      }
    }
    # -------------------------------------------------------------
    # case3: bounded from right
    if (is.finite(domain[2])&is.infinite(domain[1]))
    {
      x_min <- domain[2]
      x_max <- domain[2]
      hprime_min <- grad(h,x_min)
      hprime_max <- grad(h,x_max)
      hvalue_max <- h(x_min)
      hvalue_min <- h(x_max)
      num_iter <-0
      # Find a xk such that the slope is bounded by -tol2 and -tol1
      while(!xk_found)
      {
        while((abs(hprime_max) <= tol1|| abs(hprime_max) >= tol2 ||is.na(hprime_max)||is.infinite(hvalue_max))& num_iter <1000)
        {
          x_max <- x_max - step2
          hprime_max <- grad(h,x_max)
          hvalue_max <- h(x_max)
        }
        if (!is.na(hprime_max) & ! is.infinite(hvalue_max)& abs(hprime_max) > tol1 & abs(hprime_max) < tol2)
        {
          xk_found <- TRUE
        } else{
          num_iter <-0
          step2 <-max(step2/relax_factor,min_step)
          x_max <-domain[2]
        }
      }
      num_iter <-0
      while(!x1_found)
      {
        while((hprime_min <= tol1 || hprime_min >= tol2 ||is.na(hprime_min)||is.infinite(hvalue_min))& num_iter <1000)
        {
          x_min <- x_min -step2
          hprime_min <- grad(h,x_min)
          hvalue_min <- h(x_min)
        }
        if (!is.na(hprime_min) & ! is.infinite(hvalue_min)& hprime_min > tol1 & hprime_min < tol2)
        {
          x1_found <- TRUE
        } else{
          num_iter <-0
          step2 <-max(step2/relax_factor,min_step)
          x_min <-domain[2]
        }
      }
      num_iter <-0
      x_min_max <- x_min
      while(!is.na(hprime_min) & is.finite(hvalue_min)& hprime_min > tol1 & hprime_min < tol2& num_iter <1000)
      {
        x_min <- x_min -1
        hprime_min <- grad(h,x_min)
        hvalue_min <- h(x_min)
        num_iter <-num_iter+1
      }
      x_min <- (x_min+x_min_max)/2
    }
    # -------------------------------------------------------------
    # case4 :bounded from both sides
    if (is.finite(domain[1])&is.finite(domain[2]))
    {
      x_min <- domain[1]
      x_max <- domain[2]
      hprime_min <- grad(h,x_min)
      hprime_max <- grad(h,x_max)
      hvalue_max <- h(x_min)
      hvalue_min <- h(x_max)
      num_iter <-0
      step2 <- (domain[2]-domain[1])/100
      # find a xk such that the slope is bounded by -tol2 and -tol1
      while(!xk_found)
      {
        while((abs(hprime_max) <=tol1|| abs(hprime_max) >= tol2 ||is.na(hprime_max)||is.infinite(hvalue_max))& num_iter <1000)
        {
          x_max <- x_max - step2
          hprime_max <- grad(h,x_max)
          hvalue_max <- h(x_max)
        }
        if (!is.na(hprime_max) & ! is.infinite(hvalue_max)& abs(hprime_max) > tol1 & abs(hprime_max) < tol2)
        {
          xk_found <- TRUE
        } else{
          num_iter <-0
          step2 <-max(step2/relax_factor,min_step)
          x_max <-domain[2]
        }
      }
      # Find a x1 such that the slope not too close to zero
      num_iter <-0
      while(!x1_found)
      {
        while((abs(hprime_min) <= tol1 || abs(hprime_min) >= tol2 ||is.na(hprime_min)||is.infinite(hvalue_min))& num_iter <1000)
        {
          x_min <- x_min +step2
          hprime_min <- grad(h,x_min)
          hvalue_min <- h(x_min)
        }
        if (!is.na(hprime_min) & ! is.infinite(hvalue_min)& abs(hprime_min) > tol1 & abs(hprime_min) < tol2)
        {
          x1_found <- TRUE
        } else{
          num_iter <-0
          step2 <-max(step2/relax_factor,min_step)
          x_min <-domain[1]
        }
      }
    }
    x <- x_min+(c(1:nmesh)-1)*(x_max-x_min)/(nmesh-1)
    return(x)
  }  
}