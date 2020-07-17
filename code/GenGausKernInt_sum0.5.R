integrate_kernel_sum0.5 <- function(d, k, theta) { 

    z = exp(k)
    disp = exp(-(z*d)^(theta))*(theta/(2*(1/z)*gamma(1/theta)))
    return(disp)
  
  }

