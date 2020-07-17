integrate_kernel_sum1 <- function(d, k, theta){ 
    
    disp <-  d*exp(k)*theta*exp(-(exp(k)*d)^theta)/gamma(1/theta)
    return(disp)

}

