
cdf <- function(d, theta= k_eval, k=theta_eval){
    if(d >0){ #we don't want to try anything with a negative distance...
        cdf <- cubintegrate(f = integrate_kernel_sum0.5, lower = 0, upper = d, k=k, theta=theta, method = "pcubature")$integral
        return(cdf)
    } else {
        return(NA)
    }
}
cdf_solve <- function(d, theta = theta_eval, k = k_eval){
    return(cdf(d, theta, k)-0.250)$x
}
#k