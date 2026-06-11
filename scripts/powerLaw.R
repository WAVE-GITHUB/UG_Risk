#============================================
# Power Law for whiteflies wandering 
#============================================

power_kernel <- function(d, sigma, alpha){
  
  (1 + d/sigma)^(-alpha)
  
}