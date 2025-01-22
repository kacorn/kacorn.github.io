


ripple.fixed <- function(
    phy,
    data,
    landscapes,
    weights,
    b1s,
    b0s,
    P = NULL,
    E = NULL,
    landscape_format = c("PCs","distribution","function"),
    k = NULL,
    kx = NULL,
    ky = NULL
){
  #make list of parameters and ensure nothing is missing
  
  print("ripple is checking data inputs.")
  
  formatted_inputs <- make.ripple.data(phy = phy, data = data, 
                                       landscapes = landscapes, 
                                       landscape_format = landscape_format,
                                       E = E, k = k, kx = kx, ky = ky,
                                       weights = weights ,b1s = b1s, b0s = b0s,
                                       P = P)
  
  params <- formatted_inputs$params
  params_to_opt <- formatted_inputs$params_to_opt
  
  if(is.null(params$data)){stop("No data provided. You must provide data.")}
  
  print("Data processed. Starting to estimate likelihood.")
  
  param_vec <- c(params_to_opt$b1s, params_to_opt$b0s, params_to_opt$xs, 
                 params_to_opt$P)
  
  n_ws <- (params$E*params$L) - params$E
  
  cache_out <- list()
  #get likelihood 
  cache_out$lnl <- likfn(param_vec = c(param_vec), params)
  print(paste("lnl =", cache_out$lnl))

  cache_out$n_free_params <- (params$E*2) + n_ws + ((params$E^2) - params$E)
  cache_out$opt_params$AIC <- (2*cache_out$n_free_params) - (2*cache_out$lnl)
  spp_k <- length(params$tree$tip.label)
  cache_out$opt_params$AICc <- cache_out$opt_params$AIC -((2*spp_k*(spp_k+1))/(cache_out$n_free_params-spp_k-1))
  cache_out$input_pars <- params
  
  
  #make q matrix from fixed parameters
  cache_out$transition_matrices <- Rs_and_Q(params = params, 
                                            params_to_opt = params_to_opt)
  
  class(cache_out) <- c("list", "ripple")
  
  return(cache_out)
  
}  
#debugonce(make.ripple.data)
check <- make.ripple.data(phy = topology_tree,
                 data = data.frame(fastBM(topology_tree, bounds = c(0,1)),
                                   fastBM(topology_tree, bounds = c(0,1)),
                                   sample(seq(1:2), 100, replace = T)),
                 landscapes = list(gaussian_landscape, bimodal_landscape),
                 landscape_format = "distribution",
                 E = 2, k = 5, kx = NULL, ky = NULL,
                 weights = c(0.5,0.5), b1s = c(0.05, 0.05), 
                 b0s = c(1, 1), P = c(0,1,1,0))  

#debugonce(Rs_and_Q)
ripple_test <- ripple.fixed(phy = topology_tree,
                      data = data.frame(fastBM(topology_tree, bounds = c(0,1)),
                                        fastBM(topology_tree, bounds = c(0,1)),
                                        sample(seq(1:2), 100, replace = T)),
                      landscapes = list(gaussian_landscape, bimodal_landscape),
                      landscape_format = "distribution",
                      E = 2, k = 5, kx = NULL, ky = NULL,
                      weights = c(0.5,0.5), b1s = c(0.05, 0.05), 
<<<<<<< HEAD
                      b0s = c(1, 1), P = c(0,1,1,0))


=======
                      b0s = c(1, 1), P = c(0,1,1,0)
)
>>>>>>> fc7230deb32e4998eadd76b0efdb1287177a2712
