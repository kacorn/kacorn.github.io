
make.ripple.data <- function(phy,
                             data = NULL,
                             landscapes,
                             weights,
                             b1s,
                             b0s,
                             P = NULL,
                             E = NULL,
                             landscape_format = c("PCs","distribution","function"),
                             k = NULL,
                             kx = NULL,
                             ky = NULL){
    params <- list(tree = phy, data = data, 
                   raw_landscapes_from_user = landscapes, 
                   landscape_format = landscape_format,
                   E = E, k = k, kx = kx, ky = ky)
    params_to_opt <- list(ws = weights, b1s = b1s, b0s = b0s, P = P)
    
    if(is.null(params$data) & is.null(params$E)){
      stop("No data or regimes found or provided. You must provide data.")
      #params$E <- 1
    } else if (is.vector(params$data) & is.null(params$E)){
      params$E <- 1
      print("Data found without regimes. Assuming one regime.")
    } else if (is.vector(params$data) & !is.null(params$E)){
      params$E <- params$E
      print(paste("Data does not contain regimes.", params$E, 
                  "regime(s) provided."))
    } else if (!is.null(params$data) & is.null(params$E)){
      params$E <- length(unique(params$data[,ncol(params$data)]))
      print(paste("Data found with ", params$E, "regimes."))
    } else if (!is.null(params$data) & !is.null(params$E)) {
      params$E <- params$E
      if(length(unique(params$data[,ncol(params$data)])) != params$E & 
         params$E != 1){
        stop("Number of regimes found in data and supplied in E do not match.")
      } else print(paste("Data found with ", params$E, "regime(s)."))
    } else {warning(paste("No data found though", params$E, 
                          "regime(s) provided. You must provide data"))}
    stopifnot(params$E > 0)
    
    #check landscape is digitized properly
    if(params$landscape_format == "PCs"){
      params$kx <- length(unique(params$raw_landscapes_from_user[,1]))
      params$ky <- length(unique(params$raw_landscapes_from_user[,2]))
      params$k2 <- params$kx * params$ky
      if( params$kx == params$ky){
        print(paste("Square discretized landscape of dimensions", 
                    params$kx, "x", params$ky, "detected."))
      } else {print(paste("Rectangular discretized landscape of dimensions", 
                          params$kx, "x", params$ky, "detected."))}
    } else if(!is.null(params$kx) & !is.null(params$ky)){
      params$k2 <- params$kx * params$ky
    } else if(!is.null(params$k)){
      params$kx <- params$k
      params$ky <- params$k
      params$k2 <- params$kx * params$ky
      print(paste("Square landscape dimensions", params$k, "x", params$ky, "provided"))
    } else {stop("No landscape discretization dimensionality provided. 
          Please provide k, kx & ky, or a pre-discretized landscape.")}
    
    params$landscapes_list <- scale_and_process_landscapes(params$raw_landscapes_from_user, params)
    #pull parameters from dataset
    params$L <- length(params$landscapes_list) #number of landscapes
    
    params$n_free_ws <- (params$E*params$L) - params$E
    #set up tile IDS
    params$phenospace_tile_ids <- tile_phenospace_ids(params)
    params$ecospace_tile_ids <- tile_ecospace_ids(params)
    
    #process data: continuous values
    if(!is.null(params$data)){
      if(length(params$data[[1]]) == 1){ 
        #data is single tile binned with no regimes provided
        params$binned_data <- Q_sim_to_xy_bins(params)
        if(all(is.na(params$binned_data$tile_eco_xy))){
          stop("Data binned incorrectly.")
        } else { print("Single tile binned data processed and ecotypes extracted.")}
      } else if(length(params$data) == 2 & params$E > 1){ 
        #data has to be single tile binned and regime provided
        params$binned_data <- pre_binned_data_to_xy(params)
        if(all(is.na(params$binned_data$tile_eco_xy))){
          stop("Data binned incorrectly.")
        } else { print("Single tile binned data and ecotypes processed.") }
      } else {if (length(params$data) == 3) 
        #ok, so we definitely have an x & y, but what are they
      { if (all(is.wholenumber(params$data[,1:2]) & params$E == 1)){
        params$binned_data <- bin_data(params) #need to write something to bin data that are xy binned already
        if(all(is.na(params$binned_data$tile_eco_xy))){
          stop("Data binned incorrectly.")
        } else {
          print("Pre-binned data processed, maybe incorrectly. 
            This is an error. Check this, Katherine.") }
      } else {params$binned_data <- bin_data(params)
      if(all(is.na(params$binned_data$tile_eco_xy))){
        stop("Data binned incorrectly.")
      } else {
        print("Continuous data binned and provided ecotypes processed.")}}}
      } 
    } else (warning("No data found. Please provide data."))
    
    #did we have to remove any outside-the-bounds data? if so, tell people to remove those from tree
    if(!is.null(params$binned_data)){
      if(nrow(params$binned_data) != nrow(as.data.frame(params$data))){
        print(paste("Data detected that is outside the bounds of the landscapes.", 
                    nrow(params$data) - nrow(params$binned_data),"species removed from dataset."))
      }
    }
    
    if(!is.null(params$binned_data)){
      if(length(params$tree$tip.label) != nrow(params$binned_data)){
        if(any(rownames(params$binned_data) %in% params$tree$tip.label)){
          params$tree <- ladderize(drop.tip(params$tree, 
                                            base::setdiff(params$tree$tip.label, rownames(params$binned_data))))
          print(paste("Data detected that is not in phylogeny.", 
                      nrow(params$data) - nrow(params$binned_data),"species removed from tree."))
        } else {print("Data in phylogeny and dataset do not match. No rownames provided in dataset to rectify.
                    You must fix this.")}
      }
    }
    
    #make deltaFs
    params$deltaF <- make_deltaF_array(params)
    
    #set class
    class(params) <- c("ripple","list")
    
    #process parameters to optimize
    
    #check number of parameters is appropriate
    if(length(params_to_opt$b1s) != params$E) {
      if(length(params_to_opt$b1s == 1)) {
        params_to_opt$b1s <- rep(params_to_opt$b1s, times = params$E)
      } else { stop(paste("Incorrect number of b1s provided.", params$E, 
                          "b1s expected, but", length(params_to_opt$b1s), "provided."))}
    }
    
    if(length(params_to_opt$b0s) != params$E) {
      if(length(params_to_opt$b0s == 1)) {
        params_to_opt$b0s <- rep(params_to_opt$b0s, times = params$E)
      } else { stop(paste("Incorrect number of b0s provided.", params$E, 
                          "b0s expected, but", length(params_to_opt$b0s), "provided."))}
    }
    
    if(length(params_to_opt$b0s) != params$E) {
      stop(paste("Incorrect number of b0s provided.", params$E, 
                 "b0s expected, but", length(params_to_opt$b0s), "provided."))
    }
    
    # if(params$E ==1){
    #   params_to_opt$P <- matrix(c(0, 1, 1, 0), nrow = 2, ncol = 2)}
    
    #is P a matrix?
    if(is.matrix(params_to_opt$P) == F & !is.null(P)){
      params_to_opt$P <- matrix(params_to_opt$P, nrow = params$E, params$E)
      diag(params_to_opt$P) <- 0
      print("don't be a lazy butt. format your freaing arguments. i just did it for u")
    }
    
    if(any(dim(params_to_opt$P) != params$E) & !is.null(P)){
      stop(paste("Incorrect dimensions of between-regime transition matrix
               provided. Expected",params$E,params$E,"but",
                 dim(params_to_opt$P),"provided."))
    }
    
    #format weights--users can submit weights for each landscape if they want
    if(length(params_to_opt$ws) == (params$E*params$L)){
      #get out each set of ws and make sure they are the right length
      this_R <- 1
      w_counter <- 1
      tot_ws <- params$E*params$L
      tot_ws_per_deltaF <- tot_ws/params$E
      free_ws_per_deltaF <- tot_ws_per_deltaF - 1
      free_ws <- numeric(0) #set up to catch weights
      for(this_R in 1:cache$E){
        these_ws <- params_to_opt$ws[(w_counter):((w_counter) + tot_ws_per_deltaF-1)]
        stopifnot(sum(these_ws) == 1)
        these_free_ws <- these_ws[1:free_ws_per_deltaF]
        free_ws <- c(free_ws, these_free_ws)
        w_counter <- w_counter + tot_ws_per_deltaF
      }
      
    } else if(length(params_to_opt$ws) == ((params$E*params$L)-params$E)){
      print("Correct number of weights provided.")
    } else{ stop(paste("Incorrect number of weights provided.",
                       ((params$E*params$L)-params$E), "weights expected, but",
                       length(params_to_opt$ws), "provided."))}
    #PUT IN TEST TO MAKE SURE WEIGHTS ARE SUMMING TO 1
    
    #get ws into a format we can work with for the transformation
    if(any(params_to_opt$ws == 1)){
      params_to_opt$ws <- replace(x <- c(params_to_opt$ws), x==1, 0.999)}
    
    if(any(params_to_opt$ws == 0)){
      params_to_opt$ws <- replace(x <- c(params_to_opt$ws), x==0, 0.001)}
    #gotta back transform ws into x
    params_to_opt$xs <- sapply(params_to_opt$ws, make_x)  
    
    ripple_pars <- list(params, params_to_opt)
    names(ripple_pars) <- c("params", "params_to_opt")
    
    return(ripple_pars)
    
}

ripple <- function(
  phy,
  data = NULL,
  landscapes,
  weights,
  b1s,
  b0s,
  P = NULL,
  E = NULL,
  landscape_format = c("PCs","distribution","function"),
  k = NULL,
  kx = NULL,
  ky = NULL,
  itnmax = NULL
){
  #make list of parameters and ensure nothing is missing
  
  print("ripple is checking data inputs.")
  
  formatted_inputs <- make.ripple.data(phy = phy, data = data, 
                 landscapes = landscapes, 
                 landscape_format = landscape_format,
                 E = E, k = k, kx = kx, ky = ky,
                 weights = weights, b1s = b1s, b0s = b0s,
                 P = P)
  
  params <- formatted_inputs$params
  params_to_opt <- formatted_inputs$params_to_opt
  
  if(is.null(params$data)){stop("No data provided. You must provide data.")}
  
   #THIS HERE
  #transformation can't handle any number that's actually larger than or equal to 1 or less than or equal to 0
  
  print("Data processed. Starting to optimize parameters")
  
  param_vec <- c(params_to_opt$b1s, params_to_opt$b0s, params_to_opt$xs, 
                 params_to_opt$P)
  
  n_ws <- (params$E*params$L) - params$E
  
  optimizing <- optimx(param_vec, params = params, fn = likfn, 
                       method = "Nelder-Mead")
  
  cache_out <- list()
  #save optimizer results
  cache_out$optimizer_pars$convcode <- optimizing$convcode
  cache_out$optimizer_pars$n_iterations <- optimizing$fevals
  cache_out$optimizer_pars$exe_time <- optimizing$xtime
  #cache_out$optimizer_pars$warnings <- NA
  
  #save optim params
  cache_out$opt_params$b1s <- (optimizing[1:params$E]) #so that this function returns b1, not p1
  cache_out$opt_params$b0s <- (optimizing[(params$E+1):(params$E*2)])
  xs <- optimizing[((params$E*2)+1):((params$E*2)+n_ws)] #E*L - E
  cache_out$opt_params$ws <- sapply(xs, make_w)  #w = exp(x) / (1 + exp(x))
  if(params$E > 1){
    cache_out$opt_params$P <- matrix(optimizing[((params$E*2)+n_ws+1):length(param_vec)],
                                     nrow = params$E, ncol = params$E, byrow = T)
  } else {cache_out$opt_params$P <- 1}
  
  print("pars opt works")
  
  #get likelihood of optimized parameters quietly
  invisible(capture.output(get_lik <- ripple.fixed(phy = phy,
                             data = data,
                             landscapes = landscapes,
                             E = params$E, k = k, kx = kx, ky = ky,
                             landscape_format = landscape_format,
                             weights = as.vector(cache_out$opt_params$ws),
                             b1s = c(unlist(cache_out$opt_params$b1s)), 
                             b0s = c(unlist(cache_out$opt_params$b0s)),
                             P = matrix(c(unlist(cache_out$opt_params$P)), 
                                        nrow = params$E))))
  
  cache_out$opt_params$lnl <- get_lik$lnl
  cache_out$transition_matrices <- list()
  cache_out$opt_params$Rs <- get_lik$transition_matrices$Rs
  cache_out$opt_params$Q <- get_lik$transition_matrices$Q
  
  cache_out$n_free_params <- (params$E*2) + n_ws + ((params$E^2) - params$E)
  cache_out$opt_params$AIC <- (2*cache_out$n_free_params) - (2*cache_out$opt_params$lnl)
  spp_k <- length(params$tree$tip.label)
  cache_out$opt_params$AICc <- cache_out$opt_params$AIC -((2*spp_k*(spp_k+1))/(cache_out$n_free_params-spp_k-1))
  cache_out$input_pars <- params
  
  class(cache_out) <- c("list", "ripple")
  
  return(cache_out)
  
}

#make landscapes
gaussian_landscape <- as.data.frame(cbind(x = rnorm(n = 1000, mean = 0, sd = 1.5),
                                          y = rnorm(n = 1000, mean = 0, sd = 1.5)))
bimodal_landscape <- as.data.frame(cbind(x = c(rnorm(n = 500, mean = 2, sd = 1), rnorm(n = 500, mean = -2, sd = 1)),
                                         y = c(rnorm(n = 500, mean = 2, sd = 1), rnorm(n = 500, mean = -2, sd = 1))))
#I expect you to come in with this
raw_landscapes_from_user <- list(gaussian_landscape, bimodal_landscape)

#make tree
tree <- geiger::drop.extinct(geiger::sim.bdtree(b=1, d=0.9, 
                                                stop="taxa", n=100, seed=1, extinct=TRUE))
topology_tree <- tree
topology_tree$edge.length <- (topology_tree$edge.length/max(ape::branching.times(topology_tree)))* 10000

require(phytools)


ripple_test <- ripple(phy = topology_tree,
                      data = data.frame(fastBM(topology_tree, bounds = c(0,1)),
                                        fastBM(topology_tree, bounds = c(0,1)),
                                        sample(seq(1:2), 100, replace = T)),
                      landscapes = list(gaussian_landscape, bimodal_landscape),
                      landscape_format = "distribution",
                      E = 2, k = NULL, kx = 5, ky = 6,
                      weights = c(0.5,0.5), b1s = c(0.05, 0.05), 
                      b0s = c(1, 1), P = c(0,1,1,0))

ripple_test <- ripple.fixed(phy = topology_tree,
                      data = data.frame(fastBM(topology_tree, bounds = c(0,1)),
                                        fastBM(topology_tree, bounds = c(0,1)),
                                        sample(seq(1:2), 100, replace = T)),
                      landscapes = list(gaussian_landscape, bimodal_landscape),
                      landscape_format = "distribution",
                      E = 2, kx = 5, ky = 6,
                      weights = c(0.5,0.5), b1s = c(0.05, 0.05), 
                      b0s = c(1, 1), P = c(0,1,1,0))

ripple_check <- ripple.fixed(phy = topology_tree,
                             data = data.frame(fastBM(topology_tree, bounds = c(0,1)),
                                               fastBM(topology_tree, bounds = c(0,1)),
                                               sample(seq(1:2), 100, replace = T)),
                             landscapes = list(gaussian_landscape, bimodal_landscape),
                             landscape_format = "distribution",
                             E = 2, k = 5, kx = NULL, ky = NULL,
                             weights = c(0.5,0.5), b1s = c(0.05, 0.05), 
                             b0s = c(1, 1), P = c(0,1,1,0))

transition_matrices_test <- Rs_and_Q(ripple_test$input_pars, ripple_test$input_pars)

ripple_test$opt_params
