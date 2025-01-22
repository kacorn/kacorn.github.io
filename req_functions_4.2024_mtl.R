###############################################
# Evolution across kinematic landscapes       #
# K Corn Mar. 2023                            #
# Corn, Holzman, Martin, Uyeda                #
# Data, package, and function setup           #
###############################################

# here::i_am("landscapes_git/scripts/req_functions_4.2024.R")
# require(here)
# here() #check wd

require(castor) ; require(MASS) ; require(ape); require(optimx)

ripple <- setClass("ripple")

#to rotate our images
rotate_image <- function(x) t(apply(x, 1, rev))

#some utilities
minpositive <- function(x){ 
  min(x[x > 0]) }

is.NAOb <- function(x) is.na(x) | all(sapply(x, is.na))

make_x <- function(w) {
  x = log(-1 * (w / (w-1)))
  return(x)}

make_w <- function(x) {
  w = exp(x) / (1 + exp(x))
  return(w)}

rescale_0.2to0.8 <- function(x) {
  x_norm <- ((0.8 - 0.2)*((x - min(x)) / (max(x) - min(x))))+0.2
}

rescale_0to1 <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

makeTransparent <- function(someColor, alpha=100){
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))

#discretizer
mass_discretizer <- function(unkerneled_list, params){
  landscape_counter <- 1
  kerneled_landscapes <- list()
  
  x_list <- unlist(lapply(unkerneled_list, function(x) x$x))
  y_list <- unlist(lapply(unkerneled_list, function(x) x$y))
  
  x_range <- diff(range(x_list)) ; y_range <- diff(range(y_list))
  x_range_plus <- x_range + 0.1*x_range; y_range_plus <- y_range + 0.1*y_range
  max_range <- max(x_range_plus, y_range_plus)
  x_add_on <- (max_range - x_range_plus)/2; y_add_on <- (max_range - y_range_plus)/2
    
  upper_x <- max(x_list) + 0.05*x_range + x_add_on
  lower_x <- min(x_list) - 0.05*x_range - x_add_on
  upper_y <- max(y_list) + 0.05*y_range + y_add_on
  lower_y <- min(y_list) - 0.05*y_range - y_add_on
  
  #check
  stopifnot(round(diff(range(upper_x, lower_x)), 10) == round(diff(range(upper_y, lower_y)), 10))
  
  #don't forget that the landscape has to be square so there has to be a bunch of 0s basically off the edges
  for(landscape_counter in 1:length(unkerneled_list)){
    kerneled_landscape <- MASS::kde2d(unkerneled_list[[landscape_counter]]$x,
                                      unkerneled_list[[landscape_counter]]$y,
                                      n = params$kx,
                                      lims = c(lower_x, upper_x, lower_y, upper_y)) #c(xl, xu, yl, yu)
    kerneled_landscapes[[landscape_counter]] <- kerneled_landscape
    landscape_counter <- landscape_counter + 1
  }
  return(kerneled_landscapes)
}

global_rescale_kerneled <- function(kerneled_list_of_landscapes){
  kx <- length(kerneled_list_of_landscapes[[1]]$x)-1
  ky <- length(kerneled_list_of_landscapes[[1]]$y)-1
  n_landscapes <- length(kerneled_list_of_landscapes)
  
  #first rescale each landscape 0-1
  landscape_counter = 1
  for(landscape_counter in 1:n_landscapes){
    kerneled_list_of_landscapes[[landscape_counter]]$z <- rescale_0to1(kerneled_list_of_landscapes[[landscape_counter]]$z)
    landscape_counter <- landscape_counter + 1
  }
  
  combined_zs <- data.frame()
  landscape_counter = 1
  for(landscape_counter in 1:n_landscapes){
    combined_zs <- rbind(combined_zs, kerneled_list_of_landscapes[[landscape_counter]]$z)
    landscape_counter <- landscape_counter + 1
  }
  rescaled_z <- rescale_0to1(combined_zs)
  
  #re separate the zs 
  separation_counter <- 1
  row_counter <- 1
  for(separation_counter in 1:n_landscapes){
    kerneled_list_of_landscapes[[separation_counter]]$z <- 
      rescaled_z[row_counter:(row_counter + ky - 1),]
    kerneled_list_of_landscapes[[separation_counter]]$x <- 
      seq(from = min(kerneled_list_of_landscapes[[separation_counter]]$x), 
          to= max(kerneled_list_of_landscapes[[separation_counter]]$x), 
          length.out = (kx+1))
    kerneled_list_of_landscapes[[separation_counter]]$y <- 
      seq(from = min(kerneled_list_of_landscapes[[separation_counter]]$y), 
          to= max(kerneled_list_of_landscapes[[separation_counter]]$y), 
          length.out = (ky+1))
    #fyi the kerneler in mass makes the x and ys count = k, but i want them to be the breakpoints for my bins so they need to be k+1
    row_counter <- row_counter + ky
    separation_counter <- separation_counter + 1
  }
  return(kerneled_list_of_landscapes)
}

make_deltaF_array <- function(params){
  
  L <- params$L
  landscapes_list <- params$landscapes_list
  n_landscapes <- L
  
  deltaF_list <- list()
  
  n_singlegrid_states = params$k2
  n_singlegrid_rates = (params$k2)^2
  n_singlegrid_tiles = params$k2

  landscape_counter = 1
  
  for(landscape_counter in 1:L){
    
    deltaF <- matrix(NA, nrow = n_singlegrid_states, ncol = n_singlegrid_states)
    
    index_tile_names <- params$phenospace_tile_ids
    colnames(deltaF) <- index_tile_names$tile_pheno_xy
    rownames(deltaF) <- index_tile_names$tile_pheno_xy
    
    this_tile_counter = 16
    comp_tile_counter = 15
    
    for(this_tile_counter in 1:n_singlegrid_tiles){
      
      this_tile_info = index_tile_names[this_tile_counter,]
      this_tile_coords <- c(this_tile_info$x_pheno, 
                            this_tile_info$y_pheno)
      this_tile_height <- landscapes_list[[landscape_counter]]$z[this_tile_info$y_pheno, 
                                                                 this_tile_info$x_pheno]
      
      for(comp_tile_counter in 1:n_singlegrid_tiles){
        
        comp_tile_info = index_tile_names[comp_tile_counter,]
        comp_tile_coords <- c(comp_tile_info$x_pheno, 
                              comp_tile_info$y_pheno)
        comp_tile_height <- landscapes_list[[landscape_counter]]$z[comp_tile_info$y_pheno, 
                                                                   comp_tile_info$x_pheno]
        
        if(all(abs(this_tile_coords - comp_tile_coords) < 2)){ 
          #ensure that tiles are adjacent
          
          if(any(this_tile_coords - comp_tile_coords) != 0){ 
            
            deltaF[this_tile_info$tile_pheno_xy, 
                   comp_tile_info$tile_pheno_xy] <- 
              this_tile_height - comp_tile_height  
            #MUST be i to j or will make peaks into calderas
            
          }
          else{
            deltaF[this_tile_info$tile_pheno_xy, 
                   comp_tile_info$tile_pheno_xy] <- NA 
            #set diagonals to NA
          }
          
        } 
        else{
          deltaF[this_tile_info$tile_pheno_xy, 
                 comp_tile_info$tile_pheno_xy] <- NA
        }
        comp_tile_counter <- comp_tile_counter + 1
      }
      this_tile_counter <- this_tile_counter + 1
    }  
    
    colnames(deltaF) <- index_tile_names$tile_pheno_xy
    rownames(deltaF) <- index_tile_names$tile_pheno_xy
    deltaF_list[[landscape_counter]] <- deltaF
    
  }
  return(deltaF_list)
  
}

make_Rs_integrated_old <- function(params, params_to_opt){
  
  R_maker <- function(X_list){
    w_names <- paste0("w", 1:length(X_list))
    function(w){
      ws <- c(w, 1-sum(w))
      names(ws) <- w_names
      w_dF <- lapply(1:length(ws), function(x) ws[x] * X_list[[x]])
      Ri <- Reduce('+', w_dF)
      return(Ri)
    }
  }
  
  weights <- matrix(params_to_opt$ws, nrow = length(params_to_opt$ws), 
                    ncol = 1, byrow = T)
  
  Rw <- R_maker(params$deltaF)
  Rs <- lapply(1:length(weights), function(x) Rw(weights[x,]))
  #Rs0 <- lapply(Rs, function(x) {x[is.na(x)] <- 0; x})
  
  return(Rs)
}

tile_phenospace_ids <- function(params){
  index_tile_names <- expand.grid(rep(seq(1:params$kx)), rep(seq(1:params$ky))) 
  colnames(index_tile_names) <- c("x_pheno", "y_pheno")
  index_tile_names$tile_pheno_num = seq(1:params$k2)
  index_tile_names$tile_pheno_xy = paste0(index_tile_names$x_pheno, "_", index_tile_names$y_pheno)
  
  return(index_tile_names)
}

tile_ecospace_ids <- function(params){
  pheno_tile_names <- expand.grid(rep(seq(1:params$kx)), rep(seq(1:params$ky)))
  colnames(pheno_tile_names) <- c("x_pheno", "y_pheno")
  pheno_tile_names$tile_pheno_num <- seq(1:params$k2)
  pheno_tile_names$tile_pheno_xy <- paste0(pheno_tile_names$x_pheno, "_", 
                                           pheno_tile_names$y_pheno)
  
  #now replicate it for 3 ecological groups
  eco_tile_names <- do.call("rbind", replicate(params$E, pheno_tile_names, simplify = F))
  eco_tile_names$tile_tot_num <- seq(1:(params$k2*params$E))
  eco_tile_names$eco_id <- rep(1:params$E, each = params$k2)
  
  eco_tile_names$x_eco <- pheno_tile_names$x_pheno + (params$kx*(eco_tile_names$eco_id - 1))
  eco_tile_names$y_eco <- pheno_tile_names$y_pheno + (params$ky*(eco_tile_names$eco_id - 1))
  
  eco_tile_names$tile_eco_xy <- paste0(eco_tile_names$x_eco, "_", 
                                       eco_tile_names$y_eco)
  
  return(eco_tile_names)
  
}

bin_data <- function(params){
  #where 'params' must already contain kerneled landscapes
  
  outside_bounds <- which(params$data[,1] > max(params$landscapes_list[[1]]$x) |
                        params$data[,1] < min(params$landscapes_list[[1]]$x) |
                        params$data[,2] > max(params$landscapes_list[[1]]$y) |
                        params$data[,2] < min(params$landscapes_list[[1]]$y) )
  if(length(outside_bounds) > 0){ params$data <- params$data[-c(outside_bounds), ]}
  
  params$data$raw_bin_x <- cut(params$data[,1], 
                               breaks = params$landscapes_list[[1]]$x, 
                               include.lowest = T,
                               labels = F)
  params$data$raw_bin_y <- cut(params$data[,2], 
                               breaks = params$landscapes_list[[1]]$y, 
                               include.lowest = T, 
                               labels = F)
    if(!is.numeric(params$data[,3])) {
      params$ecotype_key <- data.frame(eco = unique(params$data[,3]),
                                       numeric_labs = seq(1:length(unique(params$data[,3]))))
      params$data[,3] <- params$ecotype_key$numeric_labs[match(params$data[,3], params$ecotype_key$eco)]
    } 
  params$data$eco_bin_x <- params$data$raw_bin_x + (params$kx*(params$data[,3] - 1))
  params$data$eco_bin_y <- params$data$raw_bin_y + (params$ky*(params$data[,3] - 1))
  params$data$tile_eco_xy <- paste0(params$data$eco_bin_x, "_", params$data$eco_bin_y)
  
  params$data <- transform(merge(params$data, params$ecospace_tile_ids, 
                                 by = "tile_eco_xy", all.x = T, no.dups = T), 
                  row.names=row.names(params$data))
  
  #params$data <- merge(params$data, params$ecospace_tile_ids, by = "tile_eco_xy",
  #                     all.x = T, no.dups = T)
  
  return(params$data)
  
}

pre_binned_data_to_xy <- function(params){
  #where 'params' must already contain discretized landscapes
  #takes a set of data with single id bin assignments as gotten from simulating under a q
  #where the first column is the tile id and the second column is the eco id
  names(params$data) <- c("tile_pheno_num", "eco")
  
  if(all(params$data[1,] %in% params$phenospace_tile_ids$tile_pheno_num)){
    params$data <- merge(params$data, params$phenospace_tile_ids, 
                         by = "tile_pheno_num",
                         all.x = T, no.dups = T)
    
    params$data$eco_bin_x <- params$data$x_pheno + (params$kx*(params$data$eco - 1))
    params$data$eco_bin_y <- params$data$y_pheno + (params$ky*(params$data$eco - 1))
    params$data$tile_eco_xy <- paste0(params$data$eco_bin_x, "_", params$data$eco_bin_y)
    
    params$data <- merge(params$data, params$ecospace_tile_ids,
                         all.x = T, no.dups = T)
  } else{params$data <- merge(params$data, params$ecospace_tile_ids,
                              all.x = T, no.dups = T)}
  

  
  return(params$data)
  
}

scale_and_process_landscapes <- function(raw_landscapes_from_user, params){
  #check if landscape is pre-discretized
  if(params$landscape_format == "PCs"){
    L <- ncol(raw_landscapes_from_user) - 2
    x <- raw_landscapes_from_user[,1]
    y <- raw_landscapes_from_user[,2]

    #fix up x and y - set to n+1
    x_shift <- abs((unique(x)[1] - unique(x)[2])/2)
    y_shift <- abs((unique(y)[1] - unique(y)[2])/2)
    x <- c(rep((unique(x)[1] - x_shift), params$kx), x + x_shift)
    y <- c(rep((unique(y)[1] - y_shift), params$ky), y + y_shift)
    
    unscaled_kerneled <- list()
    #format the landscapes
    each_z = 3
    for(each_z in 3:(L+2)){
      processed_landscape <- list()
      processed_landscape$x <- unique(x)
      processed_landscape$y <- unique(y)
      processed_landscape$z <- matrix(data = raw_landscapes_from_user[,each_z], 
                                      nrow = (length(processed_landscape$y)-1),
                                      ncol = (length(processed_landscape$x)-1),
                                      byrow = T)
      unscaled_kerneled[[(each_z-2)]] <- processed_landscape
    }
  } else if(params$landscape_format == "distribution") {
    unscaled_kerneled <- mass_discretizer(raw_landscapes_from_user, params)
  } else {print("Landscape format not provided. Trying to discretize using MASS. Errors may occur.") 
    unscaled_kerneled <- mass_discretizer(raw_landscapes_from_user, params)}
  
    landscape_kerneled_list <- global_rescale_kerneled(kerneled_list_of_landscapes = 
                                                       unscaled_kerneled)
  return(landscape_kerneled_list)
}

# make_cache <- function(params){
#   #check if number of ecological groups exists
#   if(is.null(params$data) & is.null(params$E)){
#     print("No data or regimes found or provided. Assuming one regime.")
#     params$E <- 1
#   } else if (is.vector(params$data) & is.null(params$E)){
#     params$E <- params$E
#     print("No data or regimes found or provided. Assuming one regime.")
#   } else if (is.vector(params$data) & !is.null(params$E)){
#     params$E <- params$E
#     print(paste("Data does not contain regimes.", params$E, "regime(s) provided."))
#     } else if (!is.null(params$data) & is.null(params$E)){
#     params$E <- length(unique(params$data[,ncol(params$data)]))
#     print(paste("Data found with ", params$E, "regimes."))
#   } else if (!is.null(params$data) & !is.null(params$E)) {
#     params$E <- length(unique(params$data[,ncol(params$data)]))
#     print(paste("Data found with ", params$E,
#                 "regime(s)."))
#   } else {print(paste("No data found.", params$E,
#                       "regime(s) provided."))}
#   stopifnot(params$E > 0)
# 
#   #check landscape is digitized properly
#   if(params$landscape_format == "PCs"){
#     params$kx <- length(unique(params$raw_landscapes_from_user[,1]))
#     params$ky <- length(unique(params$raw_landscapes_from_user[,2]))
#     params$k2 <- params$kx * params$ky
#     if( params$kx == params$ky){
#       print(paste("Square landscape of dimensions", params$kx, "x", params$ky, "detected."))
#     } else {print(paste("Rectangular landscape of dimensions", params$kx, "x", params$ky, "detected."))}
#   } else if(!is.null(params$kx) & !is.null(params$ky)){
#     params$k2 <- params$kx * params$ky
#   } else if(!is.null(params$k)){
#     params$kx <- params$k
#     params$ky <- params$k
#     params$k2 <- params$kx * params$ky
#     print(paste("Square landscape dimensions", params$k, "x", params$ky, "provided"))
#   } else {stop("No landscape discretization dimensionality provided.
#           Please provide k, kx & ky, or a pre-discretized landscape.")}
# 
#   params$landscapes_list <- scale_and_process_landscapes(params$raw_landscapes_from_user, params)
#   #pull parameters from dataset
#   params$L <- length(params$landscapes_list) #number of landscapes
# 
#   params$n_free_ws <- (params$E*params$L) - params$E
#   #set up tile IDS
#   params$phenospace_tile_ids <- tile_phenospace_ids(params)
#   params$ecospace_tile_ids <- tile_ecospace_ids(params)
# 
#   #process data: continuous values
#   if(!is.null(params$data)){
#     if(length(params$data[[1]]) == 1){
#       #data is single tile binned with no regimes provided
#       params$binned_data <- Q_sim_to_xy_bins(params)
#       print("Single tile binned data processed and ecotypes extracted.")
#     } else if(length(params$data) == 2 & params$E > 1){
#       #data has to be single tile binned and regime provided
#       params$binned_data <- pre_binned_data_to_xy(params)
#       print("Single tile binned data and ecotypes processed.")
#     } else {if (length(params$data) == 3)
#       #ok, so we definitely have an x & y, but what are they
#       { if (all(is.wholenumber(params$data[,1:2]) & params$E == 1)){
#         params$binned_data <- bin_data(params) #need to write something to bin data that are xy binned already
#         print("Pre-binned data processed, maybe incorrectly. This is an error. Check this, Katherine.")
#       } else {params$binned_data <- bin_data(params)
#         print("Continuous data binned and provided ecotypes processed.")}}
#     }
#   } else (warning("No data found. Please provide data."))
# 
#     #did we have to remove any outside-the-bounds data? if so, tell people to remove those from tree
#   if(!is.null(params$binned_data)){
#     if(nrow(params$binned_data) != nrow(as.data.frame(params$data))){
#       print(paste("Data detected that is outside the bounds of the landscapes.",
#                   nrow(params$data) - nrow(params$binned_data),"species removed from dataset."))
#     }
#   }
# 
#   #take it out of the tree too
#   if(!is.null(params$binned_data)){
#     if(length(params$tree$tip.label) != nrow(params$binned_data)){
#       if(any(rownames(params$binned_data) %in% params$tree$tip.label)){
#         params$tree <- ladderize(drop.tip(params$tree,
#                                           base::setdiff(params$tree$tip.label, rownames(params$binned_data))))
#         print(paste("Data detected that is not in phylogeny.",
#                     nrow(params$data) - nrow(params$binned_data),"species removed from tree."))
#       } else {print("Data in phylogeny and dataset do not match. No rownames provided in dataset to rectify.
#                     You must fix this.")}
#     }
#   }
# 
#   #make deltaFs
#   params$deltaF <- make_deltaF_array(params)
# 
#   #set class
#   class(params) <- c("ripple","list")
# 
#   return(params)
# }

make_Rs_integrated <- function(params, params_to_opt){
  E <- params$E
  
  R_maker <- function(X_list){
    w_names <- paste0("w", 1:length(X_list))
    function(w){
      ws <- c(w, 1-sum(w))
      names(ws) <- w_names
      w_dF <- lapply(1:length(ws), function(x) ws[x] * X_list[[x]])
      Ri <- Reduce('+', w_dF)
      return(Ri)
    }
  }
  
  weights <- matrix(params_to_opt$ws, nrow = E, byrow = T)
  
  Rw <- R_maker(params$deltaF)
  Rs <- lapply(1:nrow(weights), function(x) Rw(weights[x,]))
  #Rs0 <- lapply(Rs, function(x) {x[is.na(x)] <- 0; x})
  
  return(Rs)
}

Q_amalgamator <- function(params, params_to_opt){
  #set up matrix of correct shape and size to catch amalgamators
  amalgamated_matrix <- matrix(nrow = ((params$k2)*params$E), 
                               ncol = ((params$k2)*params$E))
  
  R_counter = 1
  row_counter = 1
  for(R_counter in 1:params$E){
    
    amalgamated_matrix[row_counter:(R_counter*(params$k2)),
                       row_counter:(R_counter*(params$k2))] <- 
      params$Rs[[R_counter]]
    
    row_counter <- row_counter + params$k2
    #R_counter <- R_counter + 1
  }
  
  #Now integrate P matrices
  if(params$E > 1){
    amal_x = 1
    amal_y = 1
    
    matrix_shape <- matrix(seq(1:length(params_to_opt$P)), nrow = params$E, byrow = T)
    
    diag(params_to_opt$P) <- 0 #if this diagonal is not 0 you are in trouble
    for(amal_x in 1:(params$E)){
      
      for(amal_y in 1:(params$E)){
        
        diag(amalgamated_matrix[((params$k2)*(amal_x-1) + 1):((params$k2)*amal_x),
                                ((params$k2)*(amal_y-1) + 1):((params$k2)*amal_y)]) <- 
          params_to_opt$P[amal_x, amal_y]
        
        #amal_y <- amal_y + 1
      }
      #amal_x <- amal_x + 1
    }
  } else {print("One regime found. Ignoring any P.")}
  
  amalgamated_matrix[is.na(amalgamated_matrix)] <- 0
  diag(amalgamated_matrix) <- rowSums(amalgamated_matrix) * -1
  
  return(amalgamated_matrix)
}

apply_lm_to_R <- function(params, params_to_opt){
  
  lm_Rs <- list()
  R_counter <- 1
  
  for(R_counter in 1:length(params$Rs)){
    
    lm_Rs[[R_counter]] <- params$Rs[[R_counter]]*params_to_opt$p1s[[R_counter]] + params_to_opt$p0s[[R_counter]]
    #make sure that non-permitted transitions are still NA, because this will shield us from the
    #b0s getting added to the non-permitted transitions
    R_counter <- R_counter + 1
  }
  return(lm_Rs)
}

Rs_and_Q <- function(params, params_to_opt){ 
<<<<<<< HEAD
  #print(paste("length of params to opt is", length(params_to_opt$ws), 
  #            "length of (params$E*params$L)-params$E is", (params$E*params$L)-params$E))
  stopifnot(length(params_to_opt$ws) == ((params$E*params$L)-params$E))
=======
  print(paste("length of params_to_opt$ws is", length(params_to_opt$ws), 
              ", length of (params$E*params$L)-params$E is", (params$E*params$L)-params$E))
  print(paste("Is params_to_opt NULL?", is.null(params_to_opt)))
  
  #stopifnot(length(params_to_opt$ws) == ((params$E*params$L)-params$E))
>>>>>>> fc7230deb32e4998eadd76b0efdb1287177a2712
  #make p1s
  #params_to_opt$p1s <- exp(params_to_opt$b1s)#log(params_to_opt$b1s)
  #params_to_opt$p0s <- exp(params_to_opt$b0s)
  params_to_opt$p1s <- params_to_opt$b1s
  params_to_opt$p0s <- params_to_opt$b0s
  print("check")
  #make Rs
  params$Rs <-  make_Rs_integrated(params, params_to_opt)
  #transform Rs
  params$Rs <- apply_lm_to_R(params, params_to_opt)
  params$Rs <- lapply(params$Rs, exp)
  Q <- Q_amalgamator(params, params_to_opt)
  
  return_Rs_and_Q <- list()
  return_Rs_and_Q$Rs <- params$Rs
  return_Rs_and_Q$Q <- Q
  
  return(return_Rs_and_Q)
}

plot_Rs <- function(Rs_list, grid = T, arrows = T, arrow_scale = 1, color_palette = "#3498db", ...){
  
  k <- sqrt(dim(Rs_list[[1]])[1])
  #list_of_plots <- list()
  
  R_counter = 1
  for(R_counter in 1:length(Rs_list)){
    
    this_R <- Rs_list[[R_counter]]
    R_diag <- rowSums(this_R, na.rm = T)
    R_base_matrix <- matrix(R_diag, nrow = k, byrow = T)
    
    if (grepl("#", color_palette)){
      set_of_colors <- (sapply(seq(10, 255, length.out=(2.5*k)), 
                               function(x) makeTransparent(color_palette, x)))
    } else {set_of_colors <- hcl.colors((2.5*k), palette = color_palette, ...)}
    
    plot(0, type="n", xlim=c(0,1), ylim=c(0,1), 
         xlab="", ylab="", xaxt="n", yaxt="n") #someday want this to be real numebrs with real indices
    image(rotate_image(R_base_matrix), bty="n", 
          col = set_of_colors, add=TRUE)
    #add grid lines
    if(grid == T){
      u <- (1/k)/2
      abline(h=seq(-u, 1+u, length.out=k+1), lty=2)
      abline(v=seq(-u, 1+u, length.out=k+1), lty=2)
    }
    #add arrows
    if(arrows == T){
      centers <- expand.grid(seq(0,1, length.out=k), seq(1,0, length.out=k))
      index <- expand.grid(1:k, k:1)
      u <- c(0.02, 0.01)
      scf <- arrow_scale
      #rightwards
      for(i in 1:nrow(centers)){
        if((i)%%k!=0){
          arrows(centers[i,1]+u[1], centers[i,2]+u[2], centers[i+1,1]-u[1], centers[i+1,2]+u[2], 
                 length=0.025, lwd=this_R[i+1,i]^scf, 
                 col=makeTransparent("black", 100*this_R[i+1,i]/max(this_R,na.rm=TRUE)))
        }
      }
      #leftwards
      for(i in 1:nrow(centers)){
        if((i)%%k!=0){
          arrows(centers[i+1,1]-u[1], centers[i+1,2]-u[2], centers[i,1]+u[1], centers[i,2]-u[2], 
                 length=0.025, lwd=this_R[i,i+1]^scf, 
                 col=makeTransparent("black", 100*this_R[i+1,i]/max(this_R,na.rm=TRUE)))
        }
      }
      u <- c(0.01, 0.02)
      #downwards
      for(i in 1:nrow(centers)){
        if(i<(k*k-k)){
          arrows(centers[i,1]+u[1], centers[i,2]-u[2], centers[i+k,1]+u[1], centers[i+k,2]+u[2], 
                 length=0.025, lwd=this_R[i+k,i]^scf, 
                 col=makeTransparent("black", 100*this_R[i+k,i]/max(this_R,na.rm=TRUE)))
        }
      }
      #upwards
      for(i in 1:nrow(centers)){
        if(i<(k*k-k)){
          arrows(centers[i+k,1]-u[1], centers[i+k,2]+u[2], centers[i,1]-u[1], centers[i,2]-u[2], 
                 length=0.025, lwd=this_R[i,i+k]^scf, 
                 col=makeTransparent("black", 100*this_R[i,i+k]/max(this_R,na.rm=TRUE)))
        }
      }
    }
  }
}

plot_data <- function(params, pch = 16, alpha = 100, ...){
  axes_x <- round(params$landscapes_list[[1]]$x,2)
  axes_y <- round(params$landscapes_list[[1]]$y,2)
  
  plot(NULL, xlim=c(1,params$kx), 
       ylim=c(1,params$ky), ylab = "", xlab ="", xaxt = "n", yaxt = "n", ...)
  points(jitter(params$binned_data$x_pheno), jitter(params$binned_data$y_pheno), 
         pch = pch, col = makeTransparent("black", alpha = alpha))
  axis(1, at = seq(1, params$kx, length.out = params$kx+1), labels=axes_x) 
  axis(2, at = seq(1, params$ky, length.out = params$ky+1), labels=axes_y)
}

sim_data_from_landscapes <- function(sim_params, params_to_opt = NULL){
  catch_data <- list()
  
  stopifnot(length(params_to_opt$b1s) == length(sim_params$raw_landscapes_from_user))
  
  sim_params$L <- length(sim_params$raw_landscapes_from_user)
  sim_params$E <- nrow(sim_params$P)
  sim_params$k2 <- sim_params$kx * sim_params$ky
  sim_params$landscapes_list <- scale_and_kernel_landscapes(raw_landscapes_from_user, sim_params)
  sim_params$phenospace_tile_ids <- tile_phenospace_ids(sim_params)
  
  if(is.null(params_to_opt)){params_to_opt = list(b1s = c(rep(1, L)), b0s = c(rep(0, L)))
  } else{params_to_opt = params_to_opt}
  
  deltaFs <- make_deltaF_array(sim_params)
  
  count <- 1
  for(count in 1:sim_params$L){
    modify_deltaF <- ((exp(deltaFs[[count]]) * log(params_to_opt$b1s[count])) + params_to_opt$b0s[count])
    modify_deltaF[is.na(modify_deltaF)] <- 0
    diag(modify_deltaF) <- -1*rowSums(modify_deltaF)
    
    sim_data <- simulate_mk_model(tree = sim_params$tree, Q = modify_deltaF, 
                                  root_probabilities = "stationary", include_tips = T,
                                  include_nodes = T, Nsimulations = 1,
                                  drop_dims = T)$tip_states
    sim_regime <- simulate_mk_model(tree = sim_params$tree, Q = sim_params$P, 
                                    root_probabilities = "stationary", include_tips = T,
                                    include_nodes = T, Nsimulations = 1,
                                    drop_dims = T)$tip_states
    data <- data.frame(sim_data, sim_regime)
    colnames(data) <- c("sim_data", "sim_regime")
    
    catch_data[[count]] <- list(modify_deltaF, data)
    names(catch_data[[count]]) <- c("transition_matrix", "sim_data")
  }
  return(catch_data)
}

sim_data_from_Q <- function(no_data_params){
  sim_data <- simulate_mk_model(tree = no_data_params$tree, Q = no_data_params$Q, 
                                root_probabilities = "stationary", include_tips = T,
                                include_nodes = F, Nsimulations = 1,
                                drop_dims = T)$tip_states
  return(sim_data)
}

Q_sim_to_xy_bins <- function(params){
  
  make_breaks <- c(1)
  break_counter = 1
  for(break_counter in 1:params$E){
    this_break <- params$k2*break_counter
    make_breaks <- c(make_breaks, this_break)
  }
  
  sim_regime <- cut(params$data, breaks = make_breaks, labels = F,
                    right = T, include.lowest = T)
  
  params$data <- data.frame(params$data, sim_regime)
  names(params$data) <- c("tile_tot_num", "eco")
  
  params$data <- merge(params$data, params$ecospace_tile_ids,
                      by = "tile_tot_num",
                      all.x = T, no.dups = T)
  return(params$data)
}

plot.ripple <- function(params, opt_params, color_palette = "#3498db", 
                        legend = T, grid = T, arrows = F, plot_data = F, 
                        data_pch = 16, data_col = "black", data_alpha = 200,
                        arrow_scale = 0.3, ...){
  axes_x <- round(params$landscapes_list[[1]]$x,2)
  axes_y <- round(params$landscapes_list[[1]]$y,2)
  
  if(!is.null(params$Rs)){
    Rs_list <- params$Rs
  } else {
    Rs_list <- (Rs_and_Q(params, opt_params))$Rs
  }
  
  R_counter = 1
  for(R_counter in 1:length(Rs_list)){
    
    this_R <- Rs_list[[R_counter]]
    R_diag <- rowSums(this_R, na.rm = T)
    R_base_matrix <- matrix(R_diag, nrow = params$ky, ncol = params$kx, byrow = T)
    
    if (grepl("#", color_palette)){
      set_of_colors <- (sapply(seq(10, 255, length.out=(2.5*params$kx)), 
                               function(x) makeTransparent(color_palette, x)))
    } else {set_of_colors <- hcl.colors((2.5*params$kx), palette = color_palette, ...)}
    
    #add plot legend and set up margins
    if(legend == T){
      par(mar=c(5, 4, 4, 8), xpd=T)
    } else {par(mar=c(5, 4, 4, 2), xpd=F)}
    
    #add plot
    plot(0, type="n", xlim=c(0,1), ylim=c(0,1), 
         xlab="", ylab="", xaxt="n", yaxt="n") 
    image(rotate_image(R_base_matrix), bty="n", 
          col = set_of_colors, add=TRUE)
    
    #add axis labels
    axis(1, at = seq(0, 1, length.out = params$kx+1), labels=axes_x) 
    axis(2, at = seq(0, 1, length.out = params$ky+1), labels=axes_y)
    
    #add legend
    if(legend == T){
      legend("topright", inset=c(-.25,0), bty = "n",
             legend=c("high","low"), pch=15, xpd = T,
             title="Performance\nSurface\nHeight", pt.cex = 3, pt.lwd = 1,
             x.intersp = 2, y.intersp = 2,
             col = c(set_of_colors[length(set_of_colors)], set_of_colors[1]))
      legend("topright", inset=c(-.25,0), bty = "n",
             legend=c("high","low"), pch=0, xpd = T,
             title="Performance\nSurface\nHeight", pt.cex = 3, pt.lwd = 1,
             x.intersp = 2, y.intersp = 2,
             col = "black")
    }
    
    #add points for data
    if(plot_data == T){
      if(is.null(params$data)){
        print("No data found, data cannot be plotted.")
      } else{points(jitter(rescale_0to1(params$binned_data$x_pheno)), 
                    jitter(rescale_0to1(params$binned_data$y_pheno)), 
                    pch = data_pch, 
                    col = makeTransparent(data_col, alpha = data_alpha))}
    }
    
    #add grid lines
    if(grid == T){
      ux <- (1/params$kx)/2
      uy <- (1/params$ky)/2
      abline(h=seq(-u, 1+u, length.out=params$kx+1), lty=2, xpd = F)
      abline(v=seq(-u, 1+u, length.out=params$ky+1), lty=2, xpd = F)
    }
    
    if(arrows == T){
      centers <- expand.grid(seq(0,1, length.out=params$kx), 
                             seq(1,0, length.out=params$ky))
      index <- expand.grid(1:params$kx, params$ky:1)
      u <- c(0.02, 0.01)
      scf <- arrow_scale
      #rightwards
      for(i in 1:nrow(centers)){
        if((i)%%params$kx!=0){
          arrows(centers[i,1]+u[1], centers[i,2]+u[2], centers[i+1,1]-u[1], centers[i+1,2]+u[2], 
                 length=0.025, lwd=this_R[i+1,i]^scf, 
                 col=makeTransparent("black", 100*this_R[i+1,i]/max(this_R,na.rm=TRUE)))
        }
      }
      #leftwards
      for(i in 1:nrow(centers)){
        if((i)%%params$kx!=0){
          arrows(centers[i+1,1]-u[1], centers[i+1,2]-u[2], centers[i,1]+u[1], centers[i,2]-u[2], 
                 length=0.025, lwd=this_R[i,i+1]^scf, 
                 col=makeTransparent("black", 100*this_R[i+1,i]/max(this_R,na.rm=TRUE)))
        }
      }
      u <- c(0.01, 0.02)
      #downwards
      for(i in 1:nrow(centers)){
        if(i<(params$ky*params$ky-params$ky)){
          arrows(centers[i,1]+u[1], centers[i,2]-u[2], centers[i+k,1]+u[1], centers[i+k,2]+u[2], 
                 length=0.025, lwd=this_R[i+k,i]^scf, 
                 col=makeTransparent("black", 100*this_R[i+k,i]/max(this_R,na.rm=TRUE)))
        }
      }
      #upwards
      for(i in 1:nrow(centers)){
        if(i<(params$ky*params$ky-params$ky)){
          arrows(centers[i+k,1]-u[1], centers[i+k,2]+u[2], centers[i,1]-u[1], centers[i,2]-u[2], 
                 length=0.025, lwd=this_R[i,i+k]^scf, 
                 col=makeTransparent("black", 100*this_R[i,i+k]/max(this_R,na.rm=TRUE)))
        }
      }
    }
  }
}

# get_opt_params <- function(params, params_to_opt){
#   #assuming params_to_opt is a labeled list
#   
#   if(is.null(params_to_opt)){params_to_opt <- params$params_to_opt}
#   
#   if(params$E ==1){
#     params_to_opt$P <- matrix(c(1, 1, 1, 1), nrow = 2, ncol = 2)}
#   
#   #check number of parameters is appropriate
#   if(length(params_to_opt$b1s) != params$E) {
#     stop(paste("Incorrect number of b1s provided.", params$E, "b1s expected, but", length(params_to_opt$b1s), "provided."))
#   }
#   
#   if(length(params_to_opt$b0s) != params$E) {
#     stop(paste("Incorrect number of b0s provided.", params$E, "b0s expected, but", length(params_to_opt$b0s), "provided."))
#   }
#   
#   if(length(params_to_opt$b0s) != params$E) {
#     stop(paste0("Incorrect number of b0s provided.", params$E, "b0s expected, but", length(params_to_opt$b0s), "provided."))
#   }
#   
#   if(length(params_to_opt$ws) == (params$E*params$L)){
#     #get out each set of ws and make sure they are the right length
#     this_R <- 1
#     w_counter <- 1
#     tot_ws <- params$E*params$L
#     tot_ws_per_deltaF <- tot_ws/params$E
#     free_ws_per_deltaF <- tot_ws_per_deltaF - 1
#     free_ws <- numeric(0) #set up to catch
#     for(this_R in 1:params$E){
#       these_ws <- params_to_opt$ws[(w_counter):((w_counter) + tot_ws_per_deltaF-1)]
#       stopifnot(sum(these_ws) == 1)
#       these_free_ws <- these_ws[1:free_ws_per_deltaF]
#       free_ws <- c(free_ws, these_free_ws)
#       w_counter <- w_counter + tot_ws_per_deltaF
#     }
#     
#     #pull weights here and format like we want them -- allow users to submit all the weights
#   } else if(length(params_to_opt$ws) == ((params$E*params$L)-params$E)){
#     print("Correct number of weights provided.")
#   } else{ stop("Incorrect number of weights provided.")}
#   #PUT IN TEST TO MAKE SURE WEIGHTS ARE SUMMING TO 1
#   
#   #get ws into a format we can work with for the transformation
#   if(any(params_to_opt$ws == 1)){
#     params_to_opt$ws <- replace(x <- c(params_to_opt$ws), x==1, 0.999)}
#   
#   if(any(params_to_opt$ws == 0)){
#     params_to_opt$ws <- replace(x <- c(params_to_opt$ws), x==0, 0.001)}
#   #gotta back transform ws into x
#   params_to_opt$xs <- sapply(params_to_opt$ws, make_x) #THIS HERE
#   #transformation can't handle any number that's actually larger than or equal to 1 or less than or equal to 0
#   
#   param_vec <- c(params_to_opt$b1s, params_to_opt$b0s, params_to_opt$xs, 
#                  params_to_opt$P)
#   
#   n_ws <- (params$E*params$L) - params$E
#   
#   optimizing <- optimx(param_vec, params = params, fn = likfn, 
#                        method = "Nelder-Mead")
#   
#   cache_out <- list()
#   cache_out$opt_params$b1s <- (optimizing[1:params$E]) #so that this function returns b1, not p1
#   cache_out$opt_params$b0s <- (optimizing[(params$E+1):(params$E*2)])
#   xs <- optimizing[((params$E*2)+1):((params$E*2)+n_ws)] #E*L - E
#   cache_out$opt_params$ws <- sapply(xs, make_w)  #w = exp(x) / (1 + exp(x))
#   if(params$E > 1){
#     cache_out$opt_params$P <- matrix(optimizing[((params$E*2)+n_ws+1):length(param_vec)],
#                                      nrow = params$E, ncol = params$E, byrow = T)
#   } else {params_out$opt_params$P <- 1}
#   cache_out$opt_params$lnl <- optimizing$value*-1
#   
#   #get likelihood of optimized parameters
#   # cache_out$lnl <- likfn(param_vec = c(cache_out$opt_params$b1s, 
#   #                                      cache_out$opt_params$b0s, 
#   #                                      unlist(xs),
#   #                                      cache_out$opt_params$P), params)*-1
#   #make q matrix from optimized parameters
#   
#   transition_matrices <- Rs_and_Q(params = params, 
#                                    params_to_opt = params_out$opt_params)
#   cache_out$transition_matrices <- list()
#   cache_out$transition_matrices$Rs <- transition_matrices$Rs
#   cache_out$transition_matrices$Rs <- transition_matrices$Rs
#   
#   cache_out$n_free_params <- (params$E*2) + n_ws + ((params$E^2) - params$E)
#   cache_out$opt_params$AIC <- (2*cache_out$n_free_params) - (2*cache_out$lnl)
#   k <- length(params$tree$tip.label)
#   cache_out$opt_params$AICc <- cache_out$opt_params$AIC -((2*k*(k+1))/(cache_out$n_free_params-k-1))
#   cache_out$input_pars <- params
#   
#   class(cache_out) <- c("list", "ripple")
#   
#   return(cache_out)
# }

likfn <- function(param_vec, params){
  
  #get lengths of everything
  opt_pars <- list()
  
  opt_pars$b1s <- param_vec[1:params$E]
  opt_pars$b0s <- param_vec[(params$E+1):(params$E*2)]
  n_ws <- (params$E*params$L) - params$E
  opt_pars$ws <- param_vec[((params$E*2)+1):((params$E*2)+n_ws)] #E*L - E
  if(params$E > 1){
    opt_pars$P <- matrix(param_vec[((params$E*2)+n_ws+1):length(param_vec)],
                         nrow = params$E, ncol = params$E, byrow = T)}
  
  make_Rs_and_Q <- Rs_and_Q(params, opt_pars)
  params$Rs <- make_Rs_and_Q$Rs
  params$Q <- make_Rs_and_Q$Q
  
  res <- castor::asr_mk_model(tree = params$tree, 
                              tip_states = params$binned_data$tile_tot_num,
                              Nstates = (params$k2 * params$E),
                              transition_matrix = params$Q,
                              include_ancestral_likelihoods = F)
  
  return(res$loglikelihood*-1)
  
}





