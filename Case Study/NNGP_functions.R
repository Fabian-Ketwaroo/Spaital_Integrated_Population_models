
library(nimble)


  
  # Exponential covariance matrix 
  expcov <- nimbleFunction(     
    run = function(dists = double(2), rho = double(0), sigma2 = double(0)) {
      returnType(double(2))
      n <- dim(dists)[1]
      y <- dim(dists)[2]
      result <- matrix(nrow = n, ncol = y, init = FALSE)
      nugget <- diag(n)*1e-06 # add nugget effect for numerical stability
      for(i in 1:n)
        for(j in 1:y)
          result[i, j] <- sigma2*exp(-dists[i,j]/rho)
      return(result + nugget)
    })
  
  # rho: Length-scale: Controls how quickly the covariance decays with distance.
  # sigma2: Signal Variance: Scales the overall variation of the function.
  
  
  # computeNeighbours takes a set of coordinates and determines the neighbor set using Vecchia's approximation
  # Ordering is done on the x-axis, i.e, order(coords[,1])
  computeNeighbours <- function(coords,k , proj){
    
    N <- dim(coords)[1]
    # Sort coordinates by x axis 
    sorted_idx <- order(coords[,1]) # in ascending order
    coords_sorted =coords[sorted_idx, ]
    
    # Compute euclidean distance matrix 
    coords_sorted.sf <- st_as_sf(data.frame(x = coords_sorted[, 1], y = coords_sorted[, 2]),
                                 coords = c('x', 'y'), crs = proj)
    
    
    edist = matrix(as.numeric(st_distance(coords_sorted.sf)/1000), N, N)
    #edist = as.matrix(dist(coords_sorted, diag = TRUE, upper = TRUE))
    #edist
    
    neighbors <- matrix(0,N,N) # 1 if elements are neighbors
    neighbor_idx <- matrix(0, N, k)
    neighbors_dist <- matrix(0, N, k+1) # the k+1 is the distance between the point and it self. Diag edist
    neighbors_dist[ 1, 1:k ] = edist[1, 1:k]
    
    # the 1st location has no neighbors and then k+1 location has k neighbours to the left
    for (i in  2:(k + 1)) {
      neighbor_idx[i, 1:(i - 1) ] = as.numeric(1:(i -  1))
      neighbors[i, neighbor_idx] = 1
      neighbors_dist[i,  1:(i - 1)  ] <- edist[i, neighbor_idx[i, 1:(i - 1) ]]
      neighbors_dist[ i, i:k] <- edist[i, i:k]
    }
    
    for (i in (k+2):N) {
      curr_sorted_idx = order(edist[i, ])
      curr_sorted_idx = curr_sorted_idx[curr_sorted_idx < i] # take values to the left
      neighbor_idx[i,  ] = curr_sorted_idx[1:min(i, k)]
      neighbors[i, neighbor_idx[i,]] = 1
      neighbors_dist[i, 1:k ] <- edist[i, neighbor_idx[i,  ] ]
    }
    
    
    return(list(coords_sorted = coords_sorted, edist_sorted = edist, neighbors = neighbors, neighbor_idx = neighbor_idx, neighbors_dist = neighbors_dist))
    
  }
  
  DisplayNeighbors <- function( ind, coords.sorted, neighbor_idx, proj ){
    
    coords <-coords.sorted
    
    # Coordinates for plotting
    coords.sf <- st_as_sf(data.frame(x = coords[, 1], y = coords[, 2]),
                          coords = c('x', 'y'), crs = proj)
    
    
    if(ind ==2) nlocation = st_as_sf(data.frame(t(coords[neighbor_idx[ind,],])), coords = colnames( coords)[1:2], crs = proj) 
    else nlocation = st_as_sf(data.frame(coords[neighbor_idx[ind,],]), coords = colnames( coords)[1:2], crs = proj) 
    #nlocation = data.frame(coords[neighbor_idx[ind,],])
    
    
    result <-  ggplot() +
      geom_sf(data = coords.sf, pch = 2) +
      theme_bw(base_size = 14) +
      scale_color_viridis_d() +
      labs(x = 'x', y = 'y') +
      geom_sf(data = nlocation, col = "blue")+
      geom_point(aes(x = coords[ind,1], y = coords[ind,2]), color = "red", size = 3)+ 
      geom_vline(xintercept =coords[ind,1], linetype="dotted", linewidth =1)
    
    
    return(result)
    
  }
  
  computeC <- nimbleFunction(
    run = function(edist = double(2), nid.dist = double(1), rho = double(0), sigma2 = double(0) ){
      returnType(double(2))
      
      k1 <- length(nid.dist)-1
      N <- dim(edist)[1]
      C <- matrix(0,k1+1, k1+1 ) 
      C[1:k1, 1:k1] <- expcov(nimMatrix(edist,N,N)  , rho, sigma2) # C value between neighbors of i 
      C[k1+1, 1:(k1+1)] <- sigma2*exp(-nid.dist/rho) # C value between i and its neighbors
      C[ k1+1, k1+1 ] <- C[ k1+1, k1+1 ] + 1e-06 # add nugget for computational stability
      # the last line here can be coded better. Something to come back to.
      
      return(C)
      
      
    }
  )
  
  computeAD <- nimbleFunction(
    run = function( edist = double(2), nid.dist = double(2), neighbors.id = double(2), rho = double(0), sigma2 = double(0), k = double(0)   ){
      returnType(double(2))
      
      # calculate A and D
      N <- dim(edist)[1]
      AD  <- matrix(0, N, k+1)
      AD[1, k+1] = sigma2
      
      for (i in 1:(N-1)) {
        
        if(i < k) {
          
          
          k1 <- i
          
          # Select neighbors of i
          n.id <- neighbors.id[i+1, 1:k1 ]
          bdist <- edist[n.id, n.id] # distance between neighbors of i
          ndist <- nid.dist[i+1, 1:(k1+1)] # distance between i and its neighbors
          C <- computeC(bdist, ndist , rho, sigma2 )
          
          AD[i+1,  1:k1 ] = solve(C[1:k1, 1:k1], C[ i+1, 1:k1])
          AD[i+1, k+1] =  C[ i+1, k1+1] - inprod(C[i+1, 1:k1], AD[i+1, 1:k1])
          
          
        } else{
          
          k1 <- k
          
          # Select neighbors of i
          n.id <- neighbors.id[i+1, 1:k1 ]
          
          bdist <- edist[n.id, n.id] # distance between neighbors of i
          ndist <- nid.dist[i+1, 1:(k1+1)] # distance between i and its neighbors
          C <- computeC(bdist, ndist , rho, sigma2 )
          
          
          AD[i+1,  1:k1 ] = solve(C[1:k1, 1:k1], C[ k1+1, 1:k1])
          AD[i+1, k+1] =  C[ k1+1, k1+1] - inprod(C[k1+1, 1:k1], AD[i+1, 1:k1])
          
        } 
      }
      
      return(AD)
      
    }
  )
  
  
  computeQF <- nimbleFunction(
    run = function(u = double(1), v = double(1), AD = double(2), neighbors.id = double(2) ){
      returnType(double(0))
      
      N <- dim(AD)[1]
      k <- dim(AD)[2]-1
      
      qf = u[1]*v[1]/AD[1,k+1]
      
      for(i in 2:N) { 
        if(i <= k) k1 <- i-1 else k1 <- k
        qf = qf + (u[i] - inprod(AD[i,1:k1], u[neighbors.id[i,1:k1]]))*(v[i] - inprod(AD[i, 1:k1 ], v[neighbors.id[i,1:k1] ]))/AD[i,k+1] 
      }
      
      return(qf)
      
    }
  )
  
  dmnorm_NNGP <- nimbleFunction(
    run = function( x= double(1), mu = double(1), AD = double(2), neighbors.id = double(2), log = integer(0, default =1) ){
      returnType(double(0))
      
      
      N = dim(AD)[1]
      M = dim(AD)[2]
      u = x- mu
      v = x - mu
      qf <- computeQF(u,v, AD, neighbors.id)
      
      loglik <- -0.5*(N*1.83787706649  +  sum(log(AD[1:N, M])) + qf )   # log(2pi) = 1.8378770664
      
      if(log)return(loglik)
      return(exp(loglik))
      
    }
  ) 
  
  

