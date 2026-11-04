
# Load packages
{
  #rm(list = ls())
  library(nimble)
  library(sf)
  library(MCMCvis)
  library(coda)
  library(MASS)
  library(parallel)
}

# NIMBLE functions
{
  

  expcov <- nimbleFunction(
    run = function(dists = double(2), rho = double(0), sigma2 = double(0)) {
      returnType(double(2))
      n <- dim(dists)[1]
      y <- dim(dists)[2]
      result <- matrix(nrow = n, ncol = y, init = FALSE)
      for(i in 1:n)
        for(j in 1:y)
          result[i, j] <- sigma2*exp(-dists[i,j]/rho)
      return(result)
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
    
    
    edist = matrix(as.numeric(st_distance(coords_sorted.sf)), N, N)
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
  
  rmnorm_NNGP<- nimbleFunction(
    run = function(n = integer(0),  mu = double(1), AD = double(2), neighbors.id = double(2)) {
      stop('rmnorm_NNGP should never be run')
      returnType(double(1))
      return(0)
    })
  
  
  
  
  NNGP.pred <- function(neighbors_dist.x0, neighbors_dist.N0, x0.location, rho, sigma2, w.s.x0) {
    
    c = expcov(nimMatrix(neighbors_dist.x0), rho = rho, sigma2 = sigma2) # covariance between x0 and its k nearest neighbours
    if( any(x0.location== c(1,2) ) )  C_N0 = expcov(nimMatrix(neighbors_dist.N0) , rho = rho, sigma2 = sigma2) # covariance betwwen the k nearest neighbours
    else C_N0 = expcov(neighbors_dist.N0 , rho = rho, sigma2 = sigma2) + (diag( dim(neighbors_dist.N0 )[1] )*1e-06) # covariance betwwen the k nearest neighbours
    C_x0 = expcov(nimMatrix(0), rho = rho, sigma2 = sigma2) # covariance at prediction point x0
    
    
    m = inprod(c, solve(C_N0 , w.s.x0) ) # add nuggest to aviod singularity
    v = C_x0 - inprod(c, solve(C_N0, c))
    a = rnorm(1, 0, 1)
    gp.pred = m + sqrt(v)*a
    return(gp.pred)
  }
  
  coverage <- function(est, true){
    
    true = c(true)
    N <- dim(est)[1]
    Ncov <- numeric(N)
    for (i in 1:N) {
      if(est[i,3] <= true[i] & est[i,5] >= true[i]) Ncov[i] =1
      else Ncov[i] = 0
    }
    
    return(mean(Ncov))
    
  }
  
  
}


# Settings
NC = 400
NCR = 100
N = NC + NCR
P = 100 # number of prediction points
M <- N + P
k = 15 # number of neighbors
n.occasions= 6# 5 years of data

nsim = 100 
numbers = 1:nsim

# Objects to save outputs
key_parameters = array(0, dim = c(29,7,nsim))
colnames(key_parameters) = c("mean",  "sd",    "2.5%",  "50%",   "97.5%", "Rhat",  "n.eff")
rownames(key_parameters) = c(paste0("alpha.p[", 1:(n.occasions-1), "]" ), "alpha2.af" , "alpha2.as", paste0("beta.f[", 1:5, "]" ), paste0("beta.s[", 1:5, "]" ) 
                             ,"rho.af",  "rho.as", paste0("wft[",1:(n.occasions-1), "]" ), paste0("wst[",1:(n.occasions-1), "]" ) )

w.as.est = w.af.est = matrix(0, nsim , 5)
colnames(w.as.est) =   colnames(w.af.est)=  c("Mean RB", "2.5% RB", "97.5% RB", "RMSE", "Cov")

phi_est =  matrix(0, nsim , 10)
colnames(phi_est) =   c("Mean RB", "2.5% RB", "97.5% RB", "RMSE", "Cov", paste0("Year ", 1:(n.occasions-1), " RB" ))

rec_est =   matrix(0, nsim , 10)
colnames(rec_est) = c("Mean RB", "2.5% RB", "97.5% RB", "RMSE", "Cov", paste0("Year ", 1:(n.occasions-1), " RB" ))

N_est = matrix(0, nsim , 11)
colnames(N_est) = c("Mean RB", "2.5% RB", "97.5% RB", "RMSE", "Cov", paste0("Year ", 1:(n.occasions), " RB" ))

gr_est = matrix(0, nsim , 10)
colnames(gr_est) = c("Mean RB", "2.5% RB", "97.5% RB", "RMSE", "Cov", paste0("Year ", 1:(n.occasions-1), " RB" ))

Train_time = matrix(0, nsim, 3)
colnames(Train_time) = c("user",  "system",  "elapsed")


phi.predict.est <- rect.predict.est <-   array(0, dim = c(nsim, 6, n.occasions-1))#matrix(0, nsim , 6)
colnames(phi.predict.est) = colnames(rect.predict.est) =   c("Mean RB", "2.5% RB", "97.5% RB", "MAPE", "RMSE", "Cov")  



for (n in 1:nsim) {
  
  set.seed(numbers[n])
  
 
  # Simulate data 
  coords = matrix( 0, M, 2 )
  coords[,1] = runif(M,0,1)
  coords[,2] = runif(M,0, 1)
  
  
  coords <- data.frame(x = coords[, 1], y = coords[, 2], id = 1:M)
  source = as.factor(c(rep("C", NC), rep("CR", NCR), rep("Pred", P) ))
  coords$source = source
  
  # Coordinates for plotting
  coords.sf <- st_as_sf(coords,
                        coords = c('x', 'y'))
  
  
  edist = matrix(as.numeric(st_distance(coords.sf)), M, M)
  
  # Covariates 
  TMAX.rho = 0.4
  TMAX.alpha2 = 0.2
  X_tmax = matrix(0, M, n.occasions-1)
  
  for (t in 1:(n.occasions-1)) {
    X_tmax[ , t]= MASS::mvrnorm(1, mu =rep(0,M), Sigma = expcov(edist, TMAX.rho, TMAX.alpha2))
  }
  
  
  X_long <- matrix(0, nrow = M, ncol = n.occasions-1)
  X_long[, 1:(n.occasions-1)] <- coords$x
  
  X_lat <- matrix(0, nrow = M, ncol = n.occasions-1)
  X_lat[, 1:(n.occasions-1)] <- coords$y
  
  X_tmax_long <- X_tmax * X_long
  X_tmax_lat <- X_tmax * X_lat
  
  
  # Now let's scale values 
  X_tmax <- scale(X_tmax)
  X_long <- scale(X_long)
  X_lat <- scale(X_lat)
  X_tmax_lat <- scale(X_tmax_lat)
  X_tmax_long <- scale(X_tmax_long)
  ncov = 5
  X.cov = array(0, dim = c(M, n.occasions-1, ncov))
  X.cov[,,1] = X_tmax
  X.cov[,,2] = X_long
  X.cov[,,3] = X_lat
  X.cov[,,4] = X_tmax_long
  X.cov[,,5] = X_tmax_lat
  
  colnames(X.cov) = c("TMAX", "Long", "Lat", "TMAXlong", "TMAXlat")
  
  rel.R <- matrix(0, nrow = M, ncol = n.occasions-1)
  
  # CR data only at NCR sites 
  for (i in 1:N ) {
    rel.R[i, ] <-sample(20:50, n.occasions-1, replace = FALSE) 
  }
  
  visits.all <- matrix(0, M, n.occasions-1)
  #cr.level = 0.8 # 4/5
  cr.level = 0.6 # 3/5
  for (i in 1:NCR) {
    n.ones <- (n.occasions-1)*cr.level
    id <- sample(1:(n.occasions-1), n.ones)
    visits.all[((NC+1):N)[i] , id ] <- 1
  }
  
  win.data = list(rel.R= rel.R, n.occasions=n.occasions, M=M,  mu = rep(0, M), edist = edist, dmax = max(edist), dmin = min(edist[edist!= 0] ), dq75 = quantile(edist)[4], ncov = ncov, visits = visits.all)
  str(win.data)
  
  sipm <- nimbleCode({
    
    
    for (i in 1:ncov) {
      beta.f[i] ~ dnorm(0, var = 2) 
      beta.s[i] ~ dnorm(0, var = 2) 
    }
    
    
    
    ################### 
    # GP's
    ###################
    
    # Recruitment
    alpha2.af ~ dinvgamma(3,1)
    rho.af ~ dunif(dmin/3, dmax/3)
    
    cov.af[1:(M), 1:(M)] <- expcov(edist[1:(M), 1:(M)], rho.af, alpha2.af)
    w.af.star[1:(M)] ~ dmnorm(mu[1:(M)], cov = cov.af[1:(M), 1:(M)])
    w.af[1:(M)] <- w.af.star[1:(M)] - mean(w.af.star[1:(M)])
    
    
    
    # Survival
    alpha2.as ~ dinvgamma(3,1)
    rho.as ~ dunif(dmin/3, dmax/3)
    
    cov.as[1:(M), 1:(M)] <- expcov(edist[1:(M), 1:(M)], rho.as, alpha2.as)
    w.as.star[1:(M)] ~ dmnorm(mu[1:(M)], cov = cov.as[1:(M), 1:(M)])
    w.as[1:(M)] <-  w.as.star[1:(M)] - mean(w.as.star[1:(M)])
    
    
    for (t in 1:(n.occasions-1)) {
      wst[t] ~ dnorm(0, var = 2)
      wft[t] ~ dnorm(0, var = 2)
    }
    
    
    for (t in 1:(n.occasions-1)) {
      alpha.p[t] ~ dnorm(0, var = 2)
    }
    
    
    for (i in 1:M) {
      
      # Priors and linear models
      for (t in 1:(n.occasions-1)){
        
        logit(sa[i,t]) <- wst[t] + (Xcov[i,t, 1:ncov] %*% beta.s[1:ncov])[1,1]  + w.as[i] 
        log(f[i,t]) <- wft[t] + (Xcov[i,t, 1:ncov] %*% beta.f[1:ncov])[1,1] + w.af[i] 
        p[i,t] <- ilogit(alpha.p[t])* visits[i,t] 
      }
      
    }
    
    
    for (i in 1:M) {
      
      u[i] ~ dunif(50,150)
      N[i,1] <- round(u[i])
      
      # Process model over time: our model of population dynamics
      for (t in 1:(n.occasions-1)){
        
        S[i,t+1] ~ dbin(sa[i,t], N[i,t])
        R[i,t+1] ~ dpois(f[i,t]*N[i,t])
        N[i,t+1] <- S[i, t+1] + R[i, t+1] 
        
        # Growth rate
        growth.rate[i,t] <- (N[i, t+1]) / (N[i,t] + 0.001)
      }
      
      
      for (t in 1:n.occasions){
        log(lambda[i,t]) <- log(N[i,t])
        C[i, t] ~ dpois(lambda [i,t])
      }
      
      
      # Capture-recapture data (CJS model with multinomial likelihood)
      for (t in 1:(n.occasions-1)){
        Rm[i,t,1:n.occasions] ~ dmulti(r[i,t,1:n.occasions], rel.R[i,t])
      }
      
      
      # Define the cell probabilities of the m-array
      # Main diagonal
      for (t in 1:(n.occasions-1)){
        q[i, t] <- 1-p[i, t]                # Probability of non-recapture
        r[i, t,t] <- sa[i,t]*p[i,t]
      }
      # Above main diagonal
      for (t in 1:(n.occasions-2)){
        for (j in (t+1):(n.occasions-1)){
          r[i,t,j] <- prod(sa[i,t:j])*prod(q[i,t:(j-1)])*p[i,j]
        } #j
      } #t
      # Below main diagonal
      for (t in 2:(n.occasions-1)){
        for (j in 1:(t-1)){
          r[i,t,j] <- 0
        } #j
      } #t
      # Last column: probability of non-recapture
      for (t in 1:(n.occasions-1)){
        r[i,t,n.occasions] <- 1-sum(r[i,t,1:(n.occasions-1)])
      } #t
      
      
      
    }
    
    
  })
  
  simModel = nimbleModel(sipm, constants = win.data, data = list(Xcov = X.cov), calculate = TRUE)# Here calculate needs to be set to true so that the lower dia of pr.a is set to zero not NA
  # The not positive definite error is okay here as I haven't passed any initial values 
  
  
  nodesToSim <- simModel$getDependencies(c( "alpha.p", "rho.af", "alpha2.af", "rho.as",  "alpha2.as", "wft","wst", "u",
                                            "beta.s", "beta.f" ), self = F, downstream = T)
  
  
  # Values used to simulate data
  set.seed(numbers[n])
  simModel$u = runif(1*win.data$M, 80, 120)
  simModel$alpha.p = logit(runif(n.occasions-1, 0.4, 0.7))
  simModel$beta.f = c(-0.2, 0.1, -0.4, 0.2, 0.2)  
  simModel$beta.s = c(-0.1, 0.1, 0.5, -0.2, -0.5)
  
  simModel$rho.af = 0.35
  simModel$alpha2.af = 0.1
  simModel$rho.as = 0.3
  simModel$alpha2.as = 0.2
  
  # Here, I want the growth rate to be around 1 
  # If gr =1, the f = 1- survival
  simModel$wst = logit(runif(n.occasions-1, 0.5,0.85))
  simModel$wft = log(1 -ilogit(simModel$wst)) 
  
  simModel$simulate(nodesToSim)
  #if (length(which(simModel$growth.rate == 0)) >0 ){ test[n] <- which(simModel$growth.rate == 0) } 
  # Here, growth rate was not zero for any simulated data
  
  df.true <- data.frame(x = coords[, 1],
                        y = coords[, 2],
                        survival = simModel$sa,
                        w.as = simModel$w.as,
                        recruit = simModel$f,
                        w.af = simModel$w.af,
                        id = coords[,"id"],
                        X.cov = X.cov
  )
  
  
  #================================================================
  # Fit NNGP model to training data
  #================================================================
  
  coords.train <- (coords[1:N,]) 
  
  # Coordinates for plotting
  coords.train.sf <- st_as_sf(data.frame(x = coords.train[, 1], y = coords.train[, 2]),
                              coords = c('x', 'y'))
  
  df.train <- df.true[1:N, ]
  df.train$source = source[1:N]
  
  # Compute neighbors and sort coordinated in the training data 
  Neighbors = computeNeighbours(coords = df.train, k)
  neighbors.id <- Neighbors$neighbor_idx
  
  # sort locations and data
  coords.sorted <- Neighbors$coords_sorted[, 1:2]
  edist.train = Neighbors$edist_sorted # this is the distance matrix of the ordered coordinates.
  nid.dist =Neighbors$neighbors_dist# the last column is the distance to itself
  sorted.id <- Neighbors$coords_sorted[, c("id", "source")]
  
  
  # I now only want counts ordering so I just use the NNGP on 1:NC sites not 1:M
  coords.C <- coords.train[1:NC, ]
  
  k = 15
  # Selecting k neighbors based of Vecchia's choice. Ordered along x axis.  
  Neighbors.C = computeNeighbours(coords.C, k)
  #Neighbors.C$coords_sorted # The source 1 is counts and 2 is CR. id is the location in the orignal coordinates 
  neighbors.Cid <- Neighbors.C$neighbor_idx
  coords.sorted.C <- Neighbors.C$coords_sorted[, 1:2]
  edist.C = Neighbors.C$edist_sorted # this is the distance matrix of the ordered coordinates. 
  nid.Cdist =Neighbors.C$neighbors_dist# the last column is the distance to itself
  sorted.Cid <- Neighbors.C$coords_sorted[, 3:4]
  
  
  # Match data to ordered locations.
  
  # Covariates 
  X.cov.train <- X.cov[sorted.id$id, ,]
  dim(X.cov.train)
  
  # Data 
  
  # Counts 
  C.sorted <- which(sorted.id$source=="C") # indicator for which sites are counts sites in sorted data
  C.id <- sorted.id$id[which(sorted.id$source =="C")]  # indicator for which sites are counts sites in orignial data but also matches order in sorted data
  C.all <- matrix(NA, M, n.occasions)
  C.all[C.sorted, ] = simModel$C[C.id, ]
  
  C = matrix(NA, nrow = N, ncol = n.occasions)
  #c.level = 5/6 # high
  c.level = 3/6 # low
  
  for (i in C.sorted) {
    # Temporal Misalignment
    n.ones <- (n.occasions)*c.level
    id <- sample(1:(n.occasions), n.ones)
    C[ i, id ] = c(C.all[i, id])
    
  }
  
  C.na = C
  C = C.na[!is.na(C.na)]
  #C # NA's removed by site, columns
  n.counts = length(C)
  
  # Indicators of N for C to be scalar
  site = c()
  yr = c()
  
  for (d in 1:n.occasions) {
    
    a = which(!is.na(C.na[,d]))
    b = rep(d, length(a))
    site = c(site, a)
    yr = c(yr, b)
    
  }
  
  
  # Capture recapture ordering
  CR.sorted <- which(sorted.id$source== "CR") # indicator for which sites are CR sites in sorted data
  CR.id <- sorted.id$id[which(sorted.id$source == "CR")]  # indicator for which sites are CR sites in original data but also matches order in sorted data
  
  Rm <- array(0, dim = c(N, n.occasions-1, n.occasions ))
  Rm[ CR.sorted, ,] <- simModel$Rm[CR.id, ,]#Rmaps[CR.id-NC, ,]
  
  
  rel.R <- matrix(0, nrow = N, ncol = n.occasions-1)
  for (i in CR.sorted ) {
    rel.R[i, ] <- rowSums(Rm[i,,])
  }
  
  visits.train <- matrix(0, N, n.occasions-1)
  visits.train[CR.sorted, ] <- visits.all[CR.id, ]
  
  win.data = list( NC = NC, NCR =NCR, rel.R= rel.R, n.occasions=n.occasions, M=N,  mu = rep(0, N), 
                   edist = edist.train, dmax = max(edist.train), dmin = min(edist.train[edist.train!= 0] ), 
                   nid.dist = nid.dist, n.id = neighbors.id, k = k, C.sorted = C.sorted, CR.sorted = CR.sorted,
                   nid.Cdist = nid.Cdist, n.Cid = neighbors.Cid, edist.C = edist.C,
                   n.counts = n.counts,  site = site, yr = yr, ncov = ncov, visits = visits.train )
  
  
  sipm <- nimbleCode({
    
    
    for (i in 1:ncov) {
      beta.f[i] ~ dnorm(0, var = 2)
      beta.s[i] ~ dnorm(0, var = 2)
    }
    
    
    ################### 
    # GP's
    ###################
    
    # Recruitment
    alpha2.af ~ dinvgamma(3,1)
    rho.af ~ dunif(dmin/3, dmax/3)
    
    AD.af[1:(NC), 1:(k+1)] <- computeAD(edist.C[1:(NC),1:(NC)], nid.Cdist[1:(NC), 1:(k+1)], n.Cid[1:(NC), 1:k], rho.af, alpha2.af, k )
    w.af[1:(NC)] ~ dmnorm_NNGP(mu[1:(NC)], AD.af[1:( NC), 1:(k+1)], n.Cid[1:(NC), 1:k])
    constraint.f ~ dnorm(sum( w.af[1:NC]), sd = 1)
    
    
    # Survival
    alpha2.as ~ dinvgamma(3,1)
    rho.as ~ dunif(dmin/3, dmax/3)
    
    AD.as[1:(M), 1:(k+1)] <- computeAD(edist[1:(M),1:(M)], nid.dist[1:(M), 1:(k+1)], n.id[1:(M), 1:k], rho.as, alpha2.as, k )
    w.as[1:(M)] ~ dmnorm_NNGP(mu[1:(M)], AD.as[1:(M), 1:(k+1)], n.id[1:(M), 1:k])
    constraint.s ~ dnorm(sum( w.as[1:M]), sd = 1)
    
    
    for (t in 1:(n.occasions-1)) {
      wst[t] ~ dnorm(0, var = 2)
      wft[t] ~ dnorm(0, var = 2)
    }
    
    
    for (t in 1:(n.occasions-1)) {
      alpha.p[t] ~ dnorm(0, var = 2)
    }
    
    
    for (i in 1:M) {
      
      # Priors and linear models
      for (t in 1:(n.occasions-1)){
        
        logit(sa[i,t]) <- wst[t] + (Xcov[i,t, 1:ncov] %*% beta.s[1:ncov])[1,1] + w.as[i] 
        p[i,t] <- ilogit(alpha.p[t]) * visits[i,t] 
      }
      
    }
    
    
    for (i in 1:NC) {
      
      u[C.sorted[i]] ~ dunif(50,150)
      N[C.sorted[i],1] <- round(u[C.sorted[i]])
      
      # Process model over time: our model of population dynamics
      for (t in 1:(n.occasions-1)){
        
        log(f[C.sorted[i],t]) <- wft[t] + (Xcov[C.sorted[i],t, 1:ncov] %*% beta.f[1:ncov])[1,1] + w.af[i] 
        #log(f[C.sorted[i],t]) <- wft[t] + X_tmax[C.sorted[i],t]*beta.f.tmax + X_long[C.sorted[i],t]*beta.f.long + X_lat[C.sorted[i],t]*beta.f.lat + X_tmax_long[C.sorted[i],t]*beta.f.tmaxlong + X_tmax_lat[C.sorted[i],t]*beta.f.tmaxlat  + w.af[i]
        N[C.sorted[i],t+1] ~ dpois(sa[C.sorted[i],t]*N[C.sorted[i],t] + f[C.sorted[i],t]*N[C.sorted[i],t])
        
        # Growth rate
        growth.rate[C.sorted[i],t] <- (N[C.sorted[i], t+1]) / (N[C.sorted[i],t] + 0.001)
      }
    }  
    
    for (j in 1:n.counts) {
      
      log(lambda[j]) <- log(N[site[j], yr[j]]) 
      C[j] ~ dpois(lambda[j])
      
    }
    
    for (i in 1:NCR) {
      
      for (t in 1:(n.occasions-1)){
        Rm[CR.sorted[i],t,1:n.occasions] ~ dmulti(r[CR.sorted[i],t,1:n.occasions], rel.R[CR.sorted[i],t])
      }
      
      # Define the cell probabilities of the m-array
      # Main diagonal
      for (t in 1:(n.occasions-1)){
        q[CR.sorted[i], t] <- 1-p[CR.sorted[i], t]                # Probability of non-recapture
        r[CR.sorted[i], t,t] <- sa[CR.sorted[i],t]*p[CR.sorted[i],t]
      }
      # Above main diagonal
      for (t in 1:(n.occasions-2)){
        for (j in (t+1):(n.occasions-1)){
          r[CR.sorted[i],t,j] <- prod(sa[CR.sorted[i],t:j])*prod(q[CR.sorted[i],t:(j-1)])*p[CR.sorted[i],j]
        } #j
      } #t
      # Below main diagonal
      for (t in 2:(n.occasions-1)){
        for (j in 1:(t-1)){
          r[CR.sorted[i],t,j] <- 0
        } #j
      } #t
      # Last column: probability of non-recapture
      for (t in 1:(n.occasions-1)){
        r[CR.sorted[i],t,n.occasions] <- 1-sum(r[CR.sorted[i],t,1:(n.occasions-1)])
      } #t
      
    }
    
    
  })
  
  {
    set.seed(8)
    
    inits <- list(beta.f = rep(0, ncov) , beta.s = rep(0, ncov),   u = runif(win.data$M,50, 150), 
                  rho.af = mean(edist.C)/3, alpha2.af = rinvgamma(1,3,1),  rho.as = mean(edist)/3, alpha2.as = rinvgamma(1,3,1),
                  alpha.p = rep(0,n.occasions-1)  )
    
    
    
    inits$wst = nimble::logit(runif(n.occasions-1, 0.6,0.9) )
    inits$wft = log(runif(n.occasions-1,0.3,0.8))
    
    inits$w.af= MASS::mvrnorm(1, mu =rep(0,NC), Sigma = expcov(edist.C[1:(NC),1:(NC)], inits$rho.af, inits$alpha2.af))
    inits$w.as= MASS::mvrnorm(1, mu =rep(0,N), Sigma = expcov(edist[1:(N),1:(N)], inits$rho.as, inits$alpha2.as))
    
    
    # Initial values for R,S
    Ninits = Rinits = Sinits = array(NA, dim = c(win.data$M, win.data$n.occasions)) 
    Ninits[ ,1] = round(inits$u)
    
    for (t in 1:(win.data$n.occasions-1)) {
      Rinits[ ,t+1] =  rpois(win.data$M,  exp(inits$wft[t])*Ninits[,t])
      Sinits[ ,t+1] =  rbinom(win.data$M, Ninits[,t] , ilogit(inits$wst[t]) ) 
      Ninits[, t+1] = Rinits[,t+1] + Sinits[, t+1]
    }
    
    #Ninits
    inits$N = Ninits
    inits1 = inits
    
    set.seed(5)
    
    inits <- list(beta.f = rep(0, ncov) , beta.s = rep(0, ncov),   u = runif(win.data$M, 50, 150), 
                  rho.af = mean(edist.C)/3, alpha2.af = rinvgamma(1,3,1),  rho.as = mean(edist)/3, alpha2.as = rinvgamma(1,3,1),
                  alpha.p = rep(0,n.occasions-1)  )
    
    inits$wst = nimble::logit(runif(n.occasions-1, 0.6,0.9) )
    inits$wft = log(runif(n.occasions-1,0.3,0.8))
    
    inits$w.af= MASS::mvrnorm(1, mu =rep(0,NC), Sigma = expcov(edist.C[1:(NC),1:(NC)], inits$rho.af, inits$alpha2.af))
    inits$w.as= MASS::mvrnorm(1, mu =rep(0,N), Sigma = expcov(edist[1:(N),1:(N)], inits$rho.as, inits$alpha2.as))
    
    
    # Initial values for R,S
    Ninits = Rinits = Sinits = array(NA, dim = c(win.data$M, win.data$n.occasions)) 
    Ninits[ ,1] = round(inits$u)
    
    for (t in 1:(win.data$n.occasions-1)) {
      Rinits[ ,t+1] =  rpois(win.data$M,  exp(inits$wft[t])*Ninits[,t])
      Sinits[ ,t+1] =  rbinom(win.data$M, Ninits[,t] , ilogit(inits$wst[t]) ) 
      Ninits[, t+1] = Rinits[,t+1] + Sinits[, t+1]
    }
    
    #Ninits
    inits$N = Ninits
    
    inits2 = inits
    
    
  }
  
  
  #=========================================================================================================================================================================
  # Parallization
  #=========================================================================================================================================================================
  
  #First we create a cluster, specifying the number of cores we want the cluster to operate across.
  this_cluster <- makeCluster(2)
  
  per_chain_info <- list(
    list(seed = 8,
         inits = inits1  ),
    list(seed = 5,
         inits = inits2)
  )
  
  # Create a function with all the needed code
  run_MCMC_allcode <- function(info, Rmdata, Cdata, Xdata, win.data, code) {
    
    library(nimble)
    
    expcov <- nimbleFunction(
      run = function(dists = double(2), rho = double(0), sigma2 = double(0)) {
        returnType(double(2))
        n <- dim(dists)[1]
        y <- dim(dists)[2]
        result <- matrix(nrow = n, ncol = y, init = FALSE)
        for(i in 1:n)
          for(j in 1:y)
            result[i, j] <- sigma2*exp(-dists[i,j]/rho)
        return(result)
      })
    
    assign('expcov', expcov, .GlobalEnv)
    
    
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
    
    assign('computeC', computeC, .GlobalEnv)
    
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
    
    assign('computeAD', computeAD, .GlobalEnv)
    
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
    
    assign('computeQF', computeQF, .GlobalEnv)
    
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
    
    assign('dmnorm_NNGP', dmnorm_NNGP, .GlobalEnv)
    
    rmnorm_NNGP<- nimbleFunction(
      run = function(n = integer(0),  mu = double(1), AD = double(2), neighbors.id = double(2)) {
        stop('rmnorm_NNGP should never be run')
        returnType(double(1))
        return(rep(0))
      })
    
    #compileNimble(rmnorm_NNGP)
    assign('rmnorm_NNGP', rmnorm_NNGP, .GlobalEnv)
    
    FullModel = nimbleModel(code, 
                            data = list(Rm= Rmdata,  C = Cdata, Xcov = Xdata, constraint.s = 0, constraint.f = 0  ), 
                            constants = win.data,
                            inits = info$inits )
    
    FullModelconf = configureMCMC(FullModel)  
    FullModelconf$resetMonitors()
    FullModelconf$addMonitors( "alpha.p",  "wst", "wft", "rho.as",  "alpha2.as", "rho.af",  "alpha2.af", "beta.s", "beta.f"   )
    FullModelconf$addMonitors2("N",   "sa", "f", "growth.rate", "w.as", "w.af" ) 
    
    
    FullModelconf$removeSamplers("w.af", "w.as")
    
    for (i in 1:(win.data$M) ) {
      FullModelconf$addSampler(target = paste0("w.as[", i, "]"), type = "RW")
    }
    
    for (i in 1:win.data$NC) {
      FullModelconf$addSampler(target = paste0("w.af[", i, "]"), type = "RW")
    }
    
    FullModelconf$removeSamplers("beta.s", "beta.f" )
    beta.nodes <- c( paste0("beta.s[",1:win.data$ncov, "]" ), paste0("beta.f[",1:win.data$ncov, "]" ) )
    FullModelconf$addSampler(target = beta.nodes, type = "RW_block")
    
    FullModelmcmc = buildMCMC(FullModelconf)
    cFullModel <- compileNimble(FullModel)
    FullModelmod = compileNimble(FullModelmcmc, project = FullModel, resetFunctions = TRUE)
    
    results <- runMCMC(FullModelmod,niter = 125000, nburnin = 80000,thin = 15, thin2 = 25, setSeed = info$seed)
    
    return(results)
  }
  
  
  
  NNGP_time  <- system.time (chain_output <- parLapply(cl = this_cluster, X = per_chain_info, 
                                                       fun = run_MCMC_allcode, 
                                                       win.data = win.data,
                                                       code = sipm, 
                                                       Rmdata = Rm,
                                                       Cdata = C, 
                                                       Xdata = X.cov.train))
  
  # It's good practice to close the cluster when you're done with it.
  stopCluster(this_cluster)
  
  Train_time[n,] = NNGP_time[1:3] # training time
  save(Train_time, file = "train_time_STmis_high_sim_int.Rdata")
  
  
  c1 <- mcmc(chain_output[[1]]$samples)
  c2 <- mcmc(chain_output[[2]]$samples)
  samples_og <- mcmc.list(list(c1 = c1, c2 = c2))
  res = MCMCvis::MCMCsummary(samples_og)
  key_parameters[,,n] = as.matrix(res)
  
  save(key_parameters, file = "key_parameters_STmis_high_sim_int.Rdata")
  
  
  c1 <- mcmc(chain_output[[1]]$samples2)
  c2 <- mcmc(chain_output[[2]]$samples2)
  samples2_og <- mcmc.list(list(c1 = c1, c2 = c2))
  samples2 = mcmc(do.call(rbind, samples2_og) )
  
  
  # Gaussian Process parameters
  
  ws.est = MCMCvis::MCMCsummary(samples2[, grep("w.as", colnames(samples2)) ], Rhat = FALSE)
  rb = (ws.est$mean - df.true$w.as[sorted.id$id])/ df.true$w.as[sorted.id$id]
  w.as.est[n, 1] =  mean(rb)
  w.as.est[n, 2:3] =  quantile(rb, probs = c(0.025, 0.975))
  w.as.est[n, 4] = sqrt(mean((ws.est$mean - df.true$w.as[sorted.id$id])^2)) 
  w.as.est[n, 5 ] =  coverage(ws.est,df.true$w.as[sorted.id$id])
  
  
  wf.est = MCMCvis::MCMCsummary(samples2[, grep("w.af", colnames(samples2))], Rhat = FALSE)
  rb = (wf.est$mean - df.true$w.af[C.id])/ df.true$w.af[C.id]
  w.af.est[n, 1] =  mean(rb)
  w.af.est[n, 2:3] =  quantile(rb, probs = c(0.025, 0.975))
  w.af.est[n, 4] = sqrt(mean((wf.est$mean -  df.true$w.af[C.id])^2)) 
  w.af.est[n, 5 ] =  coverage(wf.est,  df.true$w.af[C.id])
  
  save(w.as.est, file = "w_as_est_STmis_high_sim_int.Rdata")
  save(w.af.est, file = "w_af_est_STmis_high_sim_int.Rdata")
  
  
  # Survival 
  sa = as.matrix((df.true[sorted.id$id, 3:7]))
  sa_est = MCMCvis::MCMCsummary(samples2[, grep("sa", colnames(samples2) ) ], Rhat = FALSE)
  rb = (sa_est$mean - sa)/sa
  phi_est[n, 1] = mean(rb)
  phi_est[n, 2:3] =  quantile(rb, probs = c(0.025, 0.975))
  phi_est[n, 4] = sqrt(mean((sa_est$mean - sa)^2))
  phi_est[n, 5] = coverage(sa_est, sa)
  phi_est[n, 6:10] = colMeans(rb) # last year bias
  
  
  # Recruitment
  f = as.matrix((df.true[C.id, 9:13]))
  fmcmc <-samples2[, grep("f", colnames(samples2) ) ]
  
  f.id <- c()
  for (t in 1:(n.occasions-1)) {
    for (i in C.sorted) {
      f.id <- c(f.id, paste0("f[", i, ", ", t, "]"  )    )
    }
    
  }
  
  f_est  = MCMCvis::MCMCsummary(fmcmc[, f.id], Rhat = FALSE) 
  rb = (f_est$mean - f)/f
  rec_est[n, 1] = mean(rb)
  rec_est[n, 2:3] =  quantile(rb, probs = c(0.025, 0.975))
  rec_est[n, 4] = sqrt(mean((f_est$mean - f)^2))
  rec_est[n, 5] = coverage(f_est, f)
  rec_est[n, 6:10] = colMeans(rb) # last year bias
  
  # Growth rate
  gr = as.matrix(simModel$growth.rate[C.id, ])   #as.matrix((df.true[C.id, 9:13]))
  grmcmc <-samples2[, grep("gr", colnames(samples2) ) ]
  
  gr.id <- c()
  for (t in 1:(n.occasions-1)) {
    for (i in C.sorted) {
      gr.id <- c(gr.id, paste0("growth.rate[", i, ", ", t, "]"  )    )
    }
    
  }
  
  gr_mc  = MCMCvis::MCMCsummary(grmcmc[, gr.id], Rhat = FALSE) 
  rb = (gr_mc$mean - gr)/(gr)
  gr_est[n, 1] = mean(rb)
  gr_est[n, 2:3] =  quantile(rb, probs = c(0.025, 0.975))
  gr_est[n, 4] = sqrt(mean((gr_mc$mean - gr)^2))
  gr_est[n, 5] = coverage(gr_mc, gr)
  gr_est[n, 6:10] = colMeans(rb) # last year bias
  
  
  # Abundance 
  N.true = as.matrix(simModel$N[C.id, ])   #as.matrix((df.true[C.id, 9:13]))
  Nmcmc <-samples2[, grep("N", colnames(samples2) ) ]
  
  N.id <- c()
  for (t in 1:(n.occasions)) {
    for (i in C.sorted) {
      N.id <- c(N.id, paste0("N[", i, ", ", t, "]"  )    )
    }
    
  }
  
  
  N_mc  = MCMCvis::MCMCsummary(Nmcmc[, N.id], Rhat = FALSE) 
  rb = (N_mc$mean - N.true)/(N.true)
  N_est[n, 1] = mean(rb)
  N_est[n, 2:3] =  quantile(rb, probs = c(0.025, 0.975))
  N_est[n, 4] = sqrt(mean((N_mc$mean - N.true)^2))
  N_est[n, 5] = coverage(N_mc, N.true)
  N_est[n, 6:11] = colMeans(rb) # last year bias
  
  save(phi_est, file = "phi_est_STmis_high_sim_int.Rdata")
  save(rec_est, file = "rec_est_STmis_high_sim_int.Rdata")
  save(gr_est, file = "gr_est_STmis_high_sim_int.Rdata")
  save(N_est, file = "N_est_STmis_high_sim_int.Rdata")
  
  #=========================================================================================================================================================================
  # Predictions
  #=========================================================================================================================================================================
  
  # Randomly sample location with the BCR 30 
  N.p = P
  
  coords.predict = coords[(N+1):M,]
  
  # Coordinates for plotting
  coords.sf.predict <- st_as_sf(data.frame(x = coords.predict[, 1], y = coords.predict[, 2], Name = as.factor(1:N.p)),
                                coords = c('x', 'y'))
  
  
  # Let's say I want to make prediction at these locations in the first year
  X.cov.pred <- X.cov[(N+1):M, ,]
  primary.samples <- mcmc(do.call(rbind, samples_og[2001:3000, ]))
  secondary.samples <- mcmc(do.call(rbind, samples2_og[801: 1800, ]))
  
  
  w.as.MCMC <- secondary.samples[, paste0("w.as[",1:N,"]" )]
  w.af.MCMC <- secondary.samples[, paste0("w.af[",1:NC,"]" )]
  
  
  MCMC_nsamples = dim(primary.samples)[1]
  phi.pred = array(0, dim = c(MCMC_nsamples, N.p, n.occasions-1)) #matrix(0, nrow = MCMC_nsamples, N.p)
  f.pred = array(0, dim = c(MCMC_nsamples, N.p, n.occasions-1))   #matrix(0, nrow = MCMC_nsamples, N.p)
  
  print("Computing predictions...")
  
  for (t in 1:(n.occasions-1)) {
    
    for (i in 1:N.p) {
      
      #print(paste0("Year ", t,  " Prediction location ", i))
      
      #=====================================
      # Preamble for survival
      #=====================================
      
      # Add x0 to training data to find k nearest neighors
      coords.x0 = rbind(Neighbors$coords_sorted[, c("x", "y", "id")], coords.predict[i,c("x", "y", "id" )]) # add point to training data
      
      # Sort coordinates by x axis
      sorted_idx <- order(coords.x0[,1])
      x0.location = which(sorted_idx == N+1) # record location so easy to pick out necessary location
      coords_sorted =coords.x0[sorted_idx, ]
      #trdata = trdata[sorted_idx]
      
      coords.sf <- st_as_sf(data.frame(x = coords_sorted[, 1], y = coords_sorted[, 2]),
                            coords = c('x', 'y'))
      
      # Compute distance between all points
      edist = matrix(as.numeric(st_distance(coords.sf)), N+1, N+1)
      #x0.location
      
      if(x0.location <= (k+1)){
        
        neighbors_idx.x0 <- as.numeric(1:(x0.location -  1)) # k nearest neighbors of x0 - N0
        neighbors_dist.x0 <-  edist[x0.location, neighbors_idx.x0] # dist between x0 and its neighbors
        neighbors_dist.N0 <- edist[neighbors_idx.x0, neighbors_idx.x0 ] # distance between neighbors of x0
        #s_x0 = trdata[neighbors_idx.x0 ] # training data values of neighbors of x0
        
      } else {
        
        curr_sorted_idx = order(edist[x0.location, ])
        curr_sorted_idx = curr_sorted_idx[curr_sorted_idx < x0.location] # take values to the left
        neighbors_idx.x0 = curr_sorted_idx[1:k]
        neighbors_dist.x0 <-  edist[x0.location, neighbors_idx.x0] # dist between x0 and its neighbors
        neighbors_dist.N0 <- edist[neighbors_idx.x0, neighbors_idx.x0 ] # distance between neighbors of x0
        #s_x0 = trdata[neighbors_idx.x0 ] # training data values of neighbors of x0
        
      }
      
      #=====================================
      # Preamble for recruitment
      #=====================================
      
      coords.x0.f = rbind(Neighbors.C$coords_sorted[, c("x", "y", "id")], coords.predict[i, c("x", "y", "id" )]) # add point to training data
      
      # Sort coordinates by x axis
      sorted_idx.f <- order(coords.x0.f[,1])
      x0.location.f = which(sorted_idx.f == NC+1) # record location so easy to pick out necessary location
      coords_sorted.f =coords.x0.f[sorted_idx.f, ]
      #trdata = trdata[sorted_idx]
      
      coords.sf.f <- st_as_sf(data.frame(x = coords_sorted.f[, 1], y = coords_sorted.f[, 2]),
                              coords = c('x', 'y'))
      
      # Compute distance between all points
      edist.f = matrix(as.numeric(st_distance(coords.sf.f)), NC+1, NC+1)
      #x0.location
      
      if(x0.location.f <= (k+1)){
        
        neighbors_idx.x0.f <- as.numeric(1:(x0.location.f -  1)) # k nearest neighbors of x0 - N0
        neighbors_dist.x0.f <-  edist.f[x0.location.f, neighbors_idx.x0.f] # dist between x0 and its neighbors
        neighbors_dist.N0.f <- edist.f[neighbors_idx.x0.f, neighbors_idx.x0.f ] # distance between neighbors of x0
        #s_x0 = trdata[neighbors_idx.x0 ] # training data values of neighbors of x0
        
      } else {
        
        curr_sorted_idx.f = order(edist.f[x0.location.f, ])
        curr_sorted_idx.f = curr_sorted_idx.f[curr_sorted_idx.f < x0.location.f] # take values to the left
        neighbors_idx.x0.f = curr_sorted_idx.f[1:k]
        neighbors_dist.x0.f <-  edist.f[x0.location.f, neighbors_idx.x0.f] # dist between x0 and its neighbors
        neighbors_dist.N0.f <- edist.f[neighbors_idx.x0.f, neighbors_idx.x0.f ] # distance between neighbors of x0
        #s_x0 = trdata[neighbors_idx.x0 ] # training data values of neighbors of x0
        
      }
      
      for (mc in 1:MCMC_nsamples) {
        
        rho.as = primary.samples[mc,"rho.as"]
        alpha2.as =primary.samples[mc, "alpha2.as"]
        
        rho.af = primary.samples[mc,"rho.af"]
        alpha2.af =primary.samples[mc, "alpha2.af"]
        
        # For the 1st year
        wst <-  primary.samples[mc, paste0("wst[", t, "]")   ]
        wft <-  primary.samples[mc,  paste0("wft[", t, "]") ]
        
        beta.f <- primary.samples[mc, grep("beta.f", colnames(primary.samples))]
        beta.s <- primary.samples[mc, grep("beta.s", colnames(primary.samples))]
        
        w.as.x0 <- w.as.MCMC[mc, neighbors_idx.x0] # This is spatial effects for x0 neighbours
        w.af.x0 <- w.af.MCMC[mc, neighbors_idx.x0.f] # This is spatial effects for x0 neighbours
        
        
        
        gp.pred.as <- NNGP.pred(neighbors_dist.x0, neighbors_dist.N0, x0.location, rho.as, alpha2.as, w.as.x0)
        phi.pred[mc, i, t] <- ilogit(wst + X.cov.pred[i,t ,] %*% beta.s + gp.pred.as)
        
        
        gp.pred.af <- NNGP.pred(neighbors_dist.x0.f, neighbors_dist.N0.f, x0.location.f, rho.af, alpha2.af, w.af.x0)
        f.pred[mc,i , t] <- exp( wft + X.cov.pred[i,t ,] %*% beta.f  + gp.pred.af)
        
      }
      
      
    }
    
    
  }
  
  result.predict = list(phi.pred = phi.pred, f.pred = f.pred, coords.sf.predict = coords.sf.predict)
  #====================================================================================================================================================
  # Use Mean absolute percentage error (MAPE), RB, Root mean square error (RMSE)
  #====================================================================================================================================================
  
  df.true.s.pred <- df.true[(N+1):M, 3:7]
  df.true.f.pred <- df.true[(N+1):M, 9:13]
  
  for (t in 1:(n.occasions-1)) {
    
    phi.predict <-  apply(result.predict$phi.pred[,,t], 2, mean)
    phi.true.p <- df.true.s.pred[, t]
    
    rb = ( phi.predict - phi.true.p )/phi.true.p
    phi.predict.est[n, 1,t] = mean(rb)
    phi.predict.est[n, 2:3,t] = quantile(rb, probs = c(0.025, 0.975))
    
    MAPE.s = mean((abs(phi.predict - phi.true.p)/phi.true.p)*100)
    phi.predict.est[n, 4,t] = MAPE.s
    RMSE.s = sqrt(mean((phi.predict - phi.true.p)^2))
    phi.predict.est[n, 5,t] = RMSE.s
    phi.predict.est[n, 6, t] = coverage(MCMCsummary(result.predict$phi.pred[,,t], Rhat = FALSE, n.eff = FALSE), phi.true.p)
    
    
    f.predict <-  apply(result.predict$f.pred[,,t], 2, mean)
    f.true.p <- df.true.f.pred[,t]
    
    rb = ( f.predict - f.true.p )/f.true.p
    rect.predict.est[n, 1,t] = mean(rb)
    rect.predict.est[n, 2:3,t] = quantile(rb, probs = c(0.025, 0.975))
    MAPE.f = mean((abs(f.predict - f.true.p)/f.true.p)*100)
    rect.predict.est[n, 4,t] = MAPE.f
    RMSE.f = sqrt(mean((f.predict - f.true.p)^2))
    rect.predict.est[n, 5, t] = RMSE.f
    rect.predict.est[n, 6, t] =  coverage(MCMCsummary(result.predict$f.pred[,,t], Rhat = FALSE, n.eff = FALSE), f.true.p)
    
    
  }
  
  save(phi.predict.est, file = "phi_predict_STmis_high_sim_int.Rdata")
  save(rect.predict.est, file = "rect_predict_STmis_high_sim_int.Rdata")
  
  
  
  
  
}
