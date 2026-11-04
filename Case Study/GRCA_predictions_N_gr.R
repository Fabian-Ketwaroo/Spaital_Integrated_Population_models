
#=====================================================================================================================================================================================
# Code to make predictions of relative abundance and population growth rate using Bayesian Inverse distance 
#=====================================================================================================================================================================================

rm(list = ls())
# Load packages
{
  library(nimble)
  library(sf)
  library(MCMCvis)
  library(coda)
  library(MASS)
  library(dplyr)
  library(climateR)
  library(ggplot2)
  library(lubridate)
  library(terra)
  library(stars)
  library(gstat)
  
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
    
    
    edist = matrix(as.numeric(st_distance(coords_sorted.sf)/1000), N, N)
    
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
  
  
}

# Working directory should be the project root. Open the .Rproj file, or use
# setwd() to point to the folder containing "Case Study" and "Simulation study".

# Load data 
# The 'GRCA_data_2004_2014.Rdata' file is included in the 'Case study' directory.
load("Case Study/GRCA_data_2004_2014.Rdata")

# Number of count locations 
NC = dim(GRCA_data$C)[1] 

# Number of Capture-recapture location
NCR =  dim(GRCA_data$Rm)[1] 

# Total number of locations
M <- NC + NCR

# Number of sampling occasions 
n.occasions =  dim(GRCA_data$C)[2] 


BCRid = c(13,28,29, 30, 27,14)

# Plot of locations

# Note: The 'BCR_Terrestrial_master.shp' file is included in the bcr_terrestrial_shape folder 
# of the 'Case Study' directory.
# Please ensure all associated sidecar files (.shx, .dbf, .prj) remain in the 
# same folder to allow GDAL to read the spatial geometries correctly.

shapefile_path  <- 'Case Study/bcr_terrestrial_shape/BCR_Terrestrial_master.shp'


# Read the shapefile into R
bcr_data <- st_read(shapefile_path)

# Filter for BCR 30 (assumes BCR code column is 'BCR')
bcr.all <- bcr_data %>% filter(BCR %in% BCRid)

# Get spatial co-ordinates of sampling locations: 
coords = GRCA_data$GRCA_splocations[, -3]  #GRCA_10yrs$GRCA_splocations
source = as.factor(c(rep("C", NC), rep("CR", NCR) ))
coords = cbind(coords, source)
coords.raw <- cbind(coords, id = 1:M)
proj = "+proj=longlat +datum=WGS84"
spatial_data <- st_transform(bcr.all, crs = proj)  # Transform to WGS84

# Coordinates for plotting
coords.sf <- st_as_sf(data.frame(x = coords[, 1], y = coords[, 2], Source = source),
                      coords = c('x', 'y'), crs = proj)

ggplot() +
  geom_sf(data = spatial_data, fill = "white") +
  geom_sf(data = coords.sf, aes(col = Source), size = 1.5) +
  theme_bw(base_size = 14) +
  labs(x = 'longitude', y = 'latitude') + scale_color_discrete(labels = c("BBS", "MAPS"))


# Here, the locations need to be ordered. I need to keep track of where the different data points are.
coords.raw <- cbind(coords, id = 1:M)
k = 15
# Selecting k neighbors based of Vecchia's choice. Ordered along x axis.  
Neighbors = computeNeighbours(coords.raw, k, proj)
str(Neighbors)


# plot to show neighbors
# red is ind location and blue are the neighbors
DisplayNeighbors(ind = 501, coords.sorted = Neighbors$coords_sorted, neighbor_idx = Neighbors$neighbor_idx, proj)
# This function would also be better if I set the spatial outline as an argument.


neighbors.id <- Neighbors$neighbor_idx

coords.sorted <- Neighbors$coords_sorted[, 1:2]
edist = Neighbors$edist_sorted # this is the distance matrix of the ordered coordinates. 
nid.dist =Neighbors$neighbors_dist# the last column is the distance to itself
sorted.id <- Neighbors$coords_sorted[, 3:4]
C.sorted <- which(sorted.id== "C") 

k = 15
# Selecting k neighbors based of Vecchia's choice. Ordered along x axis.  
Neighbors.C = computeNeighbours(coords.raw[1:NC,], k, proj)
str(Neighbors.C)

neighbors.Cid <- Neighbors.C$neighbor_idx

coords.sorted.C <- Neighbors.C$coords_sorted[, 1:2]
edist.C = Neighbors.C$edist_sorted # this is the distance matrix of the ordered coordinates. 

# Define grid for prediction

# Define the resolution of the grid (adjust 'cellsize' for more or fewer points)
cellsize <- 0.15  # Degrees, adjust based on desired grid density

# Create a grid within the extent of spatial_data
grid <- st_make_grid(spatial_data, cellsize = cellsize, what = "centers")

# Convert grid to sf object
grid_sf <- st_as_sf(data.frame(geometry = grid), crs = st_crs(spatial_data))

# Keep only grid points inside the polygons
grid_sf <- grid_sf[spatial_data, ]  

# Extract coordinates
grid_sf$x <- st_coordinates(grid_sf)[,1]  # Longitude
grid_sf$y <- st_coordinates(grid_sf)[,2]  # Latitude


ggplot() +
  geom_sf(data = spatial_data, fill = "lightgray", color = "black") +
  geom_sf(data = grid_sf, color = "blue", size = 0.5) +
  theme_minimal() +
  labs(title = "Generated Grid Inside spatial_data")

class(grid_sf) # sf good.


coords.predict = st_coordinates(grid_sf)
N.p <- dim(coords.predict)[1]
coords.sf.predict <- st_as_sf(data.frame(x = coords.predict[, 1], y = coords.predict[, 2], Name = as.factor(1:N.p)),
                              coords = c('x', 'y'), crs = proj)


load("Case Study/GRCA_samples2_NNGP_SV.Rdata")

secondary.samples <- mcmc(do.call(rbind, samples2))
dim(secondary.samples)

Nmcmc <-secondary.samples[, grep("N", colnames(secondary.samples) ) ]
grmcmc <-secondary.samples[, grep("growth.rate", colnames(secondary.samples) ) ]


MCMC_nsamples = dim(secondary.samples)[1]
N.pred = matrix(0, nrow = MCMC_nsamples, N.p)
gr.pred = matrix(0, nrow = MCMC_nsamples, N.p)


for (mc in 1:MCMC_nsamples) {
  
  N.mean <- matrix(0, length(C.sorted), n.occasions)
  gr.mean <- matrix(0, length(C.sorted), n.occasions-1)
  
  # Posterior mean
  for (t in 1:n.occasions) {
    
    N.id <- c()
    for (i in C.sorted) {
      N.id <- c(N.id, paste0("N[", i, ", ", t, "]"  )    )
    }
    
    N.mean[, t] <- Nmcmc[mc, N.id]
    
  }
  
  for (t in 1:(n.occasions-1)) {
    
    gr.id <- c()
    
    for (i in C.sorted) {
      gr.id <- c(gr.id, paste0("growth.rate[", i, ", ", t, "]"  )    )
    }
    
    gr.mean[, t] <- grmcmc[mc, gr.id]
    
  }
  
  head(N.mean)
  head(gr.mean)
  
  N.mean.sp <- rowMeans(N.mean)
  gr.mean.sp <- rowMeans(gr.mean)
  
  N.df <- data.frame(x = coords[C.sorted, 1],
                     y = coords[C.sorted, 2],
                     survival = N.mean.sp)
  gr.df <- data.frame(x = coords[C.sorted, 1],
                      y = coords[C.sorted, 2],
                      survival = gr.mean.sp)
  
  
  N.sf <- st_as_sf(N.df, coords = c('x', 'y'), crs = proj)
  gr.sf <- st_as_sf(gr.df, coords = c('x', 'y'), crs = proj)
  
 
  idw_model <- gstat(formula = survival ~ 1, locations = N.sf, nmax = 15, set = list(idp = 1))
  
  # Predict onto the grid
  N.pred[mc, ] <- predict(idw_model, newdata = grid_sf)$var1.pred
  
  idw_model <- gstat(formula = survival ~ 1, locations = gr.sf, nmax = 15, set = list(idp = 1))
  
  # Predict onto the grid
  gr.pred[mc, ] <- predict(idw_model, newdata = grid_sf)$var1.pred
  
  
}

# strings of [inverse distance weighted interpolation] harmless.

head(N.pred)
head(gr.pred)
result.predict = list(N.pred = N.pred, gr.pred = gr.pred, coords.sf.predict = coords.sf.predict)

#=============================================================================================================================================================================
# R code to produce spatially smooth prediction maps 
#=============================================================================================================================================================================

#=======================================================
# Relative Abundance 
#=======================================================

#Predicted locations
coords.predict = st_coordinates(result.predict$coords.sf.predict)


phi.mean <-  apply(result.predict$N.pred, 2, mean)
phi.sd <- apply(result.predict$N.pred, 2, sd)

phi.pred.df <- data.frame(x = coords.predict[, 1],
                          y = coords.predict[, 2],
                          Mean = phi.mean,
                          Sd = phi.sd)


phi.pred.sf <- st_as_sf(phi.pred.df, coords = c('x', 'y'), crs = st_crs(spatial_data))

phi.pred.plot.mean <-ggplot(data = spatial_data) + geom_sf(  ) +
  geom_sf(data = phi.pred.sf, aes(col = Mean), size= 1) +
  theme_bw(base_size = 14) +
  scale_color_gradient2(midpoint = mean(phi.pred.df$Mean, na.rm = TRUE), low = '#B2182B', mid = 'white', high = '#2166AC', na.value = NA) +
  labs(x = 'x', y = 'y', col = 'Mean') 

phi.pred.plot.mean


phi.pred.sf$x <- st_coordinates(phi.pred.sf)[,1] 
phi.pred.sf$y <- st_coordinates(phi.pred.sf)[,2] 


raster_data <- st_as_stars(phi.pred.sf, dims = c("x", "y"))


ggplot() +
  geom_stars(data = raster_data, aes(fill = Mean), interpolate = TRUE) +
  coord_sf() +
  theme_minimal() +
  geom_sf(data = st_cast(spatial_data, "MULTILINESTRING")) +
  scale_fill_gradient2(low = "red", mid = "grey90", high = "blue", midpoint = mean(phi.pred.df$Mean, na.rm = TRUE), na.value = NA) +
  labs(x = 'Longitude', y = 'Latitude') + theme(legend.text = element_text(size = 10),  legend.title = element_text(size = 14))


# Sd

ggplot() +
  geom_stars(data = raster_data, aes(fill = Sd)) +
  geom_sf(data = st_cast(spatial_data, "MULTILINESTRING")) +
  scale_fill_gradient2(low = "red", mid = "grey90", high = "blue", midpoint = mean(phi.pred.df$Sd, na.rm = TRUE), na.value = NA, limits = range(phi.pred.df$Sd, na.rm = TRUE)) + 
  coord_sf() +
  theme_minimal() +  labs(x = 'Longitude', y = 'Latitude')+ theme(legend.text = element_text(size = 10),  legend.title = element_text(size = 14))


#=======================================================
# Growth rate 
#=======================================================



f.mean <-  apply(result.predict$gr.pred, 2, mean)
f.sd <- apply(result.predict$gr.pred, 2, sd)

f.pred.df <- data.frame(x = coords.predict[, 1],
                        y = coords.predict[, 2],
                        Mean = f.mean,
                        Sd = f.sd)


f.pred.sf <- st_as_sf(f.pred.df, coords = c('x', 'y'), crs = st_crs(spatial_data))

f.pred.plot.mean <-ggplot(data = spatial_data) + geom_sf(  ) +
  geom_sf(data = f.pred.sf, aes(col = Mean), size= 1) +
  theme_bw(base_size = 14) +
  scale_color_gradient2(midpoint = mean(f.pred.df$Mean, na.rm = TRUE), low = '#B2182B', mid = 'white', high = '#2166AC', na.value = NA) +
  labs(x = 'x', y = 'y', col = 'Mean') 

f.pred.plot.mean


f.pred.sf$x <- st_coordinates(f.pred.sf)[,1] 
f.pred.sf$y <- st_coordinates(f.pred.sf)[,2] 


raster_data_f <- st_as_stars(f.pred.sf, dims = c("x", "y"))


ggplot() +
  geom_stars(data = raster_data_f, aes(fill = Mean)) +
  coord_sf() +
  theme_minimal() +
  geom_sf(data = st_cast(spatial_data, "MULTILINESTRING")) +
  scale_fill_gradient2(low = "red", mid = "grey90", high = "blue", midpoint = mean(f.pred.df$Mean, na.rm = TRUE), na.value = NA, limits = range(f.pred.df$Mean, na.rm = TRUE))  + 
  labs(x = 'Longitude', y = 'Latitude')  + theme(legend.text = element_text(size = 10),  legend.title = element_text(size = 14))


# Sd

ggplot() +
  geom_stars(data = raster_data_f, aes(fill = Sd)) +
  geom_sf(data = st_cast(spatial_data, "MULTILINESTRING")) +
  scale_fill_gradient2(low = "red", mid = "grey90", high = "blue", midpoint = mean(f.pred.df$Sd, na.rm = TRUE), na.value = NA, limits = range(f.pred.df$Sd, na.rm = TRUE)) + labs(x = 'Longitude', y = 'Latitude') +
  coord_sf() +
  theme_minimal() +theme(legend.text = element_text(size = 10),  legend.title = element_text(size = 14))
