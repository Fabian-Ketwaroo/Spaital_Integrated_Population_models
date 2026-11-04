
# Code to make predictions of Survival and Recruitment via Kriging 

{
  #rm(list = ls())
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
hist(edist)
nid.dist =Neighbors$neighbors_dist# the last column is the distance to itself
sorted.id <- Neighbors$coords_sorted[, 3:4]


k = 15
# Selecting k neighbors based of Vecchia's choice. Ordered along x axis.  
Neighbors.C = computeNeighbours(coords.raw[1:NC,], k, proj)
str(Neighbors.C)

Neighbors.C$coords_sorted # The source 1 is counts and 2 is CR. id is the location in the orignal coordinates 


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


# tmax [degC] (air_temperature)
out.sp <- getTerraClim(coords.sf.predict, varname = 'tmax',  startDate = "2004-05-01",  endDate   = "2014-08-30" )

plot(out.sp$tmax[[1]], xlab= "Longitude", ylab = "Latitude")
plot(vect(coords.sf.predict), add = TRUE, pch = 16, cex = 0.2)

tmp = extract_sites(r = out.sp, pts = coords.sf.predict, ID = "Name")
tmp[1:5, 1:5]
dim(tmp) # rows: months, column: sites 

#To make the data ‘tidy’ simply pivot on the date column:
tmax = tidyr::pivot_longer(tmp, -date)
head(tmax)
dim(tmax)


tmax = cbind(tmax, month = month(tmax$date), year = year(tmax$date))
head(tmax)

Breed_seaon <- filter(tmax, month %in% c(5, 6, 7, 8)) # May - Aug
head(Breed_seaon)
dim(Breed_seaon)


years <- 2004:2014
nyears <- length(2004:2014)
Tmax_cov = matrix(0, nrow = N.p, ncol = nyears)
colnames(Tmax_cov) = paste0("Y", 13:23)

for (i in 1:N.p) {
  
  for (t in 1:nyears) {
    
    a <- filter(Breed_seaon, name ==i & year == years[t]  )
    Tmax_cov[i, t] <- median(a$value) # take the median TMAX over May-Aug
    
  }
  
}


long <- (coords.predict)[,1] #coords$longitude
lat <- (coords.predict)[,2] #coords$latitude

# Linear model to get residual
tmax.model_yr <- lm(Tmax_cov[, 1:10] ~ lat + long)
tmax.model_yr

# Residual TMAX
Tmax.resid <- tmax.model_yr$residuals



TMAX.predict.avg = rowMeans(Tmax.resid[, 1:10]) 
coords.predict = cbind(coords.predict, id = 1:N.p,  TMAX = TMAX.predict.avg)

coords.sf.tmax <- st_as_sf(data.frame(x = coords.predict[, 1], y = coords.predict[, 2], TMAX = coords.predict[,4]),
                           coords = c('x', 'y'), crs = proj)



tmax.pred.df <- data.frame(x = coords.predict[, 1],
                           y = coords.predict[, 2],
                           TMAX = coords.predict[, 4])


tmax.pred.sf <- st_as_sf(tmax.pred.df, coords = c('x', 'y'), crs = st_crs(spatial_data))

tmax.pred.plot.mean <-ggplot(data = spatial_data) + geom_sf(  ) +
  geom_sf(data = tmax.pred.sf, aes(col = TMAX), size= 1) +
  theme_bw(base_size = 14) +
  scale_color_gradient2(midpoint = mean(tmax.pred.df$TMAX, na.rm = TRUE), low = 'blue', mid = 'grey90', high = 'red', na.value = NA) +
  labs(x = 'Longitude', y = 'Latitude', col = 'TMAX') 

tmax.pred.plot.mean


tmax.pred.sf$x <- st_coordinates(tmax.pred.sf)[,1] 
tmax.pred.sf$y <- st_coordinates(tmax.pred.sf)[,2] 
TMAX_stars <- st_as_stars(tmax.pred.sf, dims = c("x", "y"))

ggplot() +
  geom_stars(data = TMAX_stars, aes(fill = TMAX)) +
  coord_sf() +
  theme_minimal() +
  geom_sf(data = st_cast(spatial_data, "MULTILINESTRING")) +
  scale_fill_gradient2(name = "Residual\nTMAX", low = "blue", mid = "grey90", high = "red", midpoint = mean(tmax.pred.sf$TMAX, na.rm = TRUE), na.value = NA, limits = range(tmax.pred.sf$TMAX, na.rm = TRUE))  + 
  labs(x = 'Longitude', y = 'Latitude')  + theme(legend.text = element_text(size = 10),  legend.title = element_text(size = 14))


# So I have: locations and TMAX at areas of predictions

X_tmax = TMAX.predict.avg
X_long = coords.predict[, "X"]
X_lat = coords.predict[, "Y"]
X_tmax_long <- X_tmax * X_long
X_tmax_lat <- X_tmax * X_lat

X_tmax <- scale(X_tmax) # scale
X_long <- scale(X_long)
X_lat <- scale(X_lat)
X_tmax_lat <- scale(X_tmax_lat)
X_tmax_long <- scale(X_tmax_long)

TMAX <- list(X_tmax = X_tmax, X_long = X_long, X_lat = X_lat, X_tmax_lat = X_tmax_lat, X_tmax_long = X_tmax_long)
#save(TMAX, file = "TMAX_residual_GRCA_11yrs.Rdata")

#==========================================================================================================

# Load MCMC samples obtained from running GRCA_SV_NNGP.R

load("Case Study/GRCA_main_samples_NNGP_SV.Rdata")
load("Case Study/GRCA_samples2_NNGP_SV.Rdata")



primary.samples <-mcmc(do.call(rbind, samples[3501:5000, ]))
secondary.samples <- mcmc(do.call(rbind, samples2))

w.as.MCMC <- secondary.samples[, paste0("w.as[",1:M,"]" )]
w.af.MCMC <- secondary.samples[, paste0("w.af[",1:NC,"]" )]


MCMC_nsamples = dim(primary.samples)[1]
phi.pred = matrix(0, nrow = MCMC_nsamples, N.p)
f.pred = matrix(0, nrow = MCMC_nsamples, N.p)


pb <- txtProgressBar(min = 1, max = N.p, style = 3)


for (i in 1:N.p) {
  
 
  Sys.sleep(1)
  
  # Print progress
  setTxtProgressBar(pb, i)
  
  
  #=====================================
  # Preamble for survival
  #=====================================
  
  # Add x0 to training data to find k nearest neighors
  coords.x0 = rbind(Neighbors$coords_sorted[, c("longitude", "latitude", "id")], coords.predict[i,c("X", "Y", "id" )]) # add point to training data
  
  # Sort coordinates by x axis 
  sorted_idx <- order(coords.x0[,1])
  x0.location = which(sorted_idx == M+1) # record location so easy to pick out necessary location
  coords_sorted =coords.x0[sorted_idx, ]
  #trdata = trdata[sorted_idx]
  
  coords.sf <- st_as_sf(data.frame(x = coords_sorted[, 1], y = coords_sorted[, 2]), 
                        coords = c('x', 'y'), crs = proj)
  
  # Compute distance between all points
  edist = matrix(as.numeric(st_distance(coords.sf)/1000), M+1, M+1)
  #x0.location
  
  if(x0.location <= (k+1)){
    
    neighbors_idx.x0 <- as.numeric(1:(x0.location -  1)) # k nearest neighbors of x0 - N0
    neighbors_dist.x0 <-  edist[x0.location, neighbors_idx.x0] # dist between x0 and its neighbors
    neighbors_dist.N0 <- edist[neighbors_idx.x0, neighbors_idx.x0 ] # distance between neighbors of x0
    
  } else {
    
    curr_sorted_idx = order(edist[x0.location, ])
    curr_sorted_idx = curr_sorted_idx[curr_sorted_idx < x0.location] # take values to the left
    neighbors_idx.x0 = curr_sorted_idx[1:k]
    neighbors_dist.x0 <-  edist[x0.location, neighbors_idx.x0] # dist between x0 and its neighbors
    neighbors_dist.N0 <- edist[neighbors_idx.x0, neighbors_idx.x0 ] # distance between neighbors of x0
    
  }
  
  #=====================================
  # Preamble for recruitment
  #=====================================
  
  coords.x0.f = rbind(Neighbors.C$coords_sorted[, c("longitude", "latitude", "id")], coords.predict[i, c("X", "Y", "id" )]) # add point to training data
  
  # Sort coordinates by x axis 
  sorted_idx.f <- order(coords.x0.f[,1])
  x0.location.f = which(sorted_idx.f == NC+1) # record location so easy to pick out necessary location
  coords_sorted.f =coords.x0.f[sorted_idx.f, ]

  
  coords.sf.f <- st_as_sf(data.frame(x = coords_sorted.f[, 1], y = coords_sorted.f[, 2]), 
                          coords = c('x', 'y'), crs = proj)
  
  # Compute distance between all points
  edist.f = matrix(as.numeric(st_distance(coords.sf.f)/1000), NC+1, NC+1)
  #x0.location
  
  if(x0.location.f <= (k+1)){
    
    neighbors_idx.x0.f <- as.numeric(1:(x0.location.f -  1)) # k nearest neighbors of x0 - N0
    neighbors_dist.x0.f <-  edist.f[x0.location.f, neighbors_idx.x0.f] # dist between x0 and its neighbors
    neighbors_dist.N0.f <- edist.f[neighbors_idx.x0.f, neighbors_idx.x0.f ] # distance between neighbors of x0
  
    
  } else {
    
    curr_sorted_idx.f = order(edist.f[x0.location.f, ])
    curr_sorted_idx.f = curr_sorted_idx.f[curr_sorted_idx.f < x0.location.f] # take values to the left
    neighbors_idx.x0.f = curr_sorted_idx.f[1:k]
    neighbors_dist.x0.f <-  edist.f[x0.location.f, neighbors_idx.x0.f] # dist between x0 and its neighbors
    neighbors_dist.N0.f <- edist.f[neighbors_idx.x0.f, neighbors_idx.x0.f ] # distance between neighbors of x0
   
    
  }
  
  for (mc in 1:MCMC_nsamples) {
    
    rho.as = primary.samples[mc,"rho.as"]
    alpha2.as =primary.samples[mc, "alpha2.as"]
    
    rho.af = primary.samples[mc,"rho.af"]
    alpha2.af =primary.samples[mc, "alpha2.af"]
    
    wst <-  mean(primary.samples[mc, grep("wst", colnames(primary.samples))])
    wft <-  mean(primary.samples[mc, grep("wft", colnames(primary.samples))])
    
    
    beta.f.tmax <- primary.samples[mc, "beta.f.tmax"]
    beta.f.tmaxlat <- primary.samples[mc, "beta.f.tmaxlat"]
    beta.f.tmaxlong <- primary.samples[mc, "beta.f.tmaxlong"]
    beta.f.lat <- primary.samples[mc, "beta.f.lat"]
    beta.f.long <- primary.samples[mc, "beta.f.long"]
    
    beta.s.tmax <- primary.samples[mc, "beta.s.tmax"]
    beta.s.tmaxlat <- primary.samples[mc, "beta.s.tmaxlat"]
    beta.s.tmaxlong <- primary.samples[mc, "beta.s.tmaxlong"]
    beta.s.lat <- primary.samples[mc, "beta.s.lat"]
    beta.s.long <- primary.samples[mc, "beta.s.long"]
    
    
    w.as.x0 <- w.as.MCMC[mc, neighbors_idx.x0] # This is spatial effects for x0 neighbours
    w.af.x0 <- w.af.MCMC[mc, neighbors_idx.x0.f] # This is spatial effects for x0 neighbours
    
    
    
    gp.pred.as <- NNGP.pred(neighbors_dist.x0, neighbors_dist.N0, x0.location, rho.as, alpha2.as, w.as.x0)
    phi.pred[mc, i] <- ilogit(wst + X_tmax[i]*beta.s.tmax + X_long[i]*beta.s.long + X_lat[i]*beta.s.lat + X_tmax_long[i]*beta.s.tmaxlong + X_tmax_lat[i]*beta.s.tmaxlat + gp.pred.as)
    
    
    gp.pred.af <- NNGP.pred(neighbors_dist.x0.f, neighbors_dist.N0.f, x0.location.f, rho.af, alpha2.af, w.af.x0)
    f.pred[mc,i ] <- exp( wft + X_tmax[i]*beta.f.tmax + X_long[i]*beta.f.long + X_lat[i]*beta.f.lat + X_tmax_long[i]*beta.f.tmaxlong + X_tmax_lat[i]*beta.f.tmaxlat + gp.pred.af)
    
  }
  
  
}

head(phi.pred)
head(f.pred)
result.predict = list(phi.pred = phi.pred, f.pred = f.pred, coords.sf.predict = coords.sf.predict)



#=============================================================================================================================================================================
# R code to produce spatially smooth prediction maps 
#=============================================================================================================================================================================


#=======================================================
# Survival
#=======================================================


phi.mean <-  apply(result.predict$phi.pred, 2, mean)
phi.sd <- apply(result.predict$phi.pred, 2, sd)

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
# Recruitment
#=======================================================



f.mean <-  apply(result.predict$f.pred, 2, mean)
f.sd <- apply(result.predict$f.pred, 2, sd)

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
  scale_fill_gradient2(low = "red", mid = "grey90", high = "blue", midpoint = mean(f.pred.df$Mean, na.rm = TRUE), na.value = NA, limits = range(f.pred.df$Mean, na.rm = TRUE)  )  + 
  labs(x = 'Longitude', y = 'Latitude')  + theme(legend.text = element_text(size = 10),  legend.title = element_text(size = 14))


ggplot() +
  geom_stars(data = raster_data_f, aes(fill = Mean)) +
  coord_sf() +
  theme_minimal() +
  geom_sf(data = st_cast(spatial_data, "MULTILINESTRING")) +
  scale_fill_gradient2(
    low = "red",
    mid = "grey90",
    high = "blue",
    midpoint = mean(f.pred.df$Mean, na.rm = TRUE),
    limits = range(f.pred.df$Mean, na.rm = TRUE),
    breaks = scales::pretty_breaks(n = 3),
    labels = scales::label_number(accuracy = 0.01),
    na.value = NA
  )+
  labs(x = 'Longitude', y = 'Latitude')  + theme(legend.text = element_text(size = 10),  legend.title = element_text(size = 14))


# Sd

ggplot() +
  geom_stars(data = raster_data_f, aes(fill = Sd)) +
  geom_sf(data = st_cast(spatial_data, "MULTILINESTRING")) +
  scale_fill_gradient2(low = "red", mid = "grey90", high = "blue", midpoint = mean(f.pred.df$Sd, na.rm = TRUE), na.value = NA, limits = range(f.pred.df$Sd, na.rm = TRUE)) + labs(x = 'Longitude', y = 'Latitude') +
  coord_sf() +
  theme_minimal() +theme(legend.text = element_text(size = 10),  legend.title = element_text(size = 14))


