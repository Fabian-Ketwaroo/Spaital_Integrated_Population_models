# Kriging predictions of survival and recruitment, using the BayesNSGP package
# ---------------------------------------------------------------------------
# Uses NNGP.pred() and computeNeighbors() from BayesNSGP (>= the version with the
# proj / dist.scale arguments and draw-wise prediction). Distances are in km, as in
# GRCA_predictions_phi_f, so that rho from the fitted model is on the same scale.

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
  library(BayesNSGP)
}

# (The standalone NNGP functions of the original analysis are no longer needed:
#  computeNeighbours -> BayesNSGP::computeNeighbors, DisplayNeighbors -> displayNeighbors,
#  NNGP.pred (one draw, hand-written neighbor search) -> BayesNSGP::NNGP.pred.)

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
Neighbors = computeNeighbors(coords.raw, k, proj)
str(Neighbors)


# plot to show neighbors
# red is ind location and blue are the neighbors
displayNeighbors(ind = 501, coords.sorted = Neighbors$coords_sorted, neighbor_idx = Neighbors$neighbor_idx, proj)
# This function would also be better if I set the spatial outline as an argument.


neighbors.id <- Neighbors$neighbor_idx

coords.sorted <- Neighbors$coords_sorted[, 1:2]
edist = Neighbors$edist_sorted # this is the distance matrix of the ordered coordinates. 
hist(edist)
nid.dist =Neighbors$neighbors_dist# the last column is the distance to itself
sorted.id <- Neighbors$coords_sorted[, 3:4]


k = 15
# Selecting k neighbors based of Vecchia's choice. Ordered along x axis.  
Neighbors.C = computeNeighbors(coords.raw[1:NC,], k, proj)
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
# Posterior draws from NNGP_BayesNSGP.R
#==========================================================================================================

res      <- readRDS("Case study/sipm_RW_NN_GP.rds")
samples  <- res$samples     # hyperparameters, regression coefficients, ...  (thinned by `thin`)
samples2 <- res$samples2    # spatial effects w.as, w.af, ...                (thinned by `thin2`)

# Kriging needs rho, alpha2 and the spatial effects from the SAME MCMC iteration.
# - If rho / alpha2 / coefficients were monitored in samples2 as well, use them from there.
# - Otherwise the two sets have different thinning, so match rows by iteration number.
need <- c("rho.as", "alpha2.as", "rho.af", "alpha2.af", "beta.f.tmax", "beta.s.tmax")
get_chain <- function(x, c) as.matrix(x[[c]])
nchain <- length(samples2)

if (all(need %in% colnames(get_chain(samples2, 1)))) {
  
  both <- lapply(seq_len(nchain), function(c) get_chain(samples2, c))
  primary.samples   <- do.call(rbind, both)
  secondary.samples <- primary.samples
  
} else {
  
  st <- if (!is.null(res$settings)) res$settings else list(nburnin = 200000,thin = 30, thin2 = 100)
  # saved iteration number of row r: nburnin + thin * r  (nburnin is a multiple of both thins here)
  pr <- sec <- vector("list", nchain)
  for (c in seq_len(nchain)) {
    a <- get_chain(samples, c); b <- get_chain(samples2, c)
    it1 <- st$nburnin + st$thin  * seq_len(nrow(a))
    it2 <- st$nburnin + st$thin2 * seq_len(nrow(b))
    common <- intersect(it1, it2)
    pr[[c]]  <- a[match(common, it1), , drop = FALSE]
    sec[[c]] <- b[match(common, it2), , drop = FALSE]
  }
  primary.samples   <- do.call(rbind, pr)
  secondary.samples <- do.call(rbind, sec)
}

stopifnot(nrow(primary.samples) == nrow(secondary.samples))
MCMC_nsamples <- nrow(primary.samples)
MCMC_nsamples   # number of posterior draws used for prediction

# Optional: use only the second half of each chain
# keep <- (MCMC_nsamples %/% 2 + 1):MCMC_nsamples
# primary.samples <- primary.samples[keep, ]; secondary.samples <- secondary.samples[keep, ]

w.as.MCMC <- secondary.samples[, paste0("w.as[", 1:M,  "]")]   # in the sorted order used for fitting
w.af.MCMC <- secondary.samples[, paste0("w.af[", 1:NC, "]")]

# Draw-wise quantities (vectors of length MCMC_nsamples)
col <- function(nm) primary.samples[, nm]
wst <- rowMeans(primary.samples[, grep("^wst\\[", colnames(primary.samples))])
wft <- rowMeans(primary.samples[, grep("^wft\\[", colnames(primary.samples))])

# Training locations, in the x-sorted order the model used (longitude, latitude only)
sites.as <- as.matrix(Neighbors$coords_sorted[,   c("longitude", "latitude")])
sites.af <- as.matrix(Neighbors.C$coords_sorted[, c("longitude", "latitude")])

phi.pred = matrix(0, nrow = MCMC_nsamples, N.p)
f.pred   = matrix(0, nrow = MCMC_nsamples, N.p)

pb <- txtProgressBar(min = 1, max = N.p, style = 3)

for (i in 1:N.p) {
  
  setTxtProgressBar(pb, i)
  x0 <- coords.predict[i, c("X", "Y")]
  
  # Spatial effect at the new location, one value per posterior draw.
  # proj + dist.scale = 1000 reproduce the distances used for fitting (km on the ellipsoid).
  gp.pred.as <- NNGP.pred(x0, sites.as, rho = col("rho.as"), sigma2 = col("alpha2.as"),
                          w.x0.all = w.as.MCMC, k = k, proj = proj, dist.scale = 1000)
  gp.pred.af <- NNGP.pred(x0, sites.af, rho = col("rho.af"), sigma2 = col("alpha2.af"),
                          w.x0.all = w.af.MCMC, k = k, proj = proj, dist.scale = 1000)
  
  phi.pred[, i] <- plogis(wst + X_tmax[i]*col("beta.s.tmax") + X_long[i]*col("beta.s.long") + X_lat[i]*col("beta.s.lat") +
                            X_tmax_long[i]*col("beta.s.tmaxlong") + X_tmax_lat[i]*col("beta.s.tmaxlat") + gp.pred.as)
  
  f.pred[, i] <- exp(wft + X_tmax[i]*col("beta.f.tmax") + X_long[i]*col("beta.f.long") + X_lat[i]*col("beta.f.lat") +
                       X_tmax_long[i]*col("beta.f.tmaxlong") + X_tmax_lat[i]*col("beta.f.tmaxlat") + gp.pred.af)
}
close(pb)

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

f.pred.sf$x <- st_coordinates(f.pred.sf)[,1] 
f.pred.sf$y <- st_coordinates(f.pred.sf)[,2] 

raster_data_f <- st_as_stars(f.pred.sf, dims = c("x", "y"))


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
