

# Load packages
{
  #rm(list = ls())
  library(nimble)
  library(sf)
  library(MCMCvis)
  library(coda)
  library(MASS)
  library(parallel)
  library(ggplot2)
}



# ==============================================================================
# NOTE ON COMPUTATIONAL WORKFLOW & MEMORY MANAGEMENT
# ==============================================================================
# Due to the high-dimensional nature of the spatial integrated population model 
# and extensive memory tracking within NIMBLE, running the full MCMC chains 
# and generating the Goodness-of-Fit (GOF) metrics simultaneously in a single 
# session can cause R to exhaust system memory and crash.
#
# WORKAROUND: We use a modular, two-step approach:
#   1. Run the primary model and save the MCMC posterior outputs to disk.
#   2. (This Script) The entire model is re-run to save GOF metrics needed for GOF analysis. 
#
# This keeps the memory footprint manageable and ensures execution stability.
# ==============================================================================

# Working directory should be the project root. Open the .Rproj file, or use
# setwd() to point to the folder containing "Case Study" and "Simulation study".

# Load data 
# The 'GRCA_data_2004_2014.Rdata' file is included in the 'Case study' directory.
load("Case Study/GRCA_data_2004_2014.Rdata")
str(GRCA_data)


# Number of count locations 
NC = dim(GRCA_data$C)[1] 

# Number of Capture-recapture location
NCR =  dim(GRCA_data$Rm)[1] 

# Total number of locations
M <- NC + NCR

# Number of sampling occasions 
n.occasions =  dim(GRCA_data$C)[2] 


# Get spatial co-ordinates of sampling locations: 
coords = GRCA_data$GRCA_splocations[, -3]  #GRCA_10yrs$GRCA_splocations
source = as.factor(c(rep("C", NC), rep("CR", NCR) ))
coords = cbind(coords, source)
coords.raw <- cbind(coords, id = 1:M)
proj = "+proj=longlat +datum=WGS84" 

# Load functions needed for the nearest-neighbor Gaussian process (NNGP) 
source("Case Study/NNGP_functions.R")


# Number of nearest-neighbor
k = 15

# We employ NNGP on survival and recruitment. Survival has NNGP terms for all locations while
# recruitment only has NNGHP terms at count locations

## NNGP preamble for survival 

# Selecting k neighbors based of Vecchia's choice. Ordered along x axis.  
Neighbors = computeNeighbours(coords.raw, k, proj)
str(Neighbors)


# plot to show Neighbors
# red is ind location and blue are the neighbors
DisplayNeighbors(ind = 501, coords.sorted = Neighbors$coords_sorted, neighbor_idx = Neighbors$neighbor_idx, proj)

# Neighbors id 
neighbors.id <- Neighbors$neighbor_idx

# Sorted co-ordinates 
coords.sorted <- Neighbors$coords_sorted[, 1:2]

# Sorted euclidean distance of all co-ordinates
edist = Neighbors$edist_sorted 

# Euclidean distance between k nearest-neighbors. The last column is the distance to itself
nid.dist =Neighbors$neighbors_dist

# ID to match sorted locations and data
sorted.id <- Neighbors$coords_sorted[, 3:4]


## NNGP preamble for recruitment

# Selecting k neighbors based of Vecchia's choice. Ordered along x axis. 
# Only on count locations 
Neighbors.C = computeNeighbours(coords.raw[1:NC,], k, proj)

# Neighbors id 
neighbors.Cid <- Neighbors.C$neighbor_idx

# Sorted count co-ordinates 
coords.sorted.C <- Neighbors.C$coords_sorted[, 1:2]

# Sorted euclidean distance of count co-ordinates
edist.C = Neighbors.C$edist_sorted # this is the distance matrix of the ordered coordinates. 

# Euclidean distance between k nearest-neighbors. The last column is the distance to itself
nid.Cdist =Neighbors.C$neighbors_dist# the last column is the distance to itself

# ID to match sorted locations and data
sorted.Cid <- Neighbors.C$coords_sorted[, 3:4]



# Residual maximum temperature (TMAX)
X_tmax = GRCA_data$Tmax.resid[, 1:(n.occasions-1) ]#

# Re-order residual TMAX to match NNGP ordering 
X_tmax <- X_tmax[sorted.id[,2],]

# Longitude 
X_long <- matrix(0, nrow = M, ncol = n.occasions-1)
X_long[, 1:(n.occasions-1)] <- coords$longitude[sorted.id[,2]]

# Latitude 
X_lat <- matrix(0, nrow = M, ncol = n.occasions-1)
X_lat[, 1:(n.occasions-1)] <- coords$latitude[sorted.id[,2]]

# residual TMAX*Longitude
X_tmax_long <- X_tmax * X_long

# residual TMAX*Longitude
X_tmax_lat <- X_tmax * X_lat

# Scale values 
X_tmax <- scale(X_tmax)
X_long <- scale(X_long)
X_lat <- scale(X_lat)
X_tmax_lat <- scale(X_tmax_lat)
X_tmax_long <- scale(X_tmax_long)



# Re-ordering BBS counts to match NNGP ordering of spatial locations
C.sorted <- which(sorted.id[,1]== "C") # indicator for which sites are counts sites in sorted data
C.id <- sorted.id[,2][which(sorted.id[,1] == "C")]  # indicator for which sites are counts sites in orignial data but also matches order in sorted data
C.all <- matrix(NA, M, n.occasions)
C.all[C.sorted, ] = GRCA_data$C[C.id, ]
rownames(C.all) = rownames(coords.sorted)

# remove missing values 
C = C.all[!is.na(C.all)]

# Number of counts 
n.counts = length(C)

# Converting observers identity (obs) and 1st time observers (obs1) to scalar form 
site = c()
yr = c()

for (n in 1:n.occasions) {
  
  a = which(!is.na(C.all[,n]))
  b = rep(n, length(a))
  site = c(site, a)
  yr = c(yr, b)
  
}

# observers identity (obs)
obs.c = GRCA_data$obs
obs.all =  matrix(NA, M, n.occasions)
obs.all[C.sorted, ] <- obs.c[C.id,]
obs = obs.all[!is.na(C.all)]

#1st time observers (obs1)
oby1 = GRCA_data$oby1
oby1.all <- matrix(NA, M, n.occasions)
oby1.all[C.sorted, ] <- oby1[C.id,]
oby1 = oby1.all[!is.na(C.all)]



# Re-ordering MAPS counts to match NNGP ordering of spatial locations 
CR.sorted <- which(sorted.id[,1]=="CR") # indicator for which sites are CR sites in sorted data
CR.id <- sorted.id[,2][which(sorted.id[,1] =="CR")]  # indicator for which sites are CR sites in original data but also matches order in sorted data

# Residents m-array
Rmaps = GRCA_data$Rm

Rm <- array(0, dim = c(M, n.occasions-1, n.occasions ))
Rm[ CR.sorted, ,] <- Rmaps[CR.id-NC, ,]

Fm <- array(0, dim = c(M, n.occasions, n.occasions+1 ))
Fm[ CR.sorted , ,] <- GRCA_data$Fm[CR.id-NC, ,] # Mixture of possible transient and residents individuals m-array


# Number of "releases"
rel.R <- matrix(0, nrow = M, ncol = n.occasions-1)
rel.F <- matrix(0, nrow = M, ncol = n.occasions)

for (i in CR.sorted ) {
  rel.R[i, ] <- rowSums(Rm[i,,])
  rel.F[i, ] <- rowSums(Fm[i,,])
}


# visits- binary indicator whether site active or not 
visits <- GRCA_data$visits
visits.all <- matrix(0, M, n.occasions-1)
visits.all[CR.sorted, ] <- visits[CR.id-NC,]
visits.all[CR.sorted[15], c(6)] <- 1


# List of constants
win.data = list(rel.F= rel.F,rel.R= rel.R, C.sorted = C.sorted, CR.sorted = CR.sorted, n.occasions=n.occasions, M=M, NCR = NCR, NC= NC, mu = rep(0, M), 
                edist = edist, dmax = max(edist), dmin = min(edist[edist!= 0] ), dq75 = quantile(edist)[4], nobs.id = sort((unique(obs.all[!is.na(obs.all)]))), nobs  = length((unique(obs.all[!is.na(obs.all)]))),
                nid.dist = nid.dist, n.id = neighbors.id, k = k,nid.Cdist = nid.Cdist, n.Cid = neighbors.Cid, edist.C = edist.C,
                visits = visits.all, n.counts = n.counts, obs = obs, first = oby1, site = site, yr = yr, dq90 = quantile(edist, 0.90))
str(win.data)



sipm <- nimbleCode({
  
  
  # Priors
  for (i in 1:(n.occasions-1)) {
    rf[i] ~ dunif(0,1) # rf[t] ... residency probability (within first captured individuals); probability that animal captured first on occasion t is a resident	
    w[i] ~ dunif(0,1) # Probability of the confirmation of residency status for the newly captured resident (e.g. by a recapture within the occasion of the first capture
    
  }
  
  rf[n.occasions] <- rf[n.occasions-1]
  w[n.occasions] <- w[n.occasions-1]
  
  
  # Fixed effects 
  beta.f.tmax ~ dnorm(0, var = 2) 
  beta.s.tmax ~ dnorm(0, var = 2) 
  
  beta.s.long ~ dnorm(0, var = 2) 
  beta.f.long ~ dnorm(0, var = 2) 
  
  beta.s.lat ~ dnorm(0, var = 2) 
  beta.f.lat ~ dnorm(0, var = 2) 
  
  beta.s.tmaxlong ~ dnorm(0, var = 2) 
  beta.f.tmaxlong ~ dnorm(0, var = 2) 
  
  beta.s.tmaxlat ~ dnorm(0, var = 2) 
  beta.f.tmaxlat ~ dnorm(0, var = 2) 
  
  # First time observer observer effect 
  eta ~ dnorm(0, var = 2) 
  
  # Observer random effect
  for (i in 1:nobs) {
    omega[nobs.id[i]] ~ dnorm(0, var = 2) 
  }
  
  
  ################### 
  # GP's
  ###################
  
  # Recruitment
  alpha2.af ~ dinvgamma(3,1) # signal variance 
  rho.af ~ dunif(dmin/3, dmax/3) # length scale 
  
  AD.af[1:(NC), 1:(k+1)] <- computeAD(edist.C[1:(NC),1:(NC)], nid.Cdist[1:(NC), 1:(k+1)], n.Cid[1:(NC), 1:k], rho.af, alpha2.af, k )
  w.af[1:(NC)] ~ dmnorm_NNGP(mu[1:(NC)], AD.af[1:( NC), 1:(k+1)], n.Cid[1:(NC), 1:k])
  constraint.f ~ dnorm(sum( w.af[1:NC]), sd = 1) # sum to zero soft constraint 
  
  
  # Survival
  alpha2.as ~ dinvgamma(3,1)  # signal variance 
  rho.as ~ dunif(dmin/3, dmax/3) # length scale 
  
  
  AD.as[1:(M), 1:(k+1)] <- computeAD(edist[1:(M),1:(M)], nid.dist[1:(M), 1:(k+1)], n.id[1:(M), 1:k], rho.as, alpha2.as, k )
  w.as[1:(M)] ~ dmnorm_NNGP(mu[1:(M)], AD.as[1:(M), 1:(k+1)], n.id[1:(M), 1:k])
  constraint.s ~ dnorm(sum( w.as[1:M]), sd = 1) # sum to zero soft constraint 
  
  
  # Temporal intercepts 
  for (t in 1:(n.occasions-1)) {
    wst[t] ~ dnorm(0, var = 2)
    wft[t] ~ dnorm(0, var = 2)
    alpha.p[t] ~ dnorm(0, var = 2)
  }
  
  
  
  # Latent processes   
  for (i in 1:M) {
    
    # Priors and linear models
    for (t in 1:(n.occasions-2)){
      logit(sa[i,t]) <- wst[t] + X_tmax[i,t]*beta.s.tmax + X_long[i,t]*beta.s.long + X_lat[i,t]*beta.s.lat + X_tmax_long[i,t]*beta.s.tmaxlong + X_tmax_lat[i,t]*beta.s.tmaxlat + w.as[i] 
      p[i,t] <- ilogit(alpha.p[t] + p.ref[i]) * visits[i,t]  
    }
    
    # remove GP for last year
    logit(sa[i,n.occasions-1]) <- wst[n.occasions-1] + X_tmax[i,n.occasions-1]*beta.s.tmax + X_long[i,n.occasions-1]*beta.s.long + X_lat[i,n.occasions-1]*beta.s.lat + X_tmax_long[i,n.occasions-1]*beta.s.tmaxlong + X_tmax_lat[i,n.occasions-1]*beta.s.tmaxlat
    p[i,n.occasions-1] <- ilogit(alpha.p[n.occasions-1] + p.ref[i]) * visits[i,n.occasions-1]  
    
  }
  
  for (i in 1:NC) {
    
    u[C.sorted[i]] ~ dunif(1,100)
    N[C.sorted[i],1] <- round(u[C.sorted[i]])
    
    # Process model over time: our model of population dynamics
    for (t in 1:(n.occasions-2)){
      
      log(f[C.sorted[i],t]) <- wft[t] + X_tmax[C.sorted[i],t]*beta.f.tmax + X_long[C.sorted[i],t]*beta.f.long + X_lat[C.sorted[i],t]*beta.f.lat + X_tmax_long[C.sorted[i],t]*beta.f.tmaxlong + X_tmax_lat[C.sorted[i],t]*beta.f.tmaxlat  + w.af[i]
      N[C.sorted[i],t+1] ~ dpois(sa[C.sorted[i],t]*N[C.sorted[i],t] + f[C.sorted[i],t]*N[C.sorted[i],t])
      
      # Growth rate
      growth.rate[C.sorted[i],t] <- (N[C.sorted[i], t+1]) / (N[C.sorted[i],t] + 0.001)
    }
    
    
    log(f[C.sorted[i],n.occasions-1]) <- wft[n.occasions-1] + X_tmax[C.sorted[i],n.occasions-1]*beta.f.tmax + X_long[C.sorted[i],n.occasions-1]*beta.f.long + X_lat[C.sorted[i],n.occasions-1]*beta.f.lat + X_tmax_long[C.sorted[i],n.occasions-1]*beta.f.tmaxlong + X_tmax_lat[C.sorted[i],n.occasions-1]*beta.f.tmaxlat
    N[C.sorted[i],n.occasions] ~ dpois(sa[C.sorted[i],n.occasions-1]*N[C.sorted[i],n.occasions-1] + f[C.sorted[i],n.occasions-1]*N[C.sorted[i],n.occasions-1])
    growth.rate[C.sorted[i],n.occasions-1] <- (N[C.sorted[i], n.occasions]) / (N[C.sorted[i],n.occasions-1] + 0.001)
    
    
    # recapture site random effect 
    p.ref[C.sorted[i]] <- 0
  }
  
  # Count likelihood 
  for (j in 1:n.counts) {
    
    log(lambda[j]) <- log(N[site[j], yr[j]]) + eta*first[j] + omega[obs[j]]
    C[j] ~ dpois(lambda[j])
    C.pred[j] ~ dpois(lambda[j])
  }
  
  # Capture-recapture  likelihood 
  for (i in 1:NCR) { # i in CR.id not supported
    
    # recapture site random effect 
    p.ref[CR.sorted[i]] ~ dnorm(0, var = 2.72)
    
    # CJS model with multinomial likelihood for mixture of residents and transients
    for (t in 1:n.occasions) {
      Fm[CR.sorted[i],t, 1:(n.occasions+1)] ~ dmulti(fm[CR.sorted[i], t,1:(n.occasions+1)], rel.F[CR.sorted[i],t])
      Fm.pred[CR.sorted[i],t, 1:(n.occasions+1)] ~ dmulti(fm[CR.sorted[i], t,1:(n.occasions+1)], rel.F[CR.sorted[i],t])
      
      
      fm.post[i, t,1:(n.occasions+1)] <- fm[CR.sorted[i], t,1:(n.occasions+1)]
      
      # Main diagonal
      fm[CR.sorted[i],t,t] <- rf[t]*w[t]
    }
    
    # Immediately above main diagonal (j = t+1)
    for (t in 1:(n.occasions-1)) {
      fm[CR.sorted[i],t,t+1] <- rf[t]*(1-w[t])*sa[CR.sorted[i],t]*p[CR.sorted[i],t]
    }
    
    # Further above main diagonal (j >= t+2)
    for (t in 1:(n.occasions-2)) {
      for (j in (t+2):n.occasions) {
        fm[CR.sorted[i],t,j] <- rf[t]*(1-w[t])*sa[CR.sorted[i],j-1]*p[CR.sorted[i],j-1] *
          prod(sa[CR.sorted[i],t:(j-2)]*(1-p[CR.sorted[i],t:(j-2)]))
      }
    }
    
    # Below main diagonal
    for (t in 2:n.occasions) {
      for (j in 1:(t-1)) {
        fm[CR.sorted[i],t,j] <- 0
      }
    }
    
    # Last column: probability of non-recapture after 1st capture
    for (t in 1:n.occasions) {
      fm[CR.sorted[i],t,n.occasions+1] <- 1 - sum(fm[CR.sorted[i], t, 1:n.occasions])
    }
    
    
    
    # Capture-recapture data (CJS model with multinomial likelihood) for residents
    for (t in 1:(n.occasions-1)){
      Rm[CR.sorted[i],t,1:n.occasions] ~ dmulti(r[CR.sorted[i],t,1:n.occasions], rel.R[CR.sorted[i],t])
      Rm.pred[CR.sorted[i],t,1:n.occasions] ~ dmulti(r[CR.sorted[i],t,1:n.occasions], rel.R[CR.sorted[i],t])
      
      rm.post[i,t,1:n.occasions]<- r[CR.sorted[i],t,1:n.occasions]
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


# Initial values 
{
  set.seed(8)
  
  inits <- list(beta.f.tmax = runif(1,-0.3, 0), beta.s.tmax = runif(1,-0.5, 0), beta.s.long = 0, beta.f.long = 0, beta.s.lat = 0, beta.f.lat =0, beta.f.tmaxlong = 0, beta.s.tmaxlong =0,
                beta.f.tmaxlat = 0, beta.s.tmaxlat = 0, eta = runif(1, -0.5, 0),  u = runif(win.data$M, 10,50), 
                rho.af = mean(edist.C)/3, alpha2.af = rinvgamma(1,3,1),  rho.as = mean(edist)/3, alpha2.as = rinvgamma(1,3,1),
                w = runif(n.occasions,0,1), rf = runif(n.occasions, 0,1), omega = rnorm(max(win.data$nobs.id) ), p.ref = rnorm(win.data$M), alpha.p = rep(0,n.occasions-1)  )
  
  
  
  inits$wst = nimble::logit(runif(n.occasions-1, 0.6,0.9) )#rnorm(n.occasions-1, 0, 1)
  inits$wft = log(runif(n.occasions-1,0.3,0.8))#rnorm(n.occasions-1, 0, 1)
  
  inits$w.af = MASS::mvrnorm(1, mu =rep(0,NC), Sigma = expcov(edist.C[1:(NC),1:(NC)], inits$rho.af, inits$alpha2.af))
  inits$w.as = MASS::mvrnorm(1, mu =rep(0,M), Sigma = expcov(edist[1:(M),1:(M)], inits$rho.as, inits$alpha2.as))
  
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
  
  inits <- list(beta.f.tmax = runif(1,-0.3, 0), beta.s.tmax = runif(1,-0.5, 0), beta.s.long = 0, beta.f.long = 0, beta.s.lat = 0, beta.f.lat =0, beta.f.tmaxlong = 0, beta.s.tmaxlong =0,
                beta.f.tmaxlat = 0, beta.s.tmaxlat = 0, eta = runif(1, -0.5, 0),  u = runif(win.data$M, 10,50), 
                rho.af = mean(edist.C)/3, alpha2.af = rinvgamma(1,3,1),  rho.as = mean(edist)/3, alpha2.as = rinvgamma(1,3,1),
                w = runif(n.occasions,0,1), rf = runif(n.occasions, 0,1), omega = rnorm(max(win.data$nobs.id) ), p.ref = rnorm(win.data$M), alpha.p = rep(0,n.occasions-1)  )
  
  
  inits$wst = nimble::logit(runif(n.occasions-1, 0.6,0.9) )#rnorm(n.occasions-1, 0, 1)
  inits$wft = log(runif(n.occasions-1,0.3,0.8))#rnorm(n.occasions-1, 0, 1)
  
  inits$w.af = MASS::mvrnorm(1, mu =rep(0,NC), Sigma = expcov(edist.C[1:(NC),1:(NC)], inits$rho.af, inits$alpha2.af))
  inits$w.as = MASS::mvrnorm(1, mu =rep(0,M), Sigma = expcov(edist[1:(M),1:(M)], inits$rho.as, inits$alpha2.as))
  
  
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

# Build model
FullModel = nimbleModel(sipm, data = list(Rm= Rm, Fm= Fm, C = C, X_tmax = X_tmax, X_long = X_long, X_lat = X_lat, X_tmax_long = X_tmax_long, X_tmax_lat = X_tmax_lat, constraint.s = 0, constraint.f = 0 ), constants = win.data,inits = inits1 , calculate = FALSE)
# [Warning] Detected backwards indexing in `1:0 index warnings are harmless and can be safely ignored.
# [Warning] dimensions specified are larger than model specification for variable `Rm` is harmless and can be safely ignored
# [Warning] dimensions specified are larger than model specification for variable `Fm` is harmless and can be safely ignored


# Configure MCMC
FullModelconf = configureMCMC(FullModel)  

# Reset Monitors 
FullModelconf$resetMonitors()
FullModelconf$addMonitors( "C.pred",   "Fm.pred", "Rm.pred", "lambda", "fm.post", "rm.post" )

# RW Block sampling of fixed effects 
FullModelconf$removeSamplers("beta.f.tmax", "beta.s.tmax", "beta.f.tmaxlong", "beta.s.tmaxlong", "beta.f.tmaxlat", "beta.s.tmaxlat" )
beta.nodes <- c("beta.f.tmax", "beta.s.tmax", "beta.f.tmaxlong", "beta.s.tmaxlong", "beta.f.tmaxlat", "beta.s.tmaxlat" )
FullModelconf$addSampler(target = beta.nodes, type = "RW_block")


# Individual Meterpolis Hastings random walk sampking of NNGP spatial terms
FullModelconf$removeSamplers("w.af", "w.as")

for (i in 1:(NC)) {
  FullModelconf$addSampler(target = paste0("w.af[", i, "]"), type = "RW")
}

for (i in 1:(M) ) {
  FullModelconf$addSampler(target = paste0("w.as[", i, "]"), type = "RW")
}

# Check samplers are correctly assigned 
FullModelconf

# Build MCMC
FullModelmcmc = buildMCMC(FullModelconf)

# Compile nimble model 
cFullModel <- compileNimble(FullModel)

# Compile MCMC 
FullModelmod = compileNimble(FullModelmcmc, project = FullModel, resetFunctions = TRUE)

# Run MCMC
NNGP_time  <- system.time(FullModelmcmc.out <- runMCMC(FullModelmod,niter = 350000, nburnin = 200000,thin = 30, thin2 = 100,  samplesAsCodaMCMC = TRUE, nchains = 2, inits = list(inits1, inits2), summary = TRUE)) 

# MCMC samples 
samples <- FullModelmcmc.out$samples 


# Combine all chains to one 
samples = mcmc(do.call(rbind, samples) )
L = dim(samples)[1]

# Let's extract posterior 
C.pred <- samples[, paste0("C.pred[",1 :n.counts, "]")]
Rm.pred <- samples[, grep("Rm.pred", colnames(samples))]
Fm.pred <- samples[, grep("Fm.pred", colnames(samples))]
lambda.post <- samples[, grep("lambda", colnames(samples))]
fm.post <- samples[, grep("fm.post", colnames(samples))]
rm.post <- samples[, grep("rm.post", colnames(samples))]


#================================================
# Count data  
#===============================================


# Year grouping
# Grouping the data across year may help reveal whether the model fails to adequately represent variation in time
# Generally, the more ways we explore the fit of our model to the data, the more confidence we have that our model is adequately representing the data and we can draw reasonable conclusions from it.

Y = n.occasions
FT.C <- FT.C.pred <- FT.BPV <- matrix(0, L, Y)
X2.C <- X2.C.pred <- X2.BPV <- matrix(0, L, Y)

for (y in 1:Y) {
  
  # Data at each year 
  C.y <-   sum(C[yr == y])
  
  for (l in 1:L) {
    
    # Data at each year and each MCMC iteration
    lambda.y <- sum(lambda.post[l, yr ==y])
    C.pred.y <- sum(C.pred[l, yr ==y])
    
    # Freeman- Tukey statistics  
    FT.C[l, y] <-  ((sqrt(C.y) - sqrt(lambda.y))^2)
    FT.C.pred[l, y] <-  ((sqrt(C.pred.y) - sqrt(lambda.y))^2)
    
    FT.BPV[l, y] = FT.C.pred[l, y]  >  FT.C[l, y] 
    
    # Chi- square statistis
    X2.C[l,y] <- sum((C.y- lambda.y)^2/( lambda.y + 0.0001 ))
    X2.C.pred[l,y] <- sum((C.pred.y- lambda.y)^2/( lambda.y + 0.0001 ))
    X2.BPV[l,y] <- X2.C.pred[l,y] > X2.C[l, y]
    
  }
  
  
}


# per-year Bayesian p-values
colMeans(FT.BPV)
colMeans(X2.BPV)

plot(colMeans(FT.BPV), xlab = "Year", ylab = "Bayesian p-value")
plot(colMeans(X2.BPV), xlab = "Years", ylab = "X2.BPV")

# overall year-grouped Bayesian p-value
FT.BPV2 <- numeric(L)

for (l in 1:L) {
  FT.BPV2[l] <- sum(FT.C.pred[l,1:Y ]) > sum(FT.C[l,1:Y ])
}

mean(FT.BPV2)

# remove second year
FT.BPV2 <- X2.BPV2 <- numeric(L)

for (l in 1:L) {
  FT.BPV2[l] <- sum(FT.C.pred[l,2:Y ]) > sum(FT.C[l,2:Y ])
  X2.BPV2[l] <- sum(X2.C.pred[l,2:Y ]) > sum(X2.C[l,2:Y ])
}

mean(FT.BPV2) 
mean(X2.BPV2)


#==========================================================================
# Site grouping
# Grouping the data across sites may help reveal whether the model fails to adequately represent variation in space

Y = NC
FT.C <- FT.C.pred <- FT.BPV <- matrix(0, L, Y)
X2.C <- X2.C.pred <- X2.BPV <- matrix(0, L, Y)

site.id <- sort(unique(site))

for (y in 1:Y) {
  
  print(y)
  
  # Data at each year 
  C.y <-   sum(C[which(site == site.id[y])])
  
  for (l in 1:L) {
    
    # Data at each year and each MCMC iteration
    lambda.y <- sum(lambda.post[l, site == site.id[y] ])
    C.pred.y <- sum(C.pred[l, site == site.id[y]])
    
    # Freeman- Tukey statistics  
    FT.C[l, y] <-  ((sqrt(C.y) - sqrt(lambda.y))^2)
    FT.C.pred[l, y] <-  ((sqrt(C.pred.y) - sqrt(lambda.y))^2)
    
    FT.BPV[l, y] = FT.C.pred[l, y]  >  FT.C[l, y] 
    
    # Chi- square statisti
    X2.C[l,y] <- ((C.y- lambda.y)^2/( lambda.y + 0.0001 ))
    X2.C.pred[l,y] <- ((C.pred.y- lambda.y)^2/( lambda.y + 0.0001 ))
    X2.BPV[l,y] <- X2.C.pred[l,y] > X2.C[l, y]
    
  }
  
  
}


colMeans(FT.BPV)
colMeans(X2.BPV)

plot(colMeans(FT.BPV), xlab = "Sampling location", ylab = "Bayesian p-value")
plot(colMeans(X2.BPV), xlab = "Site", ylab = "X2.BPV")

plot( colMeans(FT.C)  - colMeans(FT.C.pred))

# overall site-grouped Bayesian p-value
FT.BPV2 <- X2.BPV2 <- numeric(L)

for (l in 1:L) {
  FT.BPV2[l] <- sum(FT.C.pred[l,1:Y ]) > sum(FT.C[l,1:Y ])
  X2.BPV2[l] <- sum(X2.C.pred[l,1:Y ]) > sum(X2.C[l,1:Y ])
}

mean(FT.BPV2) 
mean(X2.BPV2)


ppc.df <- data.frame(fit = rowSums(FT.C[,1:Y ]), 
                     fit.rep = rowSums(FT.C.pred[,1:Y ]), 
                     color = 'lightskyblue1')
ppc.df$color[ppc.df$fit.rep > ppc.df$fit] <- 'lightsalmon'
plot(ppc.df$fit, ppc.df$fit.rep, bg = ppc.df$color, pch = 21, 
     ylab = 'Fit', xlab = 'True', ylim = range(ppc.df$fit.rep, ppc.df$fit), main = " "  )
lines(ppc.df$fit, ppc.df$fit, col = 'black')

ppc.df <- data.frame(fit = rowSums(X2.C[,1:Y ]), 
                     fit.rep = rowSums(X2.C.pred[,1:Y ]), 
                     color = 'lightskyblue1')
ppc.df$color[ppc.df$fit.rep > ppc.df$fit] <- 'lightsalmon'
plot(ppc.df$fit, ppc.df$fit.rep, bg = ppc.df$color, pch = 21, 
     ylab = 'Fit', xlab = 'True', ylim = range(ppc.df$fit.rep, ppc.df$fit), main = "X2"  )
lines(ppc.df$fit, ppc.df$fit, col = 'black')



#================================================
# CR data 
#===============================================



# Let's extract posterior 
Rm.pred <- samples[, grep("Rm.pred", colnames(samples))]
Fm.pred <- samples[, grep("Fm.pred", colnames(samples))]
fm.post <- samples[, grep("fm.post", colnames(samples))]
rm.post <- samples[, grep("rm.post", colnames(samples))]


#============================================================
# Start with Rm 
#============================================================

#============================================================
# Temporal grouping
#============================================================


# Temporal grouping

FT.CR <- array(0, dim = c(L, (n.occasions-1), n.occasions))
FT.CR.pred <- array(0, dim = c(L, (n.occasions-1), n.occasions  ))
FT.BPV <- numeric(L)

X2.CR <- numeric(L)
X2.CR.pred <- numeric(L)
X2.BPV <- numeric(L)

for (t in 1:(n.occasions-1)) {
  
  print(t)
  
  for (t2 in 1:(n.occasions)) {
    
    CR <-  sum(Rm[CR.sorted, t, t2] ) 
    id.r <- paste0("rm.post[",1:NCR,", ", t, ", ", t2, "]"  )
    id.R <- paste0("Rm.pred[",CR.sorted,", ", t, ", ", t2, "]"  )
    
    for (l in 1:L) {
      
      lambda  <- sum(rm.post[l, id.r]*rel.R[CR.sorted,t])
      
      FT.CR[l,t, t2] <- sum((sqrt(CR) - sqrt(lambda))^2)
      
      CR.pred <- sum(Rm.pred[l,id.R]) #Rm[CR.sorted[i], t, ]  
      
      FT.CR.pred[l, t, t2] <- sum((sqrt(CR.pred) - sqrt(lambda))^2)
      
    }
    
  }
  
  
}



for (l in 1:L) {
  FT.BPV[l] = sum(FT.CR.pred[l, ,]) > sum(FT.CR[l, ,])
}

mean(FT.BPV)

mean(apply(FT.CR.pred, 1, sum) > apply(FT.CR, 1, sum))

FT.BPV_T = matrix(0, L, n.occasions-1 )

for (l in 1:L) {
  for (t in 1:(n.occasions-1)) {
    FT.BPV_T[l, t] = sum(FT.CR.pred[l,t ,]) > sum(FT.CR[l,t ,]) 
  }
}


plot(colMeans(FT.BPV_T), ylab = "Bayesian p-value", xlab = "Year", main = " ")

ppc.df <- data.frame(fit = apply(FT.CR, 1, sum), 
                     fit.rep = apply(FT.CR.pred, 1, sum), 
                     color = 'lightskyblue1')
ppc.df$color[ppc.df$fit.rep > ppc.df$fit] <- 'lightsalmon'
plot(ppc.df$fit, ppc.df$fit.rep, bg = ppc.df$color, pch = 21, 
     ylab = 'Fit', xlab = 'True', ylim = range(ppc.df$fit.rep, ppc.df$fit), main = "FT"  )
lines(ppc.df$fit, ppc.df$fit, col = 'black')


#============================================================
# # Spatial grouping
#============================================================


FT.CR <- array(0, dim = c(L, NCR, n.occasions))
FT.CR.pred <- array(0, dim = c(L, NCR, n.occasions  ))
FT.BPV <- numeric(L)

X2.CR <- numeric(L)
X2.CR.pred <- numeric(L)
X2.BPV <- numeric(L)

for (i in 1:NCR) {
  
  print(i)
  
  for (t2 in 1:(n.occasions-1)) {
    
    CR <-  (Rm[CR.sorted[i], t2, 1:n.occasions] ) 
    id.r <- paste0("rm.post[",i,", ", t2, ", ", 1:n.occasions, "]"  )
    id.R <- paste0("Rm.pred[",CR.sorted[i],", ", t2, ", ", 1:n.occasions, "]"  )
    
    for (l in 1:L) {
      
      lambda  <- (rm.post[l, id.r]*rel.R[CR.sorted[i],t2])
      
      FT.CR[l,i, t2] <- sum((sqrt(CR) - sqrt(lambda))^2)
      
      CR.pred <- (Rm.pred[l,id.R]) # Posterior predictive 
      
      FT.CR.pred[l, i, t2] <- sum((sqrt(CR.pred) - sqrt(lambda))^2)
      
    }
    
  }
  
  
}


for (l in 1:L) {
  FT.BPV[l] = sum(FT.CR.pred[l, ,]) > sum(FT.CR[l, ,])
}

mean(FT.BPV)

mean(apply(FT.CR.pred, 1, sum) > apply(FT.CR, 1, sum))

FT.BPV_T = diff = matrix(0, L, NCR )

for (l in 1:L) {
  for (i in 1:NCR) {
    FT.BPV_T[l, i] = sum(FT.CR.pred[l,i ,]) > sum(FT.CR[l,i ,])
    diff[l, i] = sum(FT.CR.pred[l,i ,]) - sum(FT.CR[l,i ,])
  }
}


colMeans(FT.BPV_T)
plot(colMeans(FT.BPV_T), xlab = "Sampling location", ylab = "Bayesian p-value", main = " ")

#============================================================
# Now FM
#============================================================

# Temporal grouping

FT.CR <- array(0, dim = c(L, (n.occasions-1), n.occasions))
FT.CR.pred <- array(0, dim = c(L, (n.occasions-1), n.occasions  ))
FT.BPV <- numeric(L)

X2.CR <- numeric(L)
X2.CR.pred <- numeric(L)
X2.BPV <- numeric(L)



for (t in 1:(n.occasions-1)) {
  
  print(t)
  
  for (t2 in 1:(n.occasions)) {
    
    CR <-  sum(Fm[CR.sorted, t, t2] ) 
    id.r <- paste0("fm.post[",1:NCR,", ", t, ", ", t2, "]"  )
    id.R <- paste0("Fm.pred[",CR.sorted,", ", t, ", ", t2, "]"  )
    
    for (l in 1:L) {
      
      lambda  <- sum(fm.post[l, id.r]*rel.F[CR.sorted,t])
      
      FT.CR[l,t, t2] <- sum((sqrt(CR) - sqrt(lambda))^2)
      
      CR.pred <- sum(Fm.pred[l,id.R]) #Rm[CR.sorted[i], t, ]  
      
      FT.CR.pred[l, t, t2] <- sum((sqrt(CR.pred) - sqrt(lambda))^2)
      
    }
    
  }
  
  
}



for (l in 1:L) {
  FT.BPV[l] = sum(FT.CR.pred[l, ,]) > sum(FT.CR[l, ,])
}

mean(FT.BPV)

mean(apply(FT.CR.pred, 1, sum) > apply(FT.CR, 1, sum))

FT.BPV_T = matrix(0, L, n.occasions-1 )

for (l in 1:L) {
  for (t in 1:(n.occasions-1)) {
    FT.BPV_T[l, t] = sum(FT.CR.pred[l,t ,]) > sum(FT.CR[l,t ,]) 
  }
}


plot(colMeans(FT.BPV_T), ylab = "Bayesian p-value", xlab = "Year", main = " ")


# Spatial grouping

FT.CR <- array(0, dim = c(L, NCR, n.occasions))
FT.CR.pred <- array(0, dim = c(L, NCR, n.occasions  ))
FT.BPV <- numeric(L)

X2.CR <- numeric(L)
X2.CR.pred <- numeric(L)
X2.BPV <- numeric(L)



for (i in 1:NCR) {
  
  print(i)
  
  for (t2 in 1:(n.occasions-1)) {
    
    CR <-  (Fm[CR.sorted[i], t2, 1:n.occasions] ) 
    id.r <- paste0("fm.post[",i,", ", t2, ", ", 1:n.occasions, "]"  )
    id.R <- paste0("Fm.pred[",CR.sorted[i],", ", t2, ", ", 1:n.occasions, "]"  )
    
    for (l in 1:L) {
      
      lambda  <- (fm.post[l, id.r]*rel.F[CR.sorted[i],t2])
      
      FT.CR[l,i, t2] <- sum((sqrt(CR) - sqrt(lambda))^2)
      
      CR.pred <- (Fm.pred[l,id.R]) # Posterior predictive 
      
      FT.CR.pred[l, i, t2] <- sum((sqrt(CR.pred) - sqrt(lambda))^2)
      
    }
    
  }
  
  
}


for (l in 1:L) {
  FT.BPV[l] = sum(FT.CR.pred[l, ,]) > sum(FT.CR[l, ,])
}

mean(FT.BPV)

mean(apply(FT.CR.pred, 1, sum) > apply(FT.CR, 1, sum))

FT.BPV_T = diff = matrix(0, L, NCR )

for (l in 1:L) {
  for (i in 1:NCR) {
    FT.BPV_T[l, i] = sum(FT.CR.pred[l,i ,]) > sum(FT.CR[l,i ,])
    diff[l, i] = sum(FT.CR.pred[l,i ,]) - sum(FT.CR[l,i ,])
  }
}


colMeans(FT.BPV_T)
plot(colMeans(FT.BPV_T), xlab = "Sampling location", ylab = "Bayesian p-value", main = " ")
