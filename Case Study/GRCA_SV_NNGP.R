
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
  
  
  
  inits$wst = nimble::logit(runif(n.occasions-1, 0.6,0.9) )
  inits$wft = log(runif(n.occasions-1,0.3,0.8))
  
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
  
  
  inits$wst = nimble::logit(runif(n.occasions-1, 0.6,0.9) )
  inits$wft = log(runif(n.occasions-1,0.3,0.8))
  
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
# [Warning] dimensions specified are larger than model specification for variable `Rm` is harmless and can be safely ignored
# [Warning] dimensions specified are larger than model specification for variable `Fm` is harmless and can be safely ignored
# The following warnings on the NNGP terms can be safely ignored:
# [Warning] Possible size/dimension mismatch amongst vectors and matrices in BUGS expression: 'w.af[1:697] ~ dmnorm_NNGP(mu = mu[1:697], AD = AD.af[1:697, 1:16], neighbors.id = n.Cid[1:697, 1:15], lower_ = -Inf, upper_ = Inf)'. Ignore this warning if the user-provided distribution has multivariate parameters with distinct sizes or if size of variable differs from sizes of parameters.
# [Warning] Possible size/dimension mismatch amongst vectors and matrices in BUGS expression: 'w.as[1:759] ~ dmnorm_NNGP(mu = mu[1:759], AD = AD.as[1:759, 1:16], neighbors.id = n.id[1:759, 1:15], lower_ = -Inf, upper_ = Inf)'. Ignore this warning if the user-provided distribution has multivariate parameters with distinct sizes or if size of variable differs from sizes of parameters.



# Configure MCMC
FullModelconf = configureMCMC(FullModel)  

# Reset Monitors 
FullModelconf$resetMonitors()

FullModelconf$addMonitors( "rf","w","alpha.p", "eta",  "wst", "wft", "rho.as",  "alpha2.as", "rho.af",  "alpha2.af",
                           "beta.f.tmax", "beta.s.tmax", "beta.f.long", "beta.s.long", "beta.f.lat", "beta.s.lat", "beta.f.tmaxlong", "beta.s.tmaxlong", "beta.f.tmaxlat", "beta.s.tmaxlat"   )

FullModelconf$addMonitors2("C.pred",   "Fm.pred", "Rm.pred", "N",   "sa", "f", "growth.rate", "w.af", "w.as")


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
# These following warning messages are expected and are not a cause for concern:
# Warning] There are X individual pWAIC values that are greater than 0.4. This may indicate that the WAIC estimate is unstable (Vehtari et al., 2017), at least in cases without grouping of data nodes or multivariate data nodes.
# Warning messages:
# 1: Could not calculate the full summary of posterior samples, possibly due to NA or NaN values present in the samples array 
# 2: Could not calculate the full summary of posterior samples, possibly due to NA or NaN values present in the samples array 
# 3: Could not calculate the full summary of posterior samples, possibly due to NA or NaN values present in the samples array 


# MCMC samples for the 1st set of parameters monitored (monitors) 
samples <- FullModelmcmc.out$samples 

# Look at the posterior summaries 
MCMCsummary(samples)

# MCMC samples for the 2nd section of parameters monitored (monitors2)
samples2 <- FullModelmcmc.out$samples2


#=====================================================================
# Code to compute Posterior predictive loss (PPL)
#=====================================================================

# Convert samples2 into one long chain 
samples = mcmc(do.call(rbind, samples2) )

# Count data 
# Posterior predictive samples of counts 
C.pred <- samples[, paste0("C.pred[",1 :n.counts, "]")]

C.pred.mean <- colMeans(C.pred)
C.pred.var <- apply(C.pred, 2, var)

C.G <- sum((C - C.pred.mean)^2)
C.P <- sum(C.pred.var)
C.PPL <- C.G + C.P
C.PPL # PPL for the count data 



# Capture-recapture data 
# Posterior predictive samples of capture-recapture 
Rm.pred <- samples[, grep("Rm.pred", colnames(samples))]

id <- c()
for (i in 1:11) {
  for (j in 1:10) {
    id <- c(id, paste0("Rm.pred[",CR.sorted,", ", j, ", ", i, "]"  ) ) 
  }
}


Rm.pred.mean <- colMeans(Rm.pred[, id], na.rm = TRUE)
Rm.pred.var <- apply(Rm.pred[, id], 2, var)
Rm.c <- c(Rm[CR.sorted, ,]) # convert to scalar to match Rm.pred


Rm.G <- sum((Rm.c - Rm.pred.mean)^2)
Rm.P <- sum(Rm.pred.var)
Rm.PPl <- Rm.G + Rm.P
Rm.PPl


# Fm
Fm.pred <- samples[, grep("Fm.pred", colnames(samples))]

id.f <- c()
for (i in 1:12) {
  for (j in 1:11) {
    id.f <- c(id.f, paste0("Fm.pred[",CR.sorted,", ", j, ", ", i, "]"  ) ) 
  }
}



Fm.pred.mean <- colMeans(Fm.pred[, id.f], na.rm = TRUE)
Fm.pred.var <- apply(Fm.pred[, id.f], 2, var)
Fm.c <- c(Fm[CR.sorted, ,]) # convert to scalar to match Fm.pred


Fm.G <- sum((Fm.c - Fm.pred.mean)^2)
Fm.P <- sum(Fm.pred.var)
Fm.PPl <- Fm.G + Fm.P
Fm.PPl

# Total PPL 
PPL.GP <- C.PPL + Rm.PPl + Fm.PPl
PPL.GP 


# Goodness of fit 
G = C.G + Rm.G + Fm.G
G

# Penalty
P = C.P + Rm.P + Fm.P
P
