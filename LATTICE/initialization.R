library(igraph)
library(smacof) # for 'MDS' in sphere
library(hydra) # for 'MDS' in Poincare disk

spherical_dist <- function( x, y ){ # distance on a sphere
  acos( x %*% y )
}

hypdistpoinc <- function(x,y){
  acosh( 1 + ( 2 * sum((x-y)^2) ) / ( (1-sum(x^2))*(1-sum(y^2))  ) )
}

log_py_spherical <- function(obs, us, alpha){ 
  
  ll = 0
  N = dim(us)[1]
  print(N)
  
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      eta = alpha - spherical_dist( us[i,], us[j,] )
      ll = ll + obs[i,j]*eta - log(1 + exp(eta))
    }
  }
  return( ll )
}

log_py_poincare <- function(obs, us, alpha){ 
  
  ll = 0
  N = dim(us)[1]
  
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      eta = alpha - hypdistpoinc( us[i,], us[j,] )
      ll = ll + obs[i,j]*eta - log(1 + exp(eta))
    }
  }
  return( ll )
}

log_py_lattice <- function(obs, us, alpha){ 
  
  ll = 0
  N = dim(us)[1]
  
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      eta = alpha - l_distance( Z[i,], Z[j,] )
      ll = ll + obs[i,j]*eta - log(1 + exp(eta))
    }
  }
  return( ll )
}

init_us_sphere <- function(obs){
  ## initialise coordinates for spherical coordinates, assumes 3d
  
  ## get graph distances
  g_obs = graph_from_adjacency_matrix( obs, mode="undirected" )
  ds = distances( g_obs )
  
  ## remove Inf if disconnected
  mx = max( ds[ ds < Inf ] ) # largest finite observed value 
  ds[ ds==Inf ] = mx + 1 
  
  ## get mds coords 
  us = smacofSphere( ds, 3, type="ratio")$conf
  us = us / sqrt(rowSums(us^2))  # to ensure lie exactly on sphere
  
  return( us )    
}

init_alpha_sphere <- function( obs, us_init, n_alph=20 ){
  ## note: assumes 3-dimensional
  
  ## simple grid search, take the most likely
  ll = rep(0, n_alph)
  avec = seq(-3, 3, length=n_alph)
  for (i in 1:n_alph){ 
    print(obs)
    print(us_init)
    ll[i] = log_py_spherical(obs, us_init, avec[i])
} 
  
  a_init = avec[ which.max(ll) ] # take 'most likely' value
  
  return( a_init )
}

init_us_poincare <- function(obs){
  ## initialise coordinates for spherical coordinates, assumes 2d
  
  ## get graph distances
  g_obs = graph_from_adjacency_matrix( obs, mode="undirected" )
  ds = distances( g_obs )
  
  n=dim(g_obs)[1]
  
  ## remove Inf if disconnected
  mx = max( ds[ ds < Inf ] ) # largest finite observed value 
  ds[ ds==Inf ] = mx + 1
  
  ## get mds coords 
  us_hyd = hydra( ds, curvature=1, alpha=1, equi.adj=0 ) 
  
  us = cbind( cos(us_hyd$theta), sin(us_hyd$theta) )
  us = sweep( us, 1, us_hyd$r, '*')
  
  ## add a small amount of noise to each coordinate (some initialised in same place)
  us = us + sample_from_RHN(n, c(0,0), .1)
  
  return( us )    
}

init_alpha_poincare <- function( obs, us_init, n_alph=10000 ){
  ## note: assumes 2-dimensional latent space
  
  ## simple grid search, take the most likely
  ll = rep(0, n_alph)
  avec = seq(-100, 100, length=n_alph)
  for (i in 1:n_alph){ ll[i] = log_py_poincare(obs, us_init, avec[i]) }
  
  a_init = avec[ which.max(ll) ] # take 'most likely' value
  
  return( a_init )
}


##Discretize 

init_us_lattice<- function(obs,k=2){
    
  ## initialise coordinates for lattice-integer coordinates, assumes 2d
    
    ## get graph distances
    g_obs = graph_from_adjacency_matrix( obs, mode="undirected" )
    n<-dim(g_obs)[1]
    ds = distances( g_obs )
    
    
    ## remove Inf if disconnected
    mx = max( ds[ ds < Inf ] ) # largest finite observed value 
    ds[ ds==Inf ] = mx + 1 # q: do i want to rescale at all?
    
    ## get mds coords 
    Z = cmdscale( ds, k=2 )
}


init_alpha_lattice <- function( obs, us_init, n_alph=1000 ){
  ## note: assumes 2-dimensional latent space
  
  ## simple grid search, take the most likely
  ll = rep(0, n_alph)
  avec = seq(-100, 100, length=n_alph)
  for (i in 1:n_alph){ ll[i] = log_py_lattice(obs, us_init, avec[i]) }
  
  a_init = avec[ which.max(ll) ] # take 'most likely' value
  
  return( a_init )
}