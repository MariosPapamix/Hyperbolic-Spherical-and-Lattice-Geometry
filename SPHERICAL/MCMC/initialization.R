library(igraph)
library(smacof) # for 'MDS' in sphere
library(hydra) # for 'MDS' in Poincare disk

spherical_dist <- function( x, y ){ # distance on a sphere
  acos( x %*% y )
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


