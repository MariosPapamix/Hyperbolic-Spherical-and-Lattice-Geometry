library(igraph)
library(smacof) # for 'MDS' in sphere
library(hydra) # for 'MDS' in Poincare disk


hypdistpoinc <- function(x,y){
  acosh( 1 + ( 2 * sum((x-y)^2) ) / ( (1-sum(x^2))*(1-sum(y^2))  ) )
}


log_py_poincare <- function(obs, us, alpha){ 
  
  ll = 0
  N = dim(us)[1]
  
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      eta = alpha - hypdistpoinc( us[i], us[j] )
      ll = ll + obs[i,j]*eta - log(1 + exp(eta))
    }
  }
  return( ll )
}


init_us_poincare <- function(obs){
  ## initialise coordinates for spherical coordinates, assumes 2d
  g_obs = graph_from_adjacency_matrix( obs, mode="undirected" )
  
  ## get graph distances
  g<-dim(Y)[1]
  Dst<-Yr<-Y
  Dst<-Y*(Y==1) + g*(Y==0)
  for(r in 2:(g-1)) {
    Yr<-Yr%*%Y
    Dst<-Dst+(r-g)*( Yr>0 & Dst==g )
  }
  
  #if(infd==T){for(i in 1:g){ for(j in 1:g) { if( Dst[i,j]==g ){ Dst[i,j]<-Inf} }
  #}}
  diag(Dst)<-0
  
  ## get mds coords 
  us_hyd = hydraPlus( Dst, dim=2 ) 
  
  us = cbind( cos(us_hyd$theta), sin(us_hyd$theta) )
  us = sweep( us, 1, us_hyd$r, '*')
  
  
  #print(us)
  
  return( us )    
}

init_alpha_poincare <- function( obs, us_init, n_alph=10 ){
  ## note: assumes 2-dimensional latent space
  
  ## simple grid search, take the most likely
  ll = rep(0, n_alph)
  avec = seq(0, 3, length=n_alph)
  for (i in 1:n_alph){ ll[i] = log_py_poincare(obs, us_init, avec[i]) }
  
  a_init = avec[ which.max(ll) ] # take 'most likely' value
  
  return( a_init )
}


