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

    for (i in 1:(N)){
        for (j in (i):N){
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
    for (i in 1:n_alph){ ll[i] = log_py_spherical(obs, us_init, avec[i]) } 

    a_init = avec[ which.max(ll) ] # take 'most likely' value

    return( a_init )
}

init_us_poincare <- function(obs){
    ## initialise coordinates for spherical coordinates, assumes 2d

    ## get graph distances
    g_obs = graph_from_adjacency_matrix( obs, mode="undirected" )
    ds = distances( g_obs )
    
    ## remove Inf if disconnected
    mx = max( ds[ ds < Inf ] ) # largest finite observed value 
    ds[ ds==Inf ] = mx + 1

    ## get mds coords 
    us_hyd = hydra( ds, curvature=1, alpha=1, equi.adj=0 ) 

    us = cbind( cos(us_hyd$theta), sin(us_hyd$theta) )
    us = sweep( us, 1, us_hyd$r, '*')

    ## add a small amount of noise to each coordinate (some initialised in same place)
    us = us + sample_from_RHN(N, c(0,0), .1)

    return( us )    
}

ident_hyp<-function(Z,n){
    
    
    z_c<-matrix(0,nrow=n,ncol=1)
    
    z_c_new<-matrix(0,nrow=n,ncol=1)
    
    a<-complex(real = Z[1,1], imaginary = Z[1,2])
    
    z_c[1]<-a
    
    z_c[2]<-complex(real = Z[2,1], imaginary = Z[2,2])
    
    b<-sqrt((cosh(h_distance(Z[1,],Z[2,]))-1)/
                (cosh(h_distance(Z[1,],Z[2,]))+1))*
        (Conj(z_c[1])*z_c[2]-1)/(z_c[2]-z_c[1])
    
    #b<-0.5*(1-Conj(z_c[1])*z_c[2])/(z_c[2]-z_c[1])
    
    #if(Re(b*(z_c[2]-a)/(1-Conj(a)*z_c[2]))<0){
    #  b<--b
    #}
    
    for(i in 1:n){
        z_c[i]<-complex(real = Z[i,1], imaginary = Z[i,2])
        z_c_new[i]<-b*(z_c[i]-a)/(1-Conj(a)*z_c[i])
        Z[i,1]<-Re(z_c_new[i])
        Z[i,2]<-Im(z_c_new[i])
    }
    
    if(Z[2,1]<0){
        for(i in 1:n){
            Z[i,1]<--Z[i,1]
        }
    }
    
    if(Z[3,2]<0){
        for(i in 1:n){
            Z[i,2]<--Z[i,2]
        }
    }
    
    #reflection<-Z[3,2]
    
    #if(reflection<0){
    #    Z[3:34,2]<--Z[3:34,2]
    
    #  }
    
    #  print(Z)
    
    return(Z)
    
}

init_alpha_poincare <- function( obs, us_init, n_alph=1000 ){
    ## note: assumes 2-dimensional latent space
    
    ## simple grid search, take the most likely
    ll = rep(0, n_alph)
    avec = seq(-10, 10, length=n_alph)
    for (i in 1:n_alph){ ll[i] = log_py_poincare(obs, us_init, avec[i]) }

    a_init = avec[ which.max(ll) ] # take 'most likely' value

    return( a_init )
}
