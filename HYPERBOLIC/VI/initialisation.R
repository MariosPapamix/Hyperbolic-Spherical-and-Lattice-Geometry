library(igraph)
library(smacof) # for 'MDS' in sphere
library(hydra) # for 'MDS' in Poincare disk


init_us_poincare <- function(obs){
    ## initialise coordinates for spherical coordinates, assumes 2d

    ## initialise coordinates for spherical coordinates, assumes 2d
    
    ## get graph distances
    g_obs = graph_from_adjacency_matrix( obs, mode="undirected" )
    #print(g_obs)
    #ds = distances( g_obs )
    
    #n=dim(obs)[1]
    
    ## remove Inf if disconnected
    #mx = max( ds[ ds < Inf ] ) # largest finite observed value 
    #ds[ ds==Inf ] = mx + 1
    
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

init_alpha_poincare <- function( obs, us_init, n_alph=10 ){
    ## note: assumes 2-dimensional latent space
    
    ## simple grid search, take the most likely
    ll = rep(0, n_alph)
    avec = seq(0, 3, length=n_alph)
    for (i in 1:n_alph){ ll[i] = log_py_poincare(obs, us_init, avec[i]) }

    a_init = avec[ which.max(ll) ] # take 'most likely' value

    return( a_init )
}
